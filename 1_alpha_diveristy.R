#===============================================================================
# Load libraries and set the environment 
#===============================================================================
load("PRM.Rdata")

load(paste0(PRM$general$dir_Rdata, "/data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Load functions 
source("R/phy_alpha.R")
source("R/da_plot_fun.R")
source("R/lms_with_tests.R")

#===============================================================================
# Variables 
#===============================================================================
DirOut <- PRM$alpha$dir_out

TaxLvl <- PRM$alpha$tax_lvl

TaxNorm <- PRM$alpha$norm

AlphaInds <- PRM$alpha$measures

AlphaIndsTrans <- PRM$alpha$tran_measures

TreatmentCols <- PRM$alpha$treatment_cols

Meta <- DataLs$meta

NcolFacetWrap <- 2

FocusCols <- c("Time", "CID", "Participant", "Group", "Age", "BMI", "Sex")


#===============================================================================
# Analysis schema 
#===============================================================================
SchemaLs <- list()

SchemaLs[["Primary"]] <- 
  list("CID" = "CID2->CID3", 
       "Method" = "lmm", 
       "Formula" = "alpha~Time*Group + Sex + Age + BMI + Base_Delta_alpha + (1|Participant)", 
       "Data_Trans" = c("sqrt", "NONE"))

SchemaLs[["PerTp"]] <- 
  list("CID" = c("Week 0", "Week 8", "Week 32"), 
       "Method" = "Wilcoxon", 
       "Formula" = "alpha~Group", 
       "Data_Trans" = c("NONE"))

# Add data parameters 
for(i in names(SchemaLs)) {
  
  SchemaLs[[i]] <- c(SchemaLs[[i]], 
                      list("Taxa_Lvl" = TaxLvl, 
                           "Count_Norm" = TaxNorm,
                           "Alpha_Index" = AlphaInds))
  
}

# Create folders for output 
FolderCreate <- c("primary", "per_time_point")

for(i in FolderCreate) {
  
  dir.create(paste0(DirOut, "/", i), recursive = TRUE, showWarnings = FALSE)

}

#===============================================================================
# Calculate alpha diversity
#===============================================================================
GridAlpha <- expand.grid("Taxa_Lvl" = TaxLvl, 
                         "Count_Norm" = TaxNorm, 
                         stringsAsFactors = FALSE) 

AlphaTabs <- list()

for(i in 1:nrow(GridAlpha)) {
  
  iTaxaLvl <- GridAlpha[i, "Taxa_Lvl"]
  
  iCountNorm <- GridAlpha[i, "Count_Norm"]
  
  InstPs <- DataLs$PS[[iTaxaLvl]][[iCountNorm]]
  
  InstAlphaTab <- phy_alpha(InstPs, measures = AlphaInds) %>% 
                            bind_cols(Meta)
  
  AlphaTabs[[iTaxaLvl]][[iCountNorm]] <- InstAlphaTab 
  
}


#===============================================================================
# Objects to fill 
#===============================================================================
StatResLs <- list()

AlphaPlots <- list()


#===============================================================================
# Primary 
#===============================================================================
PrmGrid <- expand.grid(SchemaLs[["Primary"]], stringsAsFactors = FALSE)

StatSummDf <- NULL

for(i in 1:nrow(PrmGrid)) {
  
  iCID <- PrmGrid[i, "CID"] %>% 
            str_split("->", simplify = TRUE) %>% 
            as.vector()
  
  iTaxaLvl <- PrmGrid[i, "Taxa_Lvl"]
  
  iCountNorm <- PrmGrid[i, "Count_Norm"]
  
  iIndex <- PrmGrid[i, "Alpha_Index"]
  
  iTrans <- PrmGrid[i, "Data_Trans"]
  
  iMethod <- PrmGrid[i, "Method"]
  
  iFromula <- PrmGrid[i, "Formula"] %>% 
                gsub("alpha", iIndex, .)
  
  # Prepare data 
  DataInst <- AlphaTabs[[iTaxaLvl]][[iCountNorm]] %>% 
                select(all_of(c(FocusCols, iIndex)))
  
  # Transform if needed 
  if(iTrans != "NONE") {
    
    DataInst <- DataInst %>% 
                  mutate(!!iIndex := get(iTrans)(!!sym(iIndex)))
    
  }
  
  # Delta between CID1 and CID2
  DiffSum <- DataInst %>% 
                select(all_of(c("CID", "Participant", iIndex))) %>% 
                filter(CID %in% c("Week 0", "Week 8")) %>% 
                arrange(Participant, CID) %>% 
                summarise(!!paste0("Base_Delta_", iIndex) := 
                            first(!!sym(iIndex)) - last(!!sym(iIndex)), 
                          .by = Participant)
  
  # Combine data
  DataInst <- left_join(DataInst, DiffSum, by = "Participant") %>% 
                filter(CID != "Week 0") %>% 
                droplevels()
  
  StatRes <- lms_with_tests(as.formula(iFromula), 
                             data = DataInst, 
                             type = iMethod)
  
  StatRes$prm <- PrmGrid[i, ]
  
  StatResLs[["Primary"]][[paste0(iIndex, "_", iTrans)]] <- StatRes
  
  StatSummDf <- StatRes$summary$coefficients %>% 
                  as.data.frame() %>% 
                  mutate(across(everything(), function(x){round(x, 4)})) %>% 
                  rownames_to_column(var = "Term") %>% 
                  bind_cols(PrmGrid[i, ]) %>% 
                  add_row() %>% 
                  bind_rows(StatSummDf, .)
            
}

write.csv(StatSummDf, paste0(DirOut, "/primary/stat_primary.csv"), 
          na = "", row.names = FALSE)

StatResLs[["Primary"]][["Summary"]] <- StatSummDf


#-------------------------------------------------------------------------------
# Visualizations
#-------------------------------------------------------------------------------
# How far plot will dodge 
pd <- position_dodge(0.75)

for(i in 1:nrow(GridAlpha)) {
  
  iTaxaLvl <- GridAlpha[i, "Taxa_Lvl"]
  
  iCountNorm <- GridAlpha[i, "Count_Norm"]
  
  iName <- paste0(iTaxaLvl, "_", iCountNorm)
  
  InstAlphaTab <- AlphaTabs[[iTaxaLvl]][[iCountNorm]] %>% 
                    select(all_of(c(AlphaInds, "Group", "Time", "Participant"))) %>% 
                    pivot_longer(cols = all_of(AlphaInds), 
                                 names_to = "Index")
  
  InstSummar <- InstAlphaTab %>% 
                  summarise(Median = median(value), 
                            IQR = IQR(value), .by = c(Group, Time, Index)) %>% 
                  mutate(Yend = Median + IQR/2, 
                         Y = Median - IQR/2)

  # Plots 
  # Rebuild with new layer at the beginning
  AlphaPlots[["Primary"]][[iName]] <- 
    ggplot(InstSummar, aes(x = Time, 
                           y = Median, 
                           colour = Group)) + 
          geom_line(data = InstAlphaTab, 
                              aes(x = Time, 
                                  y = value,
                                  group = Participant, 
                                  colour = Group), 
                              alpha = 0.2, inherit.aes = FALSE) + 
                          geom_line(aes(group = Group), 
                                    position=pd, linewidth = 1) + 
                          geom_point(position=pd, size = 2) +
                          geom_errorbar(aes(ymin = Y, ymax = Yend), 
                                        width=1, position=pd) +
                          facet_wrap(~ Index, scales = "free") + 
                          scale_x_continuous(breaks = sort(unique(InstSummar$Time))) + 
                          theme_bw() + 
                          scale_color_manual(values = AesLs$col$Group) + 
                          xlab("Time (Weeks)") + 
                          ylab("Index value") + 
                          theme(panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank())
  
}

# Save plots 
for(i in names(AlphaPlots$Primary)) {
  
  ggsave(filename = paste0(DirOut, "/primary/", i, ".svg"), 
         plot = AlphaPlots$Primary[[i]], width = 7, height = 2.5)
  
}


#===============================================================================
# Per time point 
#===============================================================================
PrmGrid <- expand.grid(SchemaLs[["PerTp"]], stringsAsFactors = FALSE)

StatSummDf <- NULL

for(i in 1:nrow(PrmGrid)) {
  
  iCID <- PrmGrid[i, "CID"] 
  
  iTaxaLvl <- PrmGrid[i, "Taxa_Lvl"]
  
  iCountNorm <- PrmGrid[i, "Count_Norm"]
  
  iIndex <- PrmGrid[i, "Alpha_Index"]
  
  iFromula <- PrmGrid[i, "Formula"] %>% 
                  gsub("alpha", iIndex, .)
  
  # Prepare data 
  DataInst <- AlphaTabs[[iTaxaLvl]][[iCountNorm]] %>% 
                select(all_of(c(FocusCols, iIndex))) %>% 
                filter(CID == iCID)
  
  StatRes <- wilcox.test(as.formula(iFromula), 
                     data = DataInst, exact = FALSE) 
  
  StatResLs[["PerTp"]][[paste0(iIndex, "_", iTrans)]] <- StatRes
  
  StatSummDf <- c(PrmGrid[i, ], "P_Value" = StatRes$p.value) %>% 
                    bind_rows(StatSummDf, .)
  
}

write.csv(StatSummDf, paste0(DirOut, "/per_time_point/stat_per_tp.csv"))

StatResLs[["PerTp"]][["Summary"]] <- StatSummDf


#-----------------------------------------------------------------------------
# Visualize as violin plots 
#-----------------------------------------------------------------------------
for(i in 1:nrow(GridAlpha)) {
  
  iTaxaLvl <- GridAlpha[i, "Taxa_Lvl"]
  
  iCountNorm <- GridAlpha[i, "Count_Norm"]
  
  iName <- paste0(iTaxaLvl, "_", iCountNorm)
  
  InstAlphaTab <- AlphaTabs[[iTaxaLvl]][[iCountNorm]] %>% 
                      select(all_of(c(AlphaInds, FocusCols))) %>% 
                      pivot_longer(cols = all_of(AlphaInds), 
                                   names_to = "Alpha_Index") %>% 
                      mutate(CID = factor(CID, levels = c("Week 0", 
                                                          "Week 8", 
                                                          "Week 32")))
  
  # Significance df 
  SigDfInst <- InstAlphaTab %>% 
                  summarise(Max_Value = max(value), 
                            Min_Value = min(value),
                            .by = c(CID, Alpha_Index)) %>% 
                  right_join(StatResLs$PerTp$Summary, 
                             by = c("CID", "Alpha_Index")) %>% 
                  mutate(Text = ifelse(round(P_Value, 3) < 0.001, 
                                       "P<0.001", 
                                       paste0("P=", sprintf("%.3f", 
                                                            round(P_Value, 3))))) %>% 
                  mutate(Start = levels(InstAlphaTab$Group)[1], 
                         End = levels(InstAlphaTab$Group)[2], 
                         y = Max_Value + (Max_Value - Min_Value)*0.1, 
                         y_invis = Max_Value + (Max_Value - Min_Value)*0.4, 
                         CID = factor(CID, levels = c("Week 0", 
                                                      "Week 8", 
                                                      "Week 32")))
  
  AlphaPlots[["PerTp"]][[paste0(iName)]] <- 
        ggplot(InstAlphaTab, aes(x = Group, y = value)) + 
                    geom_jitter(aes(colour = Group), 
                                width = 0.2, 
                                height = 0, size = 2, alpha = 0.5) +
                    geom_violin(fill = NA) + 
                    geom_boxplot(fill = NA, width = 0.1, outlier.shape = NA) +
                    geom_signif(data = SigDfInst,
                                aes(xmin = Start,
                                    xmax = End,
                                    annotations = Text,
                                    y_position = y),
                                textsize = 4, vjust = -0.1,
                                manual = TRUE, margin_top = 1) +
                    geom_point(data = SigDfInst,
                               aes(x = End, 
                                   y = y_invis), 
                                 x=NA)  + 
          facet_grid(Alpha_Index~CID, scales = "free") + 
          theme_bw() + 
          scale_color_manual(values = AesLs$col$Group) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
 
}


# Save plots 
for(i in names(AlphaPlots$PerTp)) {
  
  ggsave(filename = paste0(DirOut, "/per_time_point/", i, ".svg"), 
         plot = AlphaPlots$PerTp[[i]], width = 7, height = 4)
  
}


#===============================================================================
# Save and clean up
#===============================================================================
save(list = c("StatResLs", "AlphaPlots"),
     file = paste0(PRM$general$dir_Rdata, "/1_alpha.Rdata"))

rm(list = ls())
gc()
