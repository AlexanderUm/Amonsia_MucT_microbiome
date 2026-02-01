#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

load(paste0(PRM$general$dir_Rdata, "/data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {require(i, character.only = TRUE)}

# Load functions 
source("R/phy_dists_ls.R")
source("R/plot_gg_ordination.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
P <- PRM$beta

if("CLR" %in% P$CountNorm) {
  
  Dists <- c(P$Dists, "Euclidean" = "euclidean")
  
} else {
  
  Dists <- P$Dists
  
}

# Objects to fill 
BetaResLs <- list()

#===============================================================================
# Analysis schema 
#===============================================================================
SchemaLs <- list()

SchemaLs[["Primary"]] <- 
  list("CID" = c("Week 0->Week 8", "Week 8->Week 32"), 
       "Formula" = c("DistMatrix ~ Sex + Age + BMI + Time*Group", 
                     "DistMatrix ~ Sex + Age + BMI + Time"))

SchemaLs[["PerTp"]] <- 
  list("CID" = c("Week 0", "Week 8", "Week 32"), 
       "Formula" = c("DistMatrix ~ Sex + Age + BMI + Group"))

# Add data parameters 
for(i in names(SchemaLs)) {
  
  SchemaLs[[i]] <- c(SchemaLs[[i]], 
                     list("TaxaLvl" = P$TaxaLvl, 
                          "CountNorm" = P$CountNorm,
                          "Dists" = Dists))
  
}

# Create folders for output 
FolderCreate <- c("primary/plots", "per_time_point/plots")

for(i in FolderCreate) {
  
  dir.create(paste0(P$DirOut, "/", i), recursive = TRUE, showWarnings = FALSE)
  
}


#===============================================================================
# Calculate distances 
#===============================================================================
# Object with distances is calculated once and will be used though out the pipeline
Grd <- expand.grid("TaxaLvl" = P$TaxaLvl, 
                   "CountNorm" = P$CountNorm, 
                   stringsAsFactors = FALSE) 

DistsLs <- list()

for(i in 1:nrow(Grd)) {
  
  iTaxaLvl <- Grd$TaxaLvl[i]
  
  iCountNorm <- Grd$CountNorm[i]
  
  # Extract data for all samples
  iPs <- DataLs$PS[[Grd$TaxaLvl[i]]][[Grd$CountNorm[i]]] 
  
  if(iCountNorm == "CLR") {
    
  # Calculate distances 
  DistsLs[[iTaxaLvl]][[iCountNorm]] <- 
        phy_dist_ls(iPs, dists = "euclidean") %>% 
        setNames("Euclidean") 
  
  } else {
    
    DistsLs[[iTaxaLvl]][[iCountNorm]] <- 
          phy_dist_ls(iPs, dists = P$Dists) %>% 
          setNames(names(P$Dists))
    
  }
  
}


#===============================================================================
# Primary 
#===============================================================================
StatResDf <- NULL

PrmGrid <- expand.grid(SchemaLs[["Primary"]], stringsAsFactors = FALSE) %>% 
            filter(!(CountNorm == "CLR" & Dists != "euclidean")) %>% 
            filter(!(CountNorm != "CLR" & Dists == "euclidean")) %>% 
            filter(CID == "Week 0->Week 8" & 
                   Formula == "DistMatrix ~ Sex + Age + BMI + Time"|
                     CID == "Week 8->Week 32" & 
                     Formula == "DistMatrix ~ Sex + Age + BMI + Time*Group") %>% 
            arrange(CID)

for(i in 1:nrow(PrmGrid)) {
  
  iTaxaLvl <- PrmGrid$TaxaLvl[i]
  
  iCountNorm <- PrmGrid$CountNorm[i]
  
  iDist <- PrmGrid$Dists[i]
  
  iCID <- PrmGrid$CID[i] %>% 
            str_split("->", simplify = TRUE) %>% 
            as.vector()
  
  iFormula <- PrmGrid$Formula[i]
  
  # Prepare data 
  iMeta <- DataLs$meta %>% 
              filter(CID %in% iCID)
  
  DistMatrix <- usedist::dist_subset(DistsLs[[iTaxaLvl]][[iCountNorm]][[names(iDist)]], 
                                   rownames(iMeta))
  
  # Statistical testing 
  AdonsResInst <- adonis2(as.formula(iFormula), 
                          data = iMeta, 
                          by = "terms", 
                          permutations = eval(str2lang(P$Permut))) 
  
  StatResDf <- AdonsResInst %>% 
                  as.data.frame() %>% 
                  mutate(across(everything(), function(x){round(x, 4)})) %>% 
                  rownames_to_column(var = "Terms") %>% 
                  bind_cols(PrmGrid[i, ]) %>% 
                  add_row() %>% 
                  bind_rows(StatResDf, .)
  
}

write.csv(StatResDf, file = paste0(P$DirOut, "/primary/", "adonis2.csv"), na = "")

#-------------------------------------------------------------------------------
# Data visualization 
#-------------------------------------------------------------------------------
PrmGridPlot <- expand.grid("CID" = c("Week 0_Week 8_Week 32"),
                           "TaxaLvl" = P$TaxaLvl, 
                           "CountNorm" = P$CountNorm,
                           "Dists" = Dists, 
                           stringsAsFactors = FALSE)

for(i in 1:nrow(PrmGridPlot)) {
  
  iTaxaLvl <- PrmGridPlot$TaxaLvl[i]
  
  iCountNorm <- PrmGridPlot$CountNorm[i]
  
  iDist <- PrmGridPlot$Dists[i]
  
  iCID <- PrmGridPlot$CID[i]
  
  iID <- paste(PrmGridPlot[i, ], collapse = "-")
  
  
 
  # Prepare data 
  iMeta <- DataLs$meta %>% 
                filter(CID %in% str_split(iCID, "_", simplify = TRUE)) %>% 
                arrange(Participant, CID) %>% 
                mutate(Visit = factor(CID, levels = c("Week 0", 
                                                    "Week 8", 
                                                    "Week 32")))
  
  DistMatrix <- usedist::dist_subset(DistsLs[[iTaxaLvl]][[iCountNorm]][[names(iDist)]], 
                                     rownames(iMeta))
  
  if(iCID == "CID1_CID2") {
    
    iColor <- NULL 
    
    } else {
      
      iColor <- "Group"
      
    }
    
  iPlotLs <- plot_gg_ordination(dist_obj = DistMatrix, 
                                meta_data = iMeta, 
                                col_connect_samp = "Participant", 
                                col_color = iColor, 
                                col_shape = "Visit")
  
  iPlot <- iPlotLs$PCoA$p +
              scale_color_manual(values = AesLs$col$Group) + 
              scale_shape_manual(values = AesLs$shape$CID) + 
              ggtitle(label = paste0("  ", names(iDist)))
 
  
  BetaResLs[["Plots"]][["Primary"]][[iTaxaLvl]][[iCountNorm]][[iCID]][[iDist]] <- iPlot
  
  save_plot(filename = paste0(P$DirOut, "/primary/plots/MDS_", iID, ".svg"),
            plot = iPlot,
            base_width = 5.5,
            base_height = 5.5)
  
  save_plot(filename = paste0(P$DirOut, "/primary/plots/MDS_", iID, ".png"),
            plot = iPlot,
            base_width = 5.5,
            base_height = 5.5)
  
}

#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
# Combine and save plots 
PrmGridComb <- PrmGridPlot %>% 
                  select(-Dists) %>% 
                  distinct()

for(i in 1:nrow(PrmGridComb)) {
  
  iTaxaLvl <- PrmGridComb$TaxaLvl[i]
  
  iCountNorm <- PrmGridComb$CountNorm[i]
  
  iCID <- PrmGridComb$CID[i]
  
  iID <- paste(PrmGridComb[i, ], collapse = "-")
  
  # Extract plot to combine 
  iPlotLs <- BetaResLs[["Plots"]][["Primary"]][[iTaxaLvl]][[iCountNorm]][[iCID]]
  
  BetaPlotLs <- lapply(iPlotLs,
                       function(x){x + theme(legend.position = "none")})

  PlotLegend <- get_legend(iPlotLs[[1]])

  PlotGrid0 <- plot_grid(plotlist = BetaPlotLs)

  FullGrid <- plot_grid(PlotGrid0, NULL,
                        PlotLegend,
                        ncol = 3,
                        rel_widths = c(1, 0.0001, 0.15))

  save_plot(filename = paste0(P$DirOut, "/primary/MDS_", iID, ".svg"),
            plot = FullGrid,
            base_width = 10.5,
            base_height = 6.5)
  
  save_plot(filename = paste0(P$DirOut, "/primary/MDS_", iID, ".png"),
            plot = FullGrid,
            base_width = 10.5,
            base_height = 6.5)

}


#===============================================================================
# Per time point 
#===============================================================================
StatResDf <- NULL

PrmGrid <- expand.grid(SchemaLs[["PerTp"]], stringsAsFactors = FALSE) %>% 
              filter(!(CountNorm == "CLR" & Dists != "euclidean")) %>% 
              filter(!(CountNorm != "CLR" & Dists == "euclidean"))

for(i in 1:nrow(PrmGrid)) {
  
  iTaxaLvl <- PrmGrid$TaxaLvl[i]
  
  iCountNorm <- PrmGrid$CountNorm[i]
  
  iDist <- PrmGrid$Dists[i]
  
  iCID <- PrmGrid$CID[i] %>% 
              str_split("->", simplify = TRUE) %>% 
              as.vector()
  
  iFormula <- PrmGrid$Formula[i]
  
  # Prepare data 
  iMeta <- DataLs$meta %>% 
              filter(CID %in% iCID) %>% 
              mutate(CID = factor(CID, levels = c("Week 0", 
                                                  "Week 8", 
                                                  "Week 32")))
  
  DistMatrix <- usedist::dist_subset(DistsLs[[iTaxaLvl]][[iCountNorm]][[names(iDist)]], 
                                     rownames(iMeta))
  
  # Statistical testing 
  AdonsResInst <- adonis2(as.formula(iFormula), 
                          data = iMeta, 
                          by = "terms", 
                          permutations = eval(str2lang(P$Permut))) 
  
  StatResDf <- AdonsResInst %>% 
                    as.data.frame() %>% 
                    mutate(across(everything(), function(x){round(x, 4)})) %>% 
                    rownames_to_column(var = "Terms") %>% 
                    bind_cols(PrmGrid[i, ]) %>% 
                    add_row() %>% 
                    bind_rows(StatResDf, .)
  
}

write.csv(StatResDf, file = paste0(P$DirOut, "/per_time_point/", "adonis2.csv"), na = "")


#-------------------------------------------------------------------------------
# Data visualization 
#-------------------------------------------------------------------------------
PrmGridPlot <- PrmGrid %>% 
                  select(-c(Formula)) %>% 
                  distinct()

for(i in 1:nrow(PrmGridPlot)) {
  
  iTaxaLvl <- PrmGridPlot$TaxaLvl[i]
  
  iCountNorm <- PrmGridPlot$CountNorm[i]
  
  iDist <- PrmGridPlot$Dists[i]
  
  iCID <- PrmGridPlot$CID[i]
  
  iID <- paste(PrmGridPlot[i, ], collapse = "-")
  
  
  # Prepare data 
  iMeta <- DataLs$meta %>% 
              filter(CID %in% iCID)
  
  DistMatrix <- usedist::dist_subset(DistsLs[[iTaxaLvl]][[iCountNorm]][[names(iDist)]], 
                                     rownames(iMeta))
  
  iPlotLs <- plot_gg_ordination(dist_obj = DistMatrix, 
                                meta_data = iMeta, 
                                col_color = "Group", 
                                col_shape = "CID")
  
  iPlot <- iPlotLs$PCoA$p +
              scale_color_manual(values = AesLs$col$Group) + 
              scale_shape_manual(values = AesLs$shape$CID) + 
              ggtitle(label = paste0("  ", iCID, ": ", names(iDist)))
  
  BetaResLs[["Plots"]][["PerTp"]][[iTaxaLvl]][[iCountNorm]][[paste0(iDist, "-", iCID)]] <- iPlot
  
  save_plot(filename = paste0(P$DirOut, "/per_time_point/plots/MDS_", iID, ".svg"),
            plot = iPlot,
            base_width = 5.5,
            base_height = 5.5)
  
}

#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
# Combine and save plots 
PrmGridComb <- PrmGridPlot %>% 
                  select(-c(Dists, CID)) %>% 
                  distinct()

for(i in 1:nrow(PrmGridComb)) {
  
  iTaxaLvl <- PrmGridComb$TaxaLvl[i]
  
  iCountNorm <- PrmGridComb$CountNorm[i]
  
  iID <- paste(PrmGridComb[i, ], collapse = "-")
  
  # Extract plot to combine 
  iPlotLs <- BetaResLs[["Plots"]][["PerTp"]][[iTaxaLvl]][[iCountNorm]]
  
  BetaPlotLs <- lapply(iPlotLs,
                       function(x){x + theme(legend.position = "none")})
  
  PlotGrid0 <- plot_grid(plotlist = BetaPlotLs, ncol = 3)
  
  FullGrid <- plot_grid(PlotGrid0, NULL,
                        PlotLegend,
                        ncol = 3,
                        rel_widths = c(1, 0.0001, 0.15))
  
  save_plot(filename = paste0(P$DirOut, "/per_time_point/MDS_", iID, ".svg"),
            plot = FullGrid,
            base_width = 12.5,
            base_height = 15)
  
}


#-------------------------------------------------------------------------------
# Save and clean up 
#-------------------------------------------------------------------------------
save(list = c("BetaResLs"), 
     file = paste0(PRM$general$dir_Rdata, "/2_beta.Rdata"))

rm(list = ls())
gc()
