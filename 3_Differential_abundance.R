#===============================================================================
# Load libraries and set the environment 
#===============================================================================
load("PRM.Rdata")

load(paste0(PRM$general$dir_Rdata, "/data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Functions
source("R/prev_filt_phy.R")

#===============================================================================
# Analysis schema 
#===============================================================================
P <- PRM$da

SchemaLs <- list()

SchemaLs[["Shift"]] <- 
  list("CID" = c("Week 0->Week 8", "Week 8->Week 32"), 
       "Method" = "linda", 
       "Formula" = c("~Time*Group + Sex + Age + BMI + (1|Participant)", 
                     "~Time + Sex + Age + BMI + (1|Participant)"), 
       "CountNorm" = "Raw")

SchemaLs[["PerTp"]] <- 
  list("CID" = c("Week 0", "Week 8", "Week 32"), 
       "Method" = "Wilcoxon", 
       "Formula" = "~Group", 
       "CountNorm" = "TSS")

# Add data parameters 
for(i in names(SchemaLs)) {
  
  SchemaLs[[i]] <- c(SchemaLs[[i]], 
                     list("TaxaLvl" = P$TaxaLvl))
  
}

# Create folders for output 
FolderCreate <- c("per_time_point/tabs", "shift/tabs")

for(i in FolderCreate) {
  
  dir.create(paste0(P$DirOut, "/", i), recursive = TRUE, showWarnings = FALSE)
  
}

#===============================================================================
# Objects to fill 
#===============================================================================
StatResLs <- list()

StatResSigLs <- list()

#===============================================================================
# Extract OTU tables that will be used for DA 
#===============================================================================
Norms <- unlist(SchemaLs) %>% 
          .[grep("CountNorm", names(.))]

Grd <- expand.grid("TaxaLvl" = P$TaxaLvl, 
                   "CountNorm" = Norms, 
                   stringsAsFactors = FALSE)

OtuTabs <- NULL

for(i in 1:nrow(Grd)) {
  
  iTaxaLvl <- Grd[i, "TaxaLvl"]
  
  iCountNorm <- Grd[i, "CountNorm"]
  
  InstPs <- DataLs$PS[[iTaxaLvl]][[iCountNorm]] 
  
  if(!taxa_are_rows(InstPs)){
    
    InstPs <- phyloseq::t(InstPs)
    
  }
  
  iOtuTab <-  InstPs %>% 
                otu_table() %>% 
                as.matrix() %>% 
                as.data.frame()
  
  # Remove low prevalent taxa
  iMeta <- DataLs$meta %>% 
            filter(CID %in% P$PrevByCDSs)
  
  iOtuTabFilt <- iOtuTab[, rownames(iMeta)] 
  
  iTaxa <- rowSums(iOtuTabFilt != 0)/ncol(iOtuTabFilt)
  
  iTaxa <- iTaxa[iTaxa >= P$MinPrev]
  
  # Collect filtered OTU table 
  OtuTabs[[iTaxaLvl]][[iCountNorm]] <- iOtuTab[names(iTaxa), ]
  
}


#===============================================================================
# Check trend in shift 
#===============================================================================
PrmGrid <- expand.grid(SchemaLs[["Shift"]], stringsAsFactors = FALSE) %>% 
            filter(CID == "Week 0->Week 8" & 
                     Formula == "~Time + Sex + Age + BMI + (1|Participant)" | 
                   CID == "Week 8->Week 32" & 
                     Formula == "~Time*Group + Sex + Age + BMI + (1|Participant)")

StatSigDf <- NULL

for(i in 1:nrow(PrmGrid)) {
  
  iCidId <- PrmGrid[i, "CID"]
  
  iCID <- PrmGrid[i, "CID"] %>% 
              str_split("->", simplify = TRUE) %>% 
              as.vector()
  
  iTaxaLvl <- PrmGrid[i, "TaxaLvl"]
  
  iCountNorm <- PrmGrid[i, "CountNorm"]
  
  iFormula <- PrmGrid[i, "Formula"]
  
  iID <- paste0(P$DirOut, "/shift/tabs/", iTaxaLvl, "/", 
                paste0(iCID, collapse = "_"))
  
  dir.create(path = iID, recursive = TRUE)
  
  # Prepare data 
  iMeta <- DataLs$meta %>% 
              filter(CID %in% iCID)
  
  iOtuTab <- OtuTabs[[iTaxaLvl]][[iCountNorm]][, rownames(iMeta)]
  
  # LinDa test 
  iLinDaRes <- linda(feature.dat = iOtuTab, 
                      meta.dat = iMeta, 
                      formula = iFormula, 
                      feature.dat.type = "count", 
                      pseudo.cnt = 1)
  
  StatResLs[["Shift"]][[iTaxaLvl]][[iCountNorm]][[iCidId]] <- iLinDaRes
  
  # Select significant results 
  for(j in names(iLinDaRes$output)) {
    
    
    write.csv(iLinDaRes$output[[j]], 
              file = paste0(iID, "/", 
                            iCountNorm, "_", 
                            gsub(":", "-", j), 
                            ".csv"))
    
    iSigRes <- iLinDaRes$output[[j]] %>% 
                    filter(padj <= 0.1)
    
    if(nrow(iSigRes) > 0) {
      
      StatSigDf <- iSigRes %>% 
                    bind_cols(PrmGrid[i, ], "Term" = j)
      
    }
    
  }
    
}

StatResSigLs[["Shift"]] <- StatSigDf


#===============================================================================
# Per time point 
#===============================================================================
PrmGrid <- expand.grid(SchemaLs[["PerTp"]], stringsAsFactors = FALSE)

StatResPerTpDf <- NULL

for(i in 1:nrow(PrmGrid)) {
  
  iCidId <- PrmGrid[i, "CID"]
  
  iCID <- PrmGrid[i, "CID"] %>% 
              str_split("->", simplify = TRUE) %>% 
              as.vector()
  
  iTaxaLvl <- PrmGrid[i, "TaxaLvl"]
  
  iCountNorm <- PrmGrid[i, "CountNorm"]
  
  iFormula <- PrmGrid[i, "Formula"]
  
  iID <- paste0(P$DirOut, "/per_time_point/tabs/", iTaxaLvl, "/")
  
  dir.create(iID, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data 
  iMeta <- DataLs$meta %>% 
                filter(CID %in% iCID)
  
  iOtuTab <- OtuTabs[[iTaxaLvl]][[iCountNorm]][, rownames(iMeta)] %>% 
                t() %>% 
                as.data.frame() 
  
  iCombData <- bind_cols(iMeta, iOtuTab)  
  
  iStatResComb <- NULL
  
  for(j in colnames(iOtuTab)) {
    
    jFormula <- paste0(j, iFormula)
    
    jRes <- wilcox.test(as.formula(jFormula), 
                           data = iCombData, exact = FALSE) 
    
    iStatResComb <- c("Taxa" = j, "P_Value" = jRes$p.value, PrmGrid[i, ]) %>% 
                      bind_rows(iStatResComb, .)
    
  }
  
  iStatResComb$Q_value <- p.adjust(iStatResComb$P_Value, method = "fdr")
  
  write.csv(iStatResComb, paste0(iID, "/", iCountNorm, "_", iCID, ".csv"))
  
  StatResPerTpDf <- bind_rows(StatResPerTpDf, iStatResComb)
  
}

# Clean environment 
rm(list = ls())
gc()
