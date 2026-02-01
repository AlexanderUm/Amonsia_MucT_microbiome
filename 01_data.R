#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

set.seed(PRM$general$seed)

for(i in PRM$general$libs) {
  
  require(i, character.only = TRUE)
  
}


#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_shorten_tax_names.R")
source("R/phy_norm_count.R")
source("R/phy_filter_taxa_by_name.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
DietCols <- PRM$general$diet_col

DirRdata <- PRM$general$dir_Rdata

DirOut <- PRM$general$dir_out

# Data preparation
qPath <- PRM$data$q_path

MetaPath <- PRM$data$m_path

TaxLvls <- PRM$data$tax_lvls

MinReads <- PRM$data$min_read_tax

MinReadsSample <- PRM$data$min_read_sample

FiltTaxaPrm <- PRM$data$tax_filt

FiltSamples <- PRM$data$samples_to_remove

# Cols 
ColsToFactors <- PRM$data$col_as_factor

ColsToRename <- PRM$data$col_to_rename

# Create directories 
sapply(c(DirRdata, DirOut), 
       function(x) {dir.create(x, showWarnings = FALSE, recursive = TRUE)})


# Object to fill up 
DataLs <- list()

ReportData <-list()


################################################################################
# Read data/ Create phyloseq
#-------------------------------------------------------------------------------
ps <- qza_to_phyloseq(features = paste0(qPath, "asv_table.qza"), 
                      tree = paste0(qPath, "rooted-tree.qza"), 
                      taxonomy = paste0(qPath, "taxonomy_07.qza"))

# Read in metadata 
Meta <- read.csv(MetaPath) %>% 
          mutate(RownNamesID = SeqID) %>% 
          column_to_rownames("RownNamesID")

# Adjust columns type
Meta <- Meta %>% 
          rename(all_of(ColsToRename)) %>% 
          mutate(Sex = case_match(Sex, 1 ~ "M", 2 ~ "F"), 
                 Time = case_match(CID, "CID1" ~ 0, "CID2" ~ 8, "CID3" ~ 32), 
                 CID = case_match(CID, 
                                  "CID1" ~ "Week 0", 
                                  "CID2" ~ "Week 8", 
                                  "CID3" ~ "Week 32")) %>% 
          filter(!Participant %in% FiltSamples) %>% 
          filter(n() == 3, .by = Participant) %>% 
          mutate(across(all_of(c(ColsToFactors, "Sex")), 
                        function(x){factor(x, levels = sort(unique(x)))}))

# Combine with phyloseq
OverSamp <- intersect(rownames(Meta), sample_names(ps))

ps <- prune_samples(OverSamp, ps)

Meta <- Meta[sample_names(ps), ]

sample_data(ps) <- Meta

DataLs[["meta"]] <- Meta

ReportData[["ps"]][["original"]] <- ps

# Write out metadata 
write.csv(Meta, file = paste0(DirOut, "/metadata_used_for_analysis.csv"))

#-------------------------------------------------------------------------------
# Taxa - filter taxa with with less than X reads in total   
#-------------------------------------------------------------------------------
ps <- prune_taxa(taxa_sums(ps) >= MinReads, ps)

#-------------------------------------------------------------------------------
# Taxa - Remove taxa that can be considerate contaminants or not belonging to 
#           microbial community AND keep only relevant taxa
#-------------------------------------------------------------------------------
ps <- phy_filter_taxa_by_name(phy = ps, 
                              keep_ls = FiltTaxaPrm$keep, 
                              kick_ls = FiltTaxaPrm$kick)


################################################################################
# Prepare list with data 
# Glom to higher taxonomic level ->
# Normalize count ->
# Make strata
#-------------------------------------------------------------------------------
for(i in TaxLvls)  {
  
  if(i == "ASV") { 
    
    InstPs <- ps
    
    taxa_names(InstPs) <- phy_shorten_tax_names(InstPs) %>% 
                              as.data.frame() %>% 
                              setNames("feature") %>% 
                              mutate(feature2 = if(n() > 1) {
                                paste0(feature, "__asv", row_number())} else {
                                  paste0(feature, "__asv")}, .by = "feature") %>% 
                              pull(feature2)
    
  } else { 
    
    InstPs <- tax_glom(ps, i)
    
    }
  
  # Adjust taxa names
  taxa_names(InstPs) <- phy_shorten_tax_names(InstPs) %>% 
                                        as.data.frame() %>% 
                                        setNames("feature") %>% 
                                        mutate(feature2 = if(n() > 1) {
                                          paste0(feature, "__", 
                                                 tolower(str_sub(i, 1, 1)), 
                                                 row_number())} else {
                                                   feature}, .by = "feature") %>% 
                                        pull(feature2)
  
  InstPsRare <- rarefy_even_depth(InstPs, rngseed = PRM$general$seed)
  
  DataLs[["PS"]][[i]] <- list("Raw" = InstPs, 
                              "Rare" = phy_norm_count(InstPs, 
                                                      norm_type = "Rare", 
                                                      seed = 34957943), 
                              "CSS" = phy_norm_count(InstPs, 
                                                     norm_type = "CSS_log2"), 
                              "CLR" = phy_norm_count(InstPs, 
                                                     norm_type = "CLR"),
                              "TSS" = phy_norm_count(InstPs, 
                                                     norm_type = "TSS"))

}


#===============================================================================
# Write out taxonomy table (will be used by collaborators)
TaxTab <- DataLs$PS$ASV$Raw %>% 
          tax_table() %>% 
          as.data.frame()

write.csv(TaxTab, file = paste0(DirOut, "/taxonomy_table.csv"))



################################################################################
# Aesthetics 
#-------------------------------------------------------------------------------
AesLs <- list()

ColVec <-  c("#999999", "#7FC97F", "#BEAED4", "#FDC086", "#386CB0", # 5
             "#BF5B17", "#1B9E77", "#7570B3", "#66A61E", "#E6AB02", # 10
             "#A6761D", "#A6CEE3", "#33A02C", "#FB9A99", "#FDBF6F", # 15
             "#FF7F00", "#6A3D9A", "#FFFF99", "#B15928", "#CCEBC5", # 20
             "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", # 25
             "#B3E2CD", "#FDCDAC", "#CBD5E8", "#E6F5C9", "#FFF2AE", # 30
             "#F1E2CC", "#CCCCCC", "#E41A1C", "#984EA3", "#FF7F00", # 35
             "#FFFF33", "#A65628", "#66C2A5", "#FC8D62", "#8DA0CB", # 40
             "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#8DD3C7", # 45
             "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", # 50
             "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F") # 55

ShapeVec <- c(16, 17, 18, 15, 8, 1:7)

AesLs[["col"]][["Group"]] <- ColVec[c(5, 1)] %>% 
                                    setNames(levels(Meta[["Group"]]))

AesLs[["col"]][["CID"]] <- ColVec[c(7, 8, 10)] %>% 
                                    setNames(levels(Meta[["CID"]]))

AesLs[["shape"]][["CID"]] <- ShapeVec[c(1, 2, 3)] %>% 
                                      setNames(levels(Meta[["CID"]]))


#-------------------------------------------------------------------------------
save(list = c("DataLs", "AesLs", "ReportData"), 
     file = paste0(DirRdata, "/data.Rdata"))

# Clean environment 
rm(list = ls())
gc()
