
# source("01_setup.R")

#####################################
#    Read DADA2 output & process    #
#####################################
# exclude Archaea, Eukaryota, & NA
my_table <- read.table("myTable.txt", header=TRUE, sep="\t") # my_table is the raw output OTU table from DADA2 --> 7,129 x 96
rows_to_rmv <- c(which(is.na(my_table["Kingdom"])),which(my_table["Kingdom"]=="Archaea"),which(my_table["Kingdom"]=="Eukaryota"))
my_table <- my_table[-rows_to_rmv,] # 6696 x 96

# change "None" to "Unassigned"
for (i in 92:96){
  levels(my_table[,i]) <- c(levels(my_table[,i]),"None")
  my_table[is.na(my_table[,i]),i] <- "None"
  my_table[,i] <- as.character(my_table[,i])
  my_table[my_table[,i]=="None",i] <- "Unassigned"
  my_table[,i] <- factor(my_table[,i])
}

#################################
#    Read rarefied OTU table    #
#################################
my_table_2 <- read.table("my_otu_table_even_10_12637.txt", header=TRUE, sep="\t") # 4,681 x 92
# expand 'taxonomy' column (separated by ';') to six cloumns
my_table_2[,taxa_levels] <- do.call(rbind, strsplit(as.character(my_table_2$taxonomy), ';'))
my_table_2[, c("OTU_ID", "taxonomy")] <- NULL # 4,681 x 96
my_table_2 <- my_table_2[,production_tracks] # change the col order to management unit order
  
apply(my_table_2[,1:90], 2, function(x) sum(x != 0))/MIN_DEPTH     # non zero ASVs ratio for each sample

###############################################################
#    Good's Coverage --> ratio of non-singleton abundances    #
###############################################################
GoodsCov_bf <- vector() # 87 samples > 0.97
GoodsCov_af <- vector() # 89 samples > 0.95
for (i in 1:ncol(my_table_2[,1:90])){
  GoodsCov_af[i] <- 1 - sum(my_table_2[which(my_table_2[,i]==1),i])/sum(my_table_2[,i])
}
for (i in 1:ncol(my_table[,1:90])){
  GoodsCov_bf[i] <- 1 - sum(my_table[which(my_table[,i]==1),i])/sum(my_table[,i])
}

###########################
#    Rarefaction Curve    #
###########################
my_table <- my_table[,production_tracks]
colors_to_use <- rep(my_colors, lengths(col_index_list))

# without background lines
png(file="Fig.S1_rarecurve_wo_sample.png", width=1200, height=800)
par(mar = c(5, 5, 4, 4))
rarecurve(t(my_table[,1:90]), step=500, col=colors_to_use, cex=0.5, lwd=3, label=FALSE,
          xlab="Number of Sequencing Reads", ylab="Number of Observed OTUs",cex.lab=2, cex.axis=2) # from vegan
dev.off()

###########################################################
#    Naively change "Unassigned" to highest known taxa    #
###########################################################
my_table_2[,taxa_levels] <- lapply(my_table_2[,taxa_levels], as.character)

# genus
for (i in 1:nrow(my_table_2)){
  if (my_table_2[i,"Genus"] == "Unassigned")
    if (my_table_2[i, "Family"] != "Unassigned")
    {my_table_2[i, "Genus"] <- paste0("[f]",my_table_2[i, "Family"])
    } else if (my_table_2[i, "Order"] != "Unassigned")
    {my_table_2[i, "Genus"] <- paste0("[o]",my_table_2[i, "Order"])
    } else if (my_table_2[i, "Class"] != "Unassigned")
    {my_table_2[i, "Genus"] <- paste0("[c]",my_table_2[i, "Class"])
    } else if (my_table_2[i, "Phylum"] != "Unassigned")
    {my_table_2[i, "Genus"] <- paste0("[p]",my_table_2[i, "Phylum"])
    } else if (my_table_2[i, "Kingdom"] != "Unassigned")
      my_table_2[i, "Genus"] <- paste0("[k]",my_table_2[i, "Kingdom"])
}
# family
for (i in 1:nrow(my_table_2)){
  if (my_table_2[i, "Family"] == "Unassigned")
    if (my_table_2[i, "Order"] != "Unassigned")
    {my_table_2[i, "Family"] <- paste0("[o]",my_table_2[i, "Order"])
    } else if (my_table_2[i, "Class"] != "Unassigned")
    {my_table_2[i, "Family"] <- paste0("[c]",my_table_2[i, "Class"])
    } else if (my_table_2[i, "Phylum"] != "Unassigned")
    {my_table_2[i, "Family"] <- paste0("[p]",my_table_2[i, "Phylum"])
    } else if (my_table_2[i, "Kingdom"] != "Unassigned")
      my_table_2[i, "Family"] <- paste0("[k]",my_table_2[i, "Kingdom"])
}    
# order
for (i in 1:nrow(my_table_2)){
  if (my_table_2[i, "Order"] == "Unassigned")
    if (my_table_2[i, "Class"] != "Unassigned")
    {my_table_2[i, "Order"] <- paste0("[c]",my_table_2[i, "Class"])
    } else if (my_table_2[i, "Phylum"] != "Unassigned")
    {my_table_2[i, "Order"] <- paste0("[p]",my_table_2[i, "Phylum"])
    } else if (my_table_2[i, "Kingdom"] != "Unassigned")
      my_table_2[i, "Order"] <- paste0("[k]",my_table_2[i, "Kingdom"])
}    
# class
for (i in 1:nrow(my_table_2)){
  if (my_table_2[i, "Class"] == "Unassigned")
    if (my_table_2[i, "Phylum"] != "Unassigned")
    {my_table_2[i, "Class"] <- paste0("[p]",my_table_2[i, "Phylum"])
    } else if (my_table_2[i, "Kingdom"] != "Unassigned")
      my_table_2[i, "Class"] <- paste0("[k]",my_table_2[i, "Kingdom"])
}    
# phylum
for (i in 1:nrow(my_table_2)){
  if (my_table_2[i, "Phylum"] == "Unassigned")
    if (my_table_2[i, "Kingdom"] != "Unassigned")
      my_table_2[i, "Phylum"] <- paste0("[k]",my_table_2[i, "Kingdom"])
}


