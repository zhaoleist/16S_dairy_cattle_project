							##########################
							#    loading packages    #
							##########################
library(dplyr)
library(reshape2)
library(vegan)
library(car) # levene test
library(ggplot2)
library(ggrepel)
library(ggalt) # encicle
library(RColorBrewer)
library(superheat)
library(qgraph) # interaction network
library(ggpubr)
# library(VennDiagram)
# library(gridExtra) # layout of ggplot objects
# library(pheatmap)
# library(corrplot)

							###################
							#    set up wd    #
							###################
if (Sys.info()["sysname"] =='Windows') {
  working_dir <- 'D:\\MyZone\\Dropbox\\pool_80_for_paper\\16s_dairyCattle\\code_in_r'
} else {
  working_dir <- '/Users/lzhao13/Dropbox/pool_80_for_paper/16s_dairyCattle/'
}
setwd(working_dir)

							#############################
							#    constant & variables   #
							#############################
N <- 90
MIN_DEPTH <- 12637

taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
col_index_list <-list(col_index_HC = 1:9,    # 9
                      col_index_PWH = 10:15,    # 6
                      col_index_BH = 16:23,    # 8
                      col_index_SCU = 24:31,    # 8
                      col_index_FH = 32:39,    # 8
                      col_index_MLLH = 40:47,   # 8
                      col_index_FC = 48:55,    # 8
                      col_index_MLC = 56:63,    # 8
                      col_index_PLLC = 64:71,    # 8
                      col_index_FOC = 72:74,    # 3
                      col_index_CUC = 75:82,    # 8
                      col_index_Hos = 83:90)    # 8
stage <- c(rep("01-HC", 9),
           rep("02-PWH", 6),
           rep("03-BH", 8),
           rep("04-SCU", 8),
           rep("05-FH", 8),
           rep("06-MLLH", 8),
           rep("07-FC", 8),
           rep("08-MLC", 8),
           rep("09-PLLC", 8),
           rep("10-FOC", 3),
           rep("11-CUC", 8),
           rep("12-Hos", 8))
production_tracks <- c("HC.001", "HC.002", "HC.003", "HC.004", "HC.005", "HC.006", "HC.007", "HC.008", "HC.009",
                       "PWH.001", "PWH.002", "PWH.003", "PWH.004", "PWH.005", "PWH.006",
                       "BH.001", "BH.002", "BH.003", "BH.004", "BH.005", "BH.006", "BH.007", "BH.008",
                       "SCU.001", "SCU.002", "SCU.003", "SCU.004", "SCU.005", "SCU.006", "SCU.007", "SCU.008",
                       "FH.001", "FH.002", "FH.003", "FH.004", "FH.005", "FH.006", "FH.007", "FH.008",
                       "MLLH.001","MLLH.002", "MLLH.003", "MLLH.004", "MLLH.005", "MLLH.006", "MLLH.007", "MLLH.008",
                       "FC.001", "FC.002","FC.003", "FC.004", "FC.005", "FC.006", "FC.007", "FC.008",
                       "MLC.001", "MLC.002", "MLC.003", "MLC.004", "MLC.005", "MLC.006", "MLC.007", "MLC.008",
                       "PLLC.001","PLLC.002","PLLC.003", "PLLC.004","PLLC.005","PLLC.006","PLLC.007","PLLC.008",
                       "FOC.001", "FOC.002", "FOC.003",
                       "CUC.001", "CUC.002", "CUC.003", "CUC.004", "CUC.005","CUC.006", "CUC.007", "CUC.008",
                       "Hos.001", "Hos.002",  "Hos.003",  "Hos.004",  "Hos.005",  "Hos.006", "Hos.007", "Hos.008",
                       "Kingdom", "Phylum","Class","Order","Family","Genus")    # correct order

							################
							#    palette   #
							################
my_colors <- c("steelblue2","blue3","seagreen3","green4","dimgray","coral4",    # 12 carefully picked colors
              "pink2","gold","darkorange2","deeppink","purple","red1")

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
