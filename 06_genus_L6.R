#############################################
#    Plot no. of ASVs assigned to genera    #
#############################################
df_no_of_ASV_assigned_to_genus <- data.frame(genus=sort(table(my_table_2[,96]), decreasing = TRUE))
df_no_of_ASV_assigned_to_genus <- df_no_of_ASV_assigned_to_genus[1:20,]
g_no_of_ASV_assigned_to_genus <- ggplot() +
  geom_bar(data=df_no_of_ASV_assigned_to_genus, aes(x=genus.Var1, y=genus.Freq),fill="white", color="black", size=3, stat="identity") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 54, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 54, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Genera (partial)") + 
  ylab("No. of OTUs") +
  annotate("text", x= 1:nrow(df_no_of_ASV_assigned_to_genus),
           y=c(250,rep(max(df_no_of_ASV_assigned_to_genus$genus.Freq) * 0.5,5), rep(250,14)),
           # y=c(max(df_no_of_ASV_assigned_to_genus$genus.Freq) * 0.5, 505, 475, rep(411.6,17)),
           label=levels(df_no_of_ASV_assigned_to_genus$genus.Var1)[1:20], 
           angle=90, size=22, fontface="bold", hjust=0, color="black") +
  scale_y_continuous(limits=c(0, max(df_no_of_ASV_assigned_to_genus$genus.Freq)+200), breaks=seq(0,max(df_no_of_ASV_assigned_to_genus$genus.Freq), 250))
ggsave("Fig.2E.tiff", plot = g_no_of_ASV_assigned_to_genus, width=32, height=24, dpi=300, compression="lzw")

###########################################
#    Relative abundance at genus level    #
###########################################
genus_uniq <- unique(my_table_2["Genus"])[[1]]
myL6 <- data.frame(matrix(nrow = length(genus_uniq), ncol=91))
for (i in 1:length(genus_uniq)){
  index <- which(my_table_2["Genus"]==as.character(genus_uniq[i]))
  myL6[i,1:N] <- apply(my_table_2[index,1:N],2,sum)
  myL6[i,91] <- as.character(genus_uniq[i])
  names(myL6) <- c(names(my_table_2[1:N]),"Genus")
}

# No. of samples containing each of the genera
myL6_temp <- as.data.frame(ifelse(myL6[,1:90] == 0, 0, 1))
df_fig_3h <- data.frame(value=apply(myL6_temp,1,sum),
                        genus=myL6$Genus)
df_fig_3h$genus <- factor(df_fig_3h$genus, levels = levels(df_no_of_ASV_assigned_to_genus$genus.Var1))
df_fig_3h <- df_fig_3h[1:20,]

g_fig_temp <- ggplot() +
  geom_bar(data=df_fig_3h, aes(x=genus, y=value), fill="white", col="black", stat="identity") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 28, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 28, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Genus (patial)") + 
  ylab("No. of samples containing the genera") +
  annotate("text", x= 1:nrow(df_fig_3h), 
           y=max(df_fig_3h$value) * 0.445, 
           label=levels(df_fig_3h$genus)[1:20], angle=90, size=10, fontface="bold") +
  scale_y_continuous(limits=c(0, max(df_fig_3h$value) + 5), breaks=seq(0,max(df_fig_3h$value), 10))
# ggsave("Fig.999_genus_in_no_of_samples.tiff", plot = g_fig_temp, width=16, height=12, dpi=300, compression="lzw")

# Relative abundance
index_shared_by_all_L6 <- which(apply(myL6[,1:90], 1, function(x) all(x != 0))  == TRUE)    # which families shared by all sample?
myL6[index_shared_by_all_L6, "Genus"]
sum(apply(myL6[index_shared_by_all_L6, 1:90], 1, sum))/sum(apply(myL6[, 1:90], 1, sum))
sum_L6 <- apply(myL6[1:N],2,sum)

# taxon ratio each sample
myL6_percentage <- data.frame(matrix(nrow = length(genus_uniq), ncol=N))
for (i in 1: length(sum_L6)){
  myL6_percentage[,i] <- myL6[,i]/sum_L6[i]
}

most_dominant_stage_L6 <- vector()
for (i in 1:length(col_index_list)){
  most_dominant_stage_L6[i] <- which.max(apply(myL6_percentage[,col_index_list[[i]]], 1, mean))
}
L6_average <- apply(myL6_percentage,1,mean)
myTable_L6 <- data.frame(L6_average,myL6["Genus"])    # overal taxa ratio

myTable_L6_byStage <- data.frame(matrix(nrow=nrow(myL6_percentage), 
                                        ncol = length(col_index_list)))
for (i in 1:length(col_index_list)){
  myTable_L6_byStage[,i] <- apply(myL6_percentage[,col_index_list[[i]]], 1, mean)
}
myTable_L6_byStage[,13] <- myL6["Genus"]
apply(myTable_L6_byStage[,1:12],2,function(x) which(x == max(x)))

index_classified_genera <- as.numeric(row.names(myTable_L6_byStage[-grep("[[a-z]]", myTable_L6_byStage$Genus),]) )
1 - sum(myTable_L6_byStage[index_classified_genera,1:12])/sum(myTable_L6_byStage[,1:12])    # 39.86% reads were not classified at genus
# merge minor genera for plotting
index_to_sum_L6 <- which(apply(myTable_L6_byStage[1:12],1,max) < 0.02)
new_sum_L6 <- apply(myTable_L6_byStage[index_to_sum_L6, 1:12],2,sum)
plot_data_L6 <- rbind(myTable_L6_byStage[-(index_to_sum_L6),1:12], new_sum_L6)
names(plot_data_L6) <- stage_name_order
plot_data_L6$Genus <-c(myTable_L6_byStage[-(index_to_sum_L6),"Genus"], "Other/Unassigned")

plot_data_L6 <- plot_data_L6[,c(paste0("stage",1:12),"Genus")]

plot_data_L6$percentage <- apply(plot_data_L6[,1:12], 1, mean)
plot_data_L6 <- plot_data_L6[order(plot_data_L6$percentage),]
plot_data_L6 <- plot_data_L6[,1:13]
myTable_L6_byStage_longF <- melt(plot_data_L6, id.vars = "Genus")
myTable_L6_byStage_longF$Genus <- factor(myTable_L6_byStage_longF$Genus, levels=plot_data_L6$Genus)

p_L6 <- ggplot() + geom_bar(data=myTable_L6_byStage_longF, 
                            aes(y=value, x=variable, fill=Genus), colour="black",
                            stat="identity") +
  scale_fill_manual(values=col_vector[1:nrow(plot_data_L6)]) +
  labs(x="Management Unit ID", 
       y="Relative abundance (%)",
       main="Composition at genus") +
  scale_x_discrete(labels=1:12) +
  scale_y_continuous(labels = function(x) x*100, expand=c(0,0)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
  theme(
    axis.line = element_line(colour = 'black'),
    axis.text.x = element_text(size = 78, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 78, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(size = 84, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 84, face = "bold"),
    legend.title = element_text(size = 72, face="bold"),
    legend.text = element_text(size = 72, face="bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.direction = "vertical",
    legend.key.size =  unit(0.6, "in"),
    legend.key.height=unit(8,"line"),    
    legend.box.margin=margin(0,0,0,50),
    plot.margin=unit(c(3,0.5,0.8,0.5), "cm")) +
  guides(fill=guide_legend(ncol=2))
ggsave("Fig.2F_composition_at_genus.tiff", plot= p_L6, width=48, height=36, compression="lzw")


######################################################################
# No. of ASVs assigned to known genera and unknown genera
table(startsWith(x = my_table_2$Genus, prefix = "[")) # known: 2169, unknown: 2512
table(startsWith(x = myL6$Genus, prefix = "[")) # known: 167, unknown: 70

# check if any stage has unique genera
my_unique_genera <- data.frame(matrix(nrow=nrow(myL6), ncol=length(col_index_list)))
for (this_stage_for_genera in 1:length(col_index_list)){
  my_current_data_for_genera <- myL6[col_index_list[[this_stage_for_genera]]]
  my_unique_genera[, this_stage_for_genera] <- apply(my_current_data_for_genera, 1, sum)
}
my_unique_genera_binary <- ifelse(my_unique_genera == 0, 0, 1)

index_unique_in_one_stage_for_genera <- which(apply(my_unique_genera_binary, 1, sum) == 1)
my_unique_genera[index_unique_in_one_stage_for_genera,]
myL6[index_unique_in_one_stage_for_genera,]      # only one:  OTU#139 Prevotellaceae_NK3B31_group

index_unique_in_two_stage_for_genera <- which(apply(my_unique_genera_binary, 1, sum) == 2)
my_unique_genera[index_unique_in_two_stage_for_genera,]
myL6[index_unique_in_two_stage_for_genera,]     # only three: OTUs#140, 162, 233

hist((apply(my_unique_genera_binary, 1, sum)))

table_unique_unit_genera <- my_unique_genera[my_unique_genera$X1 > (MIN_DEPTH * N * 0.1 * 0.01),]

# # Venn Diagram for units 7-11
# apply(my_unique_genera_binary[,1:3],2,sum)
# genera_unit1 <- which(my_unique_genera_binary[,1] == 1)
# genera_unit2 <- which(my_unique_genera_binary[,2] == 1)
# genera_unit3 <- which(my_unique_genera_binary[,3] == 1)
# genera_unit1_list <- myL6[genera_unit1,"Genus"]
# genera_unit2_list <- myL6[genera_unit2,"Genus"]
# genera_unit3_list <- myL6[genera_unit3,"Genus"]
# venn.diagram(
#   x = list(genera_unit2_list , genera_unit3_list ),
#   category.names = c("Unit 2" , "Unit 3"),
#   filename = '#16_venn_diagramm.png',
#   fill = my_colors[c(2,3)])
# 
# apply(my_unique_genera_binary[,7:11],2,sum)
# genera_unit7 <- which(my_unique_genera_binary[,7] == 1)
# genera_unit8 <- which(my_unique_genera_binary[,8] == 1)
# genera_unit9 <- which(my_unique_genera_binary[,9] == 1)
# genera_unit10 <- which(my_unique_genera_binary[,10] == 1)
# genera_unit11 <- which(my_unique_genera_binary[,11] == 1)
# #
# genera_unit7_list <- myL6[genera_unit7,"Genus"]
# genera_unit8_list <- myL6[genera_unit8,"Genus"]
# genera_unit9_list <- myL6[genera_unit9,"Genus"]
# genera_unit10_list <- myL6[genera_unit10,"Genus"]
# genera_unit11_list <- myL6[genera_unit11,"Genus"]
# 
# venn.diagram(
#   x = list(genera_unit7_list , genera_unit8_list , genera_unit9_list,
#            genera_unit10_list,genera_unit11_list),
#   category.names = c("Unit 7" , "Unit 8" , "Unit 9", "Unit 10", "Unit 11"),
#   filename = '#14_venn_diagramm.png',
#   fill = my_colors[7:11]
#   # output = TRUE ,
  # imagetype="png" ,
  # height = 1200 , 
  # width = 800
  # resolution = 300,
  # compression = "lzw",
  # lwd = 2,
  # lty = 'blank',
  # cex = 1,
  # fontface = "bold",
  # fontfamily = "sans",
  # cat.cex = 0.6,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # # cat.pos = c(-27, 27, 135),
  # # cat.dist = c(0.055, 0.055, 0.085),
  # cat.fontfamily = "sans",
  # rotation = 1
# )

###########################################################################
#    20220311 --> 
#       assign 'major/minor' tag by ra & occurrences in no. of samples    #
#       average no. of reads each unit for each genera
###########################################################################

myL6$core1 <- "" # min 0.1% ra in at least 10 samples
myL6$core2 <- "" # min 0.1% ra in at least one unit 
myL6$assign <- "" # assigned vs. unassigned 

# update two cols: 'core1' & 'assign'
thrhld <- MIN_DEPTH * 0.001
for (i in 1:nrow(myL6)){
  tester <- sum(myL6[i,1:90] > thrhld)
  if (tester >= N * 0.1){
    myL6[i,"core1"] <- "major"
  } else if (tester < N * 0.1){
    myL6[i,"core1"] <- "minor"
  }
  
  if (startsWith(myL6[i,"Genus"],"[") == TRUE){
    myL6[i,"assign"] <- "unassigned"
  } else {myL6[i,"assign"] <- "assigned"}
}

# add cols of avg. by stage & update col 'core2'
add_col_names <- paste0("unit",1:12)
myL6[,add_col_names] <- 0
for (i in 1:length(col_index_list)){
  temp_df <- myL6[,col_index_list[[i]]]
  myL6[,add_col_names[i]] <- apply(temp_df,1,sum)/length(col_index_list[[i]])
}

for (i in 1:nrow(myL6)){
  this_value <- myL6[i,add_col_names]
  if (any(this_value >= thrhld)){
    myL6[i,"core2"] <- "major"
  } else if (all(this_value < thrhld)){
    myL6[i,"core2"] <- "minor"
  }
}
myL6$identifier1 <- paste0(myL6$core1,"_",myL6$assign)
myL6$identifier2 <- paste0(myL6$core2,"_",myL6$assign)
table(myL6$identifier1)
table(myL6$identifier2)

myL6[,add_col_names] <- myL6[,add_col_names] %>% mutate_if(is.numeric, round)

write.table(x = myL6,"myL6_20220311.txt",sep = "\t", quote = FALSE, row.names = FALSE)
myL6_major <- myL6[myL6$core2=="major",]
write.table(x = myL6_major,"myL6_major_20220311.txt",sep = "\t", quote = FALSE, row.names = FALSE)
myL6_minor <- myL6[myL6$core2=="minor",]
write.table(x = myL6_minor,"myL6_minor_20220311.txt",sep = "\t", quote = FALSE, row.names = FALSE)

# major
myL6_major_2 <- myL6_major[,paste0("unit",1:12)]
myL6_major_3 <- data.frame(matrix(nrow=nrow(myL6_major_2), ncol=ncol(myL6_major_2))) # myL6_major_3 is a binary version of myL6_major_2
for (i in 1:nrow(myL6_major_2)){
  myL6_major_3[i,] <- ifelse(myL6_major_2[i,]>0,1,0)
};
table(apply(myL6_major_3,1,sum)) # no. of genera in no. of stages
# 1  2  4  7  8  9 10 11 12 
# 1  7  3  7  3  3  6 12 70 

# minor
myL6_minor_2 <- myL6_minor[,paste0("unit",1:12)]
myL6_minor_3 <- data.frame(matrix(nrow=nrow(myL6_minor_2), ncol=ncol(myL6_minor_2))) # myL6_minor_3 is a binary version of myL6_minor_2
for (i in 1:nrow(myL6_minor_2)){
  myL6_minor_3[i,] <- ifelse(myL6_minor_2[i,]>0,1,0)
};
table(apply(myL6_minor_3,1,sum)) # no. of genera in no. of stages
#  0  1  2  3  4  5  6  7  8  9 10 11 12
# 22 22  7 18 12  6  9  6 11  4  2  3  3
myL6_minor_2[myL6_minor_3[which(apply(myL6_minor_3,1,sum) == 1),]$unit1 != 0,]




# sampleID_lower_cluster <- c("BH.008", "CUC.004","CUC.005","Hos.001","MLLH.003",
#                             "FC.006","FH.006","FC.005","FC.003","MLC.007",
#                             "FC.004","MLC.005","Hos.005", "FC.002")
# sampleID_left_cluster <- c("PWH.003","HC.002","HC.009","HC.004","HC.001","HC.003",
#                            "HC.006","HC.007", "HC.008", "PWH.001","HC.005")
# sampleID_upper_cluster <- setdiff(names(myL6)[1:90], c(sampleID_lower_cluster,sampleID_left_cluster))

# my_pvalues <- c()
# for (i in 1:nrow(myL6)){
#   my_upper <- myL6[i, sampleID_upper_cluster]
#   my_lower <- myL6[i, sampleID_lower_cluster]
#   my_pvalues[i] <- t.test(my_upper,my_lower)$p.value
# }
# my_pvalues <- ifelse(is.na(my_pvalues),1,my_pvalues)
# myL6$pvalues <- my_pvalues

# L6_by_unit_subs_int$pvalues <- 0
# for (i in 1:nrow(L6_by_unit_subs_int)){
#   this_taxon <- L6_by_unit_subs_int[i,"genus"]
#   row_index <- which(myL6$Genus == this_taxon)
#   L6_by_unit_subs_int[i,"pvalues"] <- myL6[row_index,"pvalues"]
# }

# L6_by_unit_subs_int$ra <- apply(L6_by_unit_subs_int[,1:12],1,mean)/MIN_DEPTH

# # plot genera that only appeared in unit 1 (not visually good, dump!)
# genera_only_unit_1 <- c("Prevotella_9","Fusobacterium","Prevotella","Subdoligranulum",
#                         "Prevotella_2","Sutterella","Oscillospira","Parabacteroides")
# plot_df_genera_unit_1 <- L6_by_unit_subs_int[L6_by_unit_subs_int$genus %in% genera_only_unit_1,]
# plot_df_genera_unit_1 <- plot_df_genera_unit_1[,1:13]
# plot_df_genera_unit_1_lf <- melt(plot_df_genera_unit_1, id.vars = "genus")
# g <- ggplot() +
#   geom_bar(data=plot_df_genera_unit_1_lf, aes(x=variable, y=log2(value+1), fill=variable), color="white" , stat="identity", size=3) +
#   theme_classic() +
#   facet_wrap(. ~ genus, nrow=2)

# End of 20220310
#############################################

