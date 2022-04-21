############################################
#    check if any stage has unique ASVs    #
############################################
my_uniq_OTUs <- data.frame(matrix(nrow=nrow(my_table_2),ncol=length(col_index_list)))
for (this_stage in 1:length(col_index_list)){
  this_data <- my_table_2[col_index_list[[this_stage]]]
  my_uniq_OTUs[, this_stage] <- apply(this_data, 1, sum)
}
# table(apply(ifelse(my_table_2[,1:90] == 0, 0, 1),1,sum))

thrhld <- 0.001
index_ASV_over_0.001 <- which(apply(my_uniq_OTUs,1,sum) > (sum(my_uniq_OTUs) * thrhld))
length(index_ASV_over_0.001)/nrow(my_table_2) # percentage of ASVs having > 0.1% abundance

my_uniq_OTUs_binary <- ifelse(my_uniq_OTUs == 0, 0, 1)
apply(my_uniq_OTUs_binary, 2, sum)    # number of OTU in each unit

index_unique_in_one_stage <- which(apply(my_uniq_OTUs_binary, 1, sum) == 1)
my_uniq_OTUs[index_unique_in_one_stage,]
my_table_2[index_unique_in_one_stage,]
my_table_2[index_unique_in_one_stage[1:4],c(1:9,91:96)]

index_unique_in_two_stage <- which(apply(my_uniq_OTUs_binary, 1, sum) == 2)
my_uniq_OTUs[index_unique_in_two_stage,]
my_table_2[index_unique_in_two_stage,]
my_table_2[index_unique_in_two_stage[1:4],c(1:17,91:96)]

# This plot is to show every ASV existed in which units
hist((apply(my_uniq_OTUs_binary, 1, sum)))
df_OTU_dist <- data.frame(matrix(table(apply(my_uniq_OTUs_binary, 1, sum))))
df_OTU_dist$X <- factor(1:12)
names(df_OTU_dist) <- c("count", "X")

g_OTU <- ggplot(df_OTU_dist, aes(x=X, y=count)) + geom_bar(stat="identity", fill="white", colour="black", size=3) +
  geom_text(aes(label=count), vjust=-0.3, size=12, fontface="bold") +
  labs(x="Number of Management Units",
       y="Count of ASVs") +
  scale_y_continuous(breaks=seq(0,580,100), limits = c(0,max(df_OTU_dist$count))) +
  theme_classic() +
  # theme(
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   panel.background = element_blank(),
  #   panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
  theme(
    axis.line = element_line(colour = 'black'),
    axis.text.x = element_text(size = 36, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 36, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 40, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    plot.margin=unit(c(0.8,0.5,0.8,0.5), "cm"))    # t,r,b,l
# ggsave(g_OTU, file="Fig.999_OTU_in_no_of_stage.tiff", width=16, height=12, compression='lzw')

# loosen the criteria of defining "unique" and found 19 unique ASVs that only detected in mu 1
table_unique_unit <- my_uniq_OTUs[my_uniq_OTUs$X1 > (MIN_DEPTH * 9 * thrhld),]
table_unique_unit_filtered <- table_unique_unit[which(apply(table_unique_unit[2:12],1,sum) < 3),]
list_of_unique_ASV <- my_table_2[as.numeric(row.names(table_unique_unit_filtered)),]
row.names(list_of_unique_ASV)
# write.table(list_of_unique_ASV, file="unique_ASVs.txt", sep="\t", row.names = FALSE)
# no such unique ASVs that can be found in any of other units

# another three ASVs that only found in units 1 & 2
table_unique_unit_1_2 <- my_uniq_OTUs[my_uniq_OTUs$X1 > (MIN_DEPTH * lengths(col_index_list)[1] * thrhld) & my_uniq_OTUs$X2 > (MIN_DEPTH * lengths(col_index_list)[2] * thrhld),] # b/c my_uniq_OTUs is sum, not avg.
table_unique_unit_filtered_1_2 <- table_unique_unit_1_2[which(apply(table_unique_unit_1_2[3:12],1,sum) < 3),]
list_of_unique_ASV_1_2 <- my_table_2[as.numeric(row.names(table_unique_unit_filtered_1_2)),] 


#################################
##    20220311 --> dynamics    ##
#################################
alpha <- 0.05
df_for_dynamic <- myL6[,1:90]
pValues <- c()
pValues_adjusted <- c()
tk_result <- list()
for (i in 1:nrow(df_for_dynamic)){
  temp_df <- data.frame(data=as.numeric(unlist(df_for_dynamic[i,])), Stage=stage)
  pValues[i] <- anova(lm(temp_df$data ~ temp_df$Stage))$"Pr(>F)"[1]
  tk_result[i] <- TukeyHSD(aov(lm(temp_df$data ~ temp_df$Stage)))
}
pValues_adjusted <- p.adjust(pValues, method="BH")
index_sig <- which(pValues_adjusted < alpha)
index_to_keep <- which(myL6$core2 == "major")

index_sig_dynamic <- intersect(index_to_keep, index_sig) # 56

df_for_dynamic <- df_for_dynamic[index_sig_dynamic,1:90]    # df_for_dynamic now is OTUs that with significant p-values, 0.1% abundance
# 56 x 90
pValues_adjusted_for_212 <- pValues_adjusted[index_sig_dynamic]

tk_result <- tk_result[index_sig_dynamic]
N_Tukey <- c() # For the 56 active OTUs, N_Tukey is the # of pairs that significantly different. Max = 66 (choose(12,2))
tk_result_sig <- list()

for (i in 1: length(tk_result)){
  N_Tukey[i] <- nrow(as.data.frame(tk_result[i])[which(as.data.frame(tk_result[i])[4] < alpha),])
  tk_result_sig[[i]] <- which(as.data.frame(tk_result[i])[4] < alpha)
}
# check if the range contains zero
check_zero <- list()
# (as.data.frame(tk_result[1])[tk_result_sig[[1]],2] < 0) & (as.data.frame(tk_result[1])[tk_result_sig[[1]],3] > 0)
for (i in 1:length(tk_result_sig)){
  check_zero[[i]] <- (as.data.frame(tk_result[i])[tk_result_sig[[i]],2] < 0) & 
    (as.data.frame(tk_result[i])[tk_result_sig[[i]],3] > 0)
}
sum(unlist(check_zero))
# factor(unlist(tk_result_sig), levels=1:66)
# table(factor(unlist(tk_result_sig), levels=1:66))

btw_unit <- data.frame(id=row.names(as.data.frame(tk_result[1])),
                       freq=table(factor(unlist(tk_result_sig), levels=1:66)))
btw_unit_final <- data.frame(matrix(nrow=12, ncol=12))
names(btw_unit_final) <- paste0("",1:12)
row.names(btw_unit_final) <- paste0("",1:12)
diag(btw_unit_final) <- 66

btw_unit_final[1,2:12] <- c(btw_unit$freq.Freq[1:11])
btw_unit_final[2:12,1] <- c(btw_unit$freq.Freq[1:11])

btw_unit_final[2,3:12] <- btw_unit$freq.Freq[12:21]
btw_unit_final[3:12,2] <- btw_unit$freq.Freq[12:21]

btw_unit_final[3,4:12] <- btw_unit$freq.Freq[22:30]
btw_unit_final[4:12,3] <- btw_unit$freq.Freq[22:30]

btw_unit_final[4,5:12] <- btw_unit$freq.Freq[31:38]
btw_unit_final[5:12,4] <- btw_unit$freq.Freq[31:38]

btw_unit_final[5,6:12] <- btw_unit$freq.Freq[39:45]
btw_unit_final[6:12,5] <- btw_unit$freq.Freq[39:45]

btw_unit_final[6,7:12] <- btw_unit$freq.Freq[46:51]
btw_unit_final[7:12,6] <- btw_unit$freq.Freq[46:51]

btw_unit_final[7,8:12] <- btw_unit$freq.Freq[52:56]
btw_unit_final[8:12,7] <- btw_unit$freq.Freq[52:56]

btw_unit_final[8,9:12] <- btw_unit$freq.Freq[57:60]
btw_unit_final[9:12,8] <- btw_unit$freq.Freq[57:60]

btw_unit_final[9,10:12] <- btw_unit$freq.Freq[61:63]
btw_unit_final[10:12,9] <- btw_unit$freq.Freq[61:63]

btw_unit_final[10,11:12] <- btw_unit$freq.Freq[64:65]
btw_unit_final[11:12,10] <- btw_unit$freq.Freq[64:65]

btw_unit_final[11,12] <- btw_unit$freq.Freq[66]
btw_unit_final[12,11] <- btw_unit$freq.Freq[66]

diag(btw_unit_final) <- NA

tiff(filename = "Fig.5A_between_units_sig_2022.tiff", width=16, height=12, units="in", res=300, compression="lzw")
# Heatmap(as.matrix(btw_unit_final), cluster_rows = FALSE, cluster_columns = FALSE,
#         heatmap_legend_param = list(legend_height=unit(8,"cm"),
#                                     grid_width=unit(2,"cm"),
#                                     title="", fontsize=20,
#                                     labels_gp=gpar(fontsize=24),
#                                     grid_border=5,
#                                     at=seq(0,70,10)),
#         border = TRUE, row_names_gp = gpar(fontsize=28),
#         column_names_gp = gpar(fontsize=28),
#         col=c("darkblue","lightgreen","yellow"))
superheat(as.matrix(btw_unit_final),
          left.label.text.size = 12,
          bottom.label.text.size = 12,
          bottom.label.text.angle = 0,
          bottom.label.size = 0.1,
          left.label.size=0.05,
          # left.label.text.alignment = "right", 
         # bottom.label.text.alignment = "right",
          legend.height = 0.4,
          legend.width = 3,
          legend.text.size = 28,
          # legend.num.ticks = 3,
          # legend.vspace = -0.1,
          row.title = "Management Unit ID",
          column.title = "Management Unit ID",
          row.title.size = 12,
          column.title.size = 12)
dev.off()

# index_for_taxa <- as.numeric(rownames(df_for_dynamic))
# df_for_dynamic[,91:96] <- my_table_2_raw[index_for_taxa,91:96]

a_df <- data.frame(matrix(nrow=nrow(btw_unit), ncol=2))
for (i in 1:nrow(btw_unit)){
  these_units <- paste0(unlist(strsplit(btw_unit[i,"id"],split="-"))[c(1,3)],collapse="-")
  a_df[i,1] <- these_units
  index_vec <- c()
  for (j in 1:length(tk_result_sig)){
    if (i %in% tk_result_sig[[j]]){
      index_vec <- c(index_vec,j)
    }
  }
  a_df[i,2] <- paste0(index_vec, collapse=";")
}

names(a_df) <- c("paired_units", "index_in_index_sig_dynamic")
a_df$row_in_myL6 <- ""
for (i in 1:nrow(a_df)){
  index_vec <- a_df[i,"index_in_index_sig_dynamic"]
  index_vec2 <- as.numeric(unlist(strsplit(index_vec,split=";")))
  rows_in_myL6 <- index_sig_dynamic[index_vec2]
  a_df[i,"row_in_myL6"] <- paste0(rows_in_myL6, collapse=";")
}

b_df <- data.frame(matrix(nrow=nrow(a_df)*2,ncol=3))
names(b_df) <- c("denom_unit", "pair_unit", "value")
for (i in 1:nrow(a_df)){
  rows_to_use <- as.numeric(unlist(strsplit(a_df[i,"row_in_myL6"],split=";")))
  temp_df <- myL6[rows_to_use,paste0("unit",1:12)]
  two_units <- as.numeric(unlist(strsplit(a_df[i,"paired_units"], split="-")))
  b_df[i*2-1,"value"] <- sum(temp_df[,two_units[1]])
  b_df[i*2-1,"denom_unit"] <- names(temp_df)[two_units[1]]
  b_df[i*2-1,"pair_unit"] <- names(temp_df)[two_units[2]]
  
  b_df[i*2,"value"] <- sum(temp_df[,two_units[2]])
  b_df[i*2,"denom_unit"] <- names(temp_df)[two_units[2]]
  b_df[i*2,"pair_unit"] <- names(temp_df)[two_units[1]]
}

b_df$denom_unit <- factor(b_df$denom_unit, levels=paste0("unit",1:12))
b_df$pair_unit <- factor(b_df$pair_unit, levels=paste0("unit",1:12))
supp.labs <- paste("unit",1:12)
names(supp.labs) <- paste0("unit",1:12)

g <- ggplot(b_df) +
  geom_bar(aes(x=pair_unit,y=value/MIN_DEPTH*100,fill=pair_unit), color="black", stat="identity") + 
  facet_wrap(. ~ denom_unit, ncol=3, labeller = labeller(denom_unit = supp.labs)) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.text.y = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 28, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 28, face = "bold"),
    legend.title = element_text(size=24, face="bold"),
    legend.text = element_text(size=24),
    legend.text.align = 0.5,
    # legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    strip.text.x = element_text(size = 24, colour = "black", face="bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(x="Management Unit ID", y="Relative Abundance (%)", fill="Mgmt \nUnit ID") +
  scale_x_discrete(labels=c("unit1" = "1",
                            "unit2" = "2",
                            "unit3" = "3",
                            "unit4" = "4",
                            "unit5" = "5",
                            "unit6" = "6",
                            "unit7" = "7",
                            "unit8" = "8",
                            "unit9" = "9",
                            "unit10" = "10",
                            "unit11" = "11",
                            "unit12" = "12")) +
  scale_fill_manual(labels = 1:12, values=my_colors)
ggsave("Fig.5B_active_genera_bar.tiff", plot = g, width=20, height=16, dpi=300, compression="lzw")

# ############################################################
# df_212_ASVs_output <- df_for_dynamic
# df_212_ASVs_output$pValues <- pValues_adjusted_for_212
# df_212_ASVs_output$N_Pairs <- N_Tukey
# df_212_ASVs_output$ASV_ID <- row.names(df_212_ASVs_output)
# row.names(df_212_ASVs_output) <- 1:nrow(df_212_ASVs_output)
# ############################################################
# 
# df_otu_sig_phylum <- as.data.frame(table(df_for_dynamic$Phylum))
# df_otu_sig_phylum <- df_otu_sig_phylum[order(df_otu_sig_phylum$Freq, decreasing = TRUE),]    # descending order
# df_otu_sig_phylum$Var1 <- factor(df_otu_sig_phylum$Var1, levels=unique(df_otu_sig_phylum$Var1))    # change order of levels
# df_otu_sig_phylum$Letter <- factor(LETTERS[1:nrow(df_otu_sig_phylum)])
# 
# g_dynamics <- ggplot(data = df_otu_sig_phylum,aes(x=Letter, y=Freq)) + 
#   geom_bar(aes(fill=Var1),stat="identity", fill="white", color="black", width=0.7, size=4) +
#   geom_text(aes(label=Freq), vjust=-0.3, size=24) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(size=64, face="bold", colour="black"),
#     axis.text.y = element_text(size = 64, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
#     axis.title.x = element_text(size = 72, face = "bold", margin = margin(b = 10)),
#     axis.title.y = element_text(size = 72, face = "bold")) +
#   annotate("text", x= nrow(df_otu_sig_phylum)-2, y=seq(130,50,length.out = 7), 
#            label=paste0(LETTERS[0:7],": ",levels(df_otu_sig_phylum$Var1)), angle=0, size=28, fontface="plain", hjust=0) +
#   xlab("Phylum") + 
#   ylab("Frequency") +
#   ylim(0,135)
# # g_dynamics
# ggsave("FIG.7A_otu_sig_phylum.tiff", plot = g_dynamics, width = 32, height = 16, compression="lzw")
# 
# # genera in phyla Firmicutes
# LETTERS_EXT <- factor(c(LETTERS, letters[1:2]), levels=c(LETTERS, letters[1:2]))
# xx <- c(11, 7) # for adjusting legend position of fig.7b, 7c
# names(xx) <- c("Firmicutes", "Bacteroidetes")
# yy <- c(42, 27)
# names(yy) <- c("Firmicutes", "Bacteroidetes")
# zz <- c(5, 8)
# names(zz) <- c("Firmicutes", "Bacteroidetes")
# myHeight <- c(48,32)
# names(myHeight) <- c("Firmicutes", "Bacteroidetes")
# for (phylum_name in c("Firmicutes", "Bacteroidetes")){
#   df_otu_sig_genus <- as.data.frame(table(df_for_dynamic[df_for_dynamic$Phylum == phylum_name,"Genus"]))
#   df_otu_sig_genus <- df_otu_sig_genus[order(df_otu_sig_genus$Freq, decreasing = TRUE),] 
#   df_otu_sig_genus$Var1 <- factor(df_otu_sig_genus$Var1, levels=unique(df_otu_sig_genus$Var1))
# # num_other_genera <- sum(df_otu_sig_genus[11:nrow(df_otu_sig_genus),"Freq"])
# # df_otu_sig_genus_plot <- data.frame(Var1=c(as.character(unlist(df_otu_sig_genus[1:10,"Var1"])), "Other"),
# #                                  Freq=c(as.numeric(unlist(df_otu_sig_genus[1:10,"Freq"])), num_other_genera))
#   df_otu_sig_genus_plot <- df_otu_sig_genus
#   df_otu_sig_genus_plot$Var1 <- factor(df_otu_sig_genus_plot$Var1, levels=unique(df_otu_sig_genus_plot$Var1))
#   df_otu_sig_genus_plot$Letter <- LETTERS_EXT[0:nrow(df_otu_sig_genus_plot)]
#   
#   g_dynamics_2 <- ggplot(data = df_otu_sig_genus_plot) + 
#     geom_bar(aes(x=Letter, y=Freq), stat="identity", fill="white", color="black", size=3) +
#     annotate("text", x=1:nrow(df_otu_sig_genus_plot),y=df_otu_sig_genus_plot$Freq,label=as.double(df_otu_sig_genus_plot$Freq), vjust=-0.2, size=24) +
#     theme_classic() +
#     theme(
#       axis.text.x = element_text(size=68, face="bold", colour="black"),
#       axis.text.y = element_text(size = 72, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
#       axis.title.x = element_text(size = 72, face = "bold", margin = margin(t = 20)),
#       axis.title.y = element_text(size = 72, face = "bold"),
#       plot.margin = unit(c(1,1,1,1), "cm")) +
#     annotate("text", x= as.numeric(xx[phylum_name]), y=seq(as.numeric(yy[phylum_name]), as.numeric(zz[phylum_name]), length.out = nrow(df_otu_sig_genus_plot)), 
#              label=paste0(LETTERS_EXT[0:nrow(df_otu_sig_genus_plot)], ": ", levels(df_otu_sig_genus_plot$Var1)), angle=0, size=26, fontface="plain", hjust=0) +
#     xlab(paste("Genus in phylum ", phylum_name, sep="")) + 
#     ylab("Frequency") +
#     scale_y_continuous(limits=c(min(df_otu_sig_genus_plot$Freq)-1, max(df_otu_sig_genus_plot$Freq) + 10))
#   ggsave(paste("otu_sig_",phylum_name, ".tiff", sep=""), plot = g_dynamics_2, width = 32, height = myHeight[phylum_name], compression="lzw")
#   assign(paste("g", phylum_name, sep=""),g_dynamics_2)
# }
# # ggsave("Fig7.png", grid.arrange(g_dynamics, gFirmicutes, gBacteroidetes, nrow=3), width=40,height=30)


# ####################
# #    Plot fig.8    #
# ####################
# # phylum
# df_fig_8 <- my_table_2[index_sig,]
# df_fig_8a_L2_unique <- unique(df_fig_8["Phylum"])[[1]]
# df_fig_8a <- data.frame(matrix(nrow=length(df_fig_8a_L2_unique), ncol=91))
# for (i in 1:length(df_fig_8a_L2_unique)){
#   index <- which(df_fig_8["Phylum"] == as.character(df_fig_8a_L2_unique[i]))
#   df_fig_8a[i,1:N] <- apply(df_fig_8[index,1:N], 2, sum)
#   df_fig_8a[i, 91] <- as.character(df_fig_8a_L2_unique[i])
#   names(df_fig_8a) <- c(names(df_fig_8[1:N]),"Phylum")
# }
# df_fig_8a[,1:90] <- df_fig_8a[,1:90]/MIN_DEPTH
# average_fig_8a <- apply(df_fig_8a[,1:90],1,mean)
# names(average_fig_8a) <- df_fig_8a$Phylum
# 
# df_fig_8a_by_stage <- data.frame(matrix(nrow=nrow(df_fig_8a), 
#                                         ncol= length(col_index_list)))
# for (i in 1:length(col_index_list)){
#   df_fig_8a_by_stage[,i] <- apply(df_fig_8a[,col_index_list[[i]]], 1, mean)
# }
# df_fig_8a_by_stage[,13] <- df_fig_8a["Phylum"]  
# index_8a_merge <- which(apply(df_fig_8a_by_stage[,1:12],1,max) < 0.02)
# other_8a <- apply(df_fig_8a_by_stage[index_8a_merge,1:12],2,sum)
# plot_df_fig_8a_by_stage <-rbind(df_fig_8a_by_stage[-index_8a_merge,1:12], other_8a)
# names(plot_df_fig_8a_by_stage) <- paste0("stage", 1:12)
# plot_df_fig_8a_by_stage$Phylum <- c(df_fig_8a_by_stage[-index_8a_merge,"Phylum"], "Other/Unassigned")
# plot_df_fig_8a_by_stage$Phylum <- factor(plot_df_fig_8a_by_stage$Phylum,
#                                          levels=rev(c("Firmicutes","Bacteroidetes","Patescibacteria",
#                                                   "Verrucomicrobia", "Other/Unassigned")))
# 
# plot_df_fig_8a_by_stage_longF <- melt(plot_df_fig_8a_by_stage, id.vars = "Phylum")
# 
# g_fig_8a <- ggplot() + 
#   geom_bar(data=plot_df_fig_8a_by_stage_longF,
#            aes(x=variable, y=value, fill=Phylum), stat="identity") +
#   scale_fill_manual(values=col_vector) +
#   labs(x="",
#        y="Relative abundance (%)",
#        main="Composition at phylum") +
#   scale_x_discrete(labels=1:12) +
#   # scale_y_continuous(labels = scales::percent) +
#   scale_y_continuous(labels = function(x) x*100, expand=c(0,0)) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
#   theme(
#     axis.line = element_line(colour = 'black'),
#     axis.text.x = element_text(size = 72, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
#     axis.text.y = element_text(size = 72, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
#     axis.title.x = element_text(size = 78, face = "bold", margin = margin(t = 25, r = 0, b = 0, l = 0)),
#     axis.title.y = element_text(size = 78, face = "bold", margin = margin(t = 0, r = 25, b = 0, l = 0)),
#     legend.title = element_text(size = 72, face="bold"),
#     legend.text = element_text(size = 72, face="bold"),
#     legend.background = element_blank(),
#     legend.key = element_blank(),
#     legend.direction = "vertical",
#     legend.key.size =  unit(0.5, "in"),
#     legend.key.height=unit(5,"line"),    # legend space
#     plot.margin=unit(c(0.8,0.5,0.8,0.8), "cm"))    # t,r,b,l
# ggsave("FIG.8A.tiff", plot=g_fig_8a, width=48, height=18, compression="lzw")
# ############################################################################
# # family
# df_fig_8b_L5_unique <- unique(df_fig_8["Family"])[[1]]
# df_fig_8b <- data.frame(matrix(nrow=length(df_fig_8b_L5_unique), ncol=91))
# for (i in 1:length(df_fig_8b_L5_unique)){
#   index <- which(df_fig_8["Family"] == as.character(df_fig_8b_L5_unique[i]))
#   df_fig_8b[i,1:N] <- apply(df_fig_8[index,1:N], 2, sum)
#   df_fig_8b[i, 91] <- as.character(df_fig_8b_L5_unique[i])
#   names(df_fig_8b) <- c(names(df_fig_8[1:N]),"Family")
# }
# df_fig_8b[,1:90] <- df_fig_8b[,1:90]/MIN_DEPTH
# average_fig_8b <- apply(df_fig_8b[,1:90],1,mean)
# names(average_fig_8b) <- df_fig_8b$Family
# 
# df_fig_8b_by_stage <- data.frame(matrix(nrow=nrow(df_fig_8b), 
#                                         ncol= length(col_index_list)))
# for (i in 1:length(col_index_list)){
#   df_fig_8b_by_stage[,i] <- apply(df_fig_8b[,col_index_list[[i]]], 1, mean)
# }
# df_fig_8b_by_stage[,13] <- df_fig_8b["Family"]  
# index_8b_merge <- which(apply(df_fig_8b_by_stage[,1:12],1,max) < 0.02)
# other_8b <- apply(df_fig_8b_by_stage[index_8b_merge,1:12],2,sum)
# plot_df_fig_8b_by_stage <-rbind(df_fig_8b_by_stage[-index_8b_merge,1:12], other_8b)
# names(plot_df_fig_8b_by_stage) <- paste0("stage", 1:12)
# plot_df_fig_8b_by_stage$Family <- c(df_fig_8b_by_stage[-index_8b_merge,"Family"], "Other/Unassigned")
# plot_df_fig_8b_by_stage <- plot_df_fig_8b_by_stage[order(apply(plot_df_fig_8b_by_stage[,1:12],1,mean)),]
# plot_df_fig_8b_by_stage$Family <- factor(plot_df_fig_8b_by_stage$Family,
#                                           levels=plot_df_fig_8b_by_stage[order(apply(plot_df_fig_8b_by_stage[,1:12],1,mean)),13])
# 
# plot_df_fig_8b_by_stage_longF <- melt(plot_df_fig_8b_by_stage, id.vars = "Family")
# 
# g_fig_8b <- ggplot() + 
#   geom_bar(data=plot_df_fig_8b_by_stage_longF,
#            aes(x=variable, y=value, fill=Family), stat="identity") +
#   scale_fill_manual(values=col_vector) +
#   labs(x="",
#        y="Relative abundance (%)",
#        main="Composition at family") +
#   scale_x_discrete(labels=1:12) +
#   # scale_y_continuous(labels = scales::percent) +
#   scale_y_continuous(labels = function(x) x*100, expand=c(0,0)) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
#   theme(
#     axis.line = element_line(colour = 'black'),
#     axis.text.x = element_text(size = 72, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
#     axis.text.y = element_text(size = 72, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
#     axis.title.x = element_text(size = 78, face = "bold", margin = margin(t = 25, r = 0, b = 0, l = 0)),
#     axis.title.y = element_text(size = 78, face = "bold", margin = margin(t = 0, r = 25, b = 0, l = 0)),
#     legend.title = element_text(size = 72, face="bold"),
#     legend.text = element_text(size = 72, face="bold"),
#     legend.background = element_blank(),
#     legend.key = element_blank(),
#     legend.direction = "vertical",
#     legend.key.size =  unit(0.5, "in"),
#     legend.key.height=unit(5,"line"),    # legend space
#     plot.margin=unit(c(0.8,0.5,0.8,0.8), "cm"))    # t,r,b,l
# ggsave("FIG.8B-.tiff", plot=g_fig_8b, width=48, height=18, compression="lzw")
# ############################################################################
# # genus
# df_fig_8c_L6_unique <- unique(df_fig_8["Genus"])[[1]]
# df_fig_8c <- data.frame(matrix(nrow=length(df_fig_8c_L6_unique), ncol=91))
# for (i in 1:length(df_fig_8c_L6_unique)){
#   index <- which(df_fig_8["Genus"] == as.character(df_fig_8c_L6_unique[i]))
#   df_fig_8c[i,1:N] <- apply(df_fig_8[index,1:N], 2, sum)
#   df_fig_8c[i, 91] <- as.character(df_fig_8c_L6_unique[i])
#   names(df_fig_8c) <- c(names(df_fig_8[1:N]),"Genus")
# }
# df_fig_8c[,1:90] <- df_fig_8c[,1:90]/MIN_DEPTH
# average_fig_8c <- apply(df_fig_8c[,1:90],1,mean)
# names(average_fig_8c) <- df_fig_8c$Genus
# 
# df_fig_8c_by_stage <- data.frame(matrix(nrow=nrow(df_fig_8c), 
#                                         ncol= length(col_index_list)))
# for (i in 1:length(col_index_list)){
#   df_fig_8c_by_stage[,i] <- apply(df_fig_8c[,col_index_list[[i]]], 1, mean)
# }
# df_fig_8c_by_stage[,13] <- df_fig_8c["Genus"]  
# index_8c_merge <- which(apply(df_fig_8c_by_stage[,1:12],1,max) < 0.02)
# other_8c <- apply(df_fig_8c_by_stage[index_8c_merge,1:12],2,sum)
# plot_df_fig_8c_by_stage <-rbind(df_fig_8c_by_stage[-index_8c_merge,1:12], other_8c)
# names(plot_df_fig_8c_by_stage) <- paste0("stage", 1:12)
# plot_df_fig_8c_by_stage$Genus <- c(df_fig_8c_by_stage[-index_8c_merge,"Genus"], "Other/Unassigned")
# plot_df_fig_8c_by_stage <- plot_df_fig_8c_by_stage[order(apply(plot_df_fig_8c_by_stage[,1:12],1,mean)),]
# plot_df_fig_8c_by_stage$Genus <- factor(plot_df_fig_8c_by_stage$Genus,
#                                          levels=plot_df_fig_8c_by_stage[order(apply(plot_df_fig_8c_by_stage[,1:12],1,mean)),13])
# 
# plot_df_fig_8c_by_stage_longF <- melt(plot_df_fig_8c_by_stage, id.vars = "Genus")
# 
# g_fig_8c <- ggplot() + 
#   geom_bar(data=plot_df_fig_8c_by_stage_longF,
#            aes(x=variable, y=value, fill=Genus), stat="identity") +
#   scale_fill_manual(values=col_vector) +
#   labs(x="Management Unit ID",
#        y="Relative abundance (%)",
#        main="Composition at genus") +
#   scale_x_discrete(labels=1:12) +
#   # scale_y_continuous(labels = scales::percent) +
#   scale_y_continuous(labels = function(x) x*100, expand=c(0,0)) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
#   theme(
#     axis.line = element_line(colour = 'black'),
#     axis.text.x = element_text(size = 72, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
#     axis.text.y = element_text(size = 72, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
#     axis.title.x = element_text(size = 78, face = "bold", margin = margin(t = 25, r = 0, b = 0, l = 0)),
#     axis.title.y = element_text(size = 78, face = "bold", margin = margin(t = 0, r = 25, b = 0, l = 0)),
#     legend.title = element_text(size = 64, face="bold"),
#     legend.text = element_text(size = 64, face="bold"),
#     legend.background = element_blank(),
#     legend.key = element_blank(),
#     legend.direction = "vertical",
#     legend.key.size =  unit(0.5, "in"),
#     legend.key.height=unit(4,"line"),    # legend space
#     plot.margin=unit(c(0.8,0.5,0.8,0.8), "cm"))    # t,r,b,l
# ggsave("FIG.8C.tiff", plot=g_fig_8c, width=48, height=18, compression="lzw")
# 
# 
# 
# # K-means troubled here b/c most OTUs clustered into one big cluster
# df_kmeans <- data.frame(matrix(nrow=nrow(df_for_dynamic), ncol=length(col_index_list)))
# for (i in 1:length(col_index_list)){
#   df_kmeans[,i] <- apply(df_for_dynamic[,col_index_list[[i]]], 1, mean)
# }
# 
# pdf(paste("OTU#","125",".pdf",sep=""))
# par(mfrow=c(5,2))
# for (a in 1:nrow(df_kmeans)){
#   plot(1:12, df_kmeans[a,], type="l")
#   }
# dev.off()

