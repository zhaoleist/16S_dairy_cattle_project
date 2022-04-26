#############################################
#    Plot no. of ASVs assigned to family    #
#############################################
df_no_of_ASV_assigned_to_family <- data.frame(family=sort(table(my_table_2[,95]), decreasing = TRUE))
df_no_of_ASV_assigned_to_family <- df_no_of_ASV_assigned_to_family[1:20,]
g_no_of_ASV_assigned_to_family <- ggplot() +
  geom_bar(data=df_no_of_ASV_assigned_to_family, aes(x=family.Var1, y=family.Freq),fill="white", color="black", stat="identity", size=3) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 54, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 54, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Families (partial)") + 
  ylab("No. of OTUs") +
  annotate("text", x= 1:nrow(df_no_of_ASV_assigned_to_family),
           y=max(df_no_of_ASV_assigned_to_family$family.Freq) * 0.36,
           label=levels(df_no_of_ASV_assigned_to_family$family.Var1)[1:20], 
           angle=90, size=24, fontface="bold", hjust=0) +
  scale_y_continuous(limits=c(0, max(df_no_of_ASV_assigned_to_family$family.Freq) + 275), breaks=seq(0,max(df_no_of_ASV_assigned_to_family$family.Freq)+59, 500))
ggsave("Fig.2C.tiff", plot = g_no_of_ASV_assigned_to_family, width=32, height=24, dpi=300, compression="lzw")

############################################
#    Relative abundance at family level    #
############################################
family_uniq <- unique(my_table_2["Family"])[[1]]
myL5 <- data.frame(matrix(nrow = length(family_uniq), ncol=91))    # family table
for (i in 1:length(family_uniq)){
  index <- which(my_table_2["Family"]==as.character(family_uniq[i]))
  myL5[i,1:N] <- apply(my_table_2[index,1:N],2,sum)
  myL5[i,91] <- as.character(family_uniq[i])
  names(myL5) <- c(names(my_table_2[1:N]),"Family")
}

# No. of samples containing each of the 99 families
myL5_temp <- as.data.frame(ifelse(myL5[,1:90] == 0, 0, 1))
df_fig_temp <- data.frame(value=apply(myL5_temp,1,sum),
                        family=myL5$Family)
df_fig_temp$family <- factor(df_fig_temp$family, levels = levels(df_no_of_ASV_assigned_to_family$family.Var1))
df_fig_temp <- df_fig_temp[1:20,]

g_fig_temp <- ggplot() +
  geom_bar(data=df_fig_temp, aes(x=family, y=value), fill="white", col="black", stat="identity") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 28, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 28, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Families (patial)") + 
  ylab("No. of samples containing the family") +
  annotate("text", x= 1:nrow(df_fig_temp), 
           y=max(df_fig_temp$value) * 0.445, 
           label=levels(df_fig_temp$family)[1:20], angle=90, size=10, fontface="bold") +
  scale_y_continuous(limits=c(0, max(df_fig_temp$value) + 5), breaks=seq(0,max(df_fig_temp$value), 10))
# ggsave("Fig.999_family_in_no_of_samples.tiff", plot = g_fig_temp, width=16, height=12, dpi=300, compression="lzw")

# Relative abundance
index_shared_by_all_L5 <- which(apply(myL5[,1:90], 1, function(x) all(x != 0))  == TRUE)    # which families shared by all sample?
myL5[index_shared_by_all_L5, "Family"]
sum(apply(myL5[index_shared_by_all_L5, 1:90], 1, sum))/sum(apply(myL5[,1:90], 1, sum))    # ratio of shared phyla to all families: 0.8723036
sum_L5 <- apply(myL5[1:N],2,sum)

myL5_percentage <- data.frame(matrix(nrow = length(family_uniq), ncol=N))    # taxa percentage for each sample 
for (i in 1: length(sum_L5)){
  myL5_percentage[,i] <- myL5[,i]/sum_L5[i]
}
table(apply(myL5_percentage[,1:90],2,function(x) which(x == max(x))))

most_dominant_stage_L5 <- vector()
for (i in 1:length(col_index_list)){
  most_dominant_stage_L5[i] <- which.max(apply(myL5_percentage[,col_index_list[[i]]], 1, mean))    # same taxon (Ruminococcaceae) dominant in all units
}

L5_average <- apply(myL5_percentage,1,mean)
myTable_L5 <- data.frame(L5_average,myL5["Family"])

myTable_L5_byStage <- data.frame(matrix(nrow=nrow(myL5_percentage), 
                                        ncol = length(col_index_list)))
for (i in 1:length(col_index_list)){
  myTable_L5_byStage[,i] <- apply(myL5_percentage[,col_index_list[[i]]], 1, mean)
}
myTable_L5_byStage[,13] <- myL5["Family"]

# merge minor taxa together as "Other/Unassigned"
index_to_sum_L5 <- which(apply(myTable_L5_byStage[1:12],1,max) < 0.02)    
new_sum_L5 <- apply(myTable_L5_byStage[index_to_sum_L5, 1:12],2,sum)
plot_data_L5 <- rbind(myTable_L5_byStage[-(index_to_sum_L5),1:12], new_sum_L5)
names(plot_data_L5) <- stage_name_order
plot_data_L5$Family <-c(myTable_L5_byStage[-(index_to_sum_L5),"Family"], "Other/Unassigned")    
plot_data_L5 <- plot_data_L5[,c(paste0("stage",1:12),"Family")]
plot_data_L5$percentage <- apply(plot_data_L5[,1:12], 1, mean)
plot_data_L5 <- plot_data_L5[order(plot_data_L5$percentage),]
plot_data_L5 <- plot_data_L5[,1:13]

# change format for plotting
myTable_L5_byStage_longF <- melt(plot_data_L5, id.vars = "Family")
myTable_L5_byStage_longF$Family <- factor(myTable_L5_byStage_longF$Family, levels=plot_data_L5$Family)

p_L5 <- ggplot() + geom_bar(data=myTable_L5_byStage_longF, 
                            aes(y=value, x=variable, fill=Family), colour="black",
                            stat="identity") +
  scale_fill_manual(values=col_vector[1:nrow(plot_data_L5)]) +
  labs(x="Management Unit ID", 
       y="Relative abundance (%)",
       main="Composition at family") +
  scale_x_discrete(labels=1:12) +
  scale_y_continuous(labels = function(x) x*100, expand=c(0,0)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
  theme(
    axis.line = element_line(colour = 'black'),
    axis.text.x = element_text(size = 64, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 64, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(size = 72, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 72, face = "bold"),
    legend.title = element_text(size = 72, face="bold"),
    legend.text = element_text(size = 72, face="bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size =  unit(1, "in"),
    legend.box.margin=margin(0,0,0,50),
    legend.key.height=unit(3,"line"),
    plot.margin=unit(c(2,0.5,0.8,0.5), "cm"))
nrow(myTable_L5_byStage[-grep("[[a-z]]", myTable_L5_byStage$Family),])    # number of classified families out of 99
ggsave("Fig.2D_composition_at_family.tiff", plot= p_L5, width=32, height=24, compression="lzw")

