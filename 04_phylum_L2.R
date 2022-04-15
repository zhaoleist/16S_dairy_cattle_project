##########################################################
#    Plot no. of ASVs assigned to phyla    #
##########################################################
df_no_of_ASV_assigned_to_phyla <- data.frame(phylum=sort(table(my_table_2[,92]), decreasing = TRUE))
g_no_of_ASV_assigned_to_phyla <- ggplot() +
  geom_bar(data=df_no_of_ASV_assigned_to_phyla, aes(x=phylum.Var1, y=phylum.Freq),fill="white", color="black", stat="identity", size=3) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 54, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 54, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Phyla") + 
  ylab("No. of OTUs") +
  annotate("text", x= 1:nrow(df_no_of_ASV_assigned_to_phyla), 
           y=max(df_no_of_ASV_assigned_to_phyla$phylum.Freq) * 0.115, 
           label=levels(df_no_of_ASV_assigned_to_phyla$phylum.Var1), 
           angle=90, size=26, fontface="bold",hjust=0) +
  scale_y_continuous(limits=c(0, max(df_no_of_ASV_assigned_to_phyla$phylum.Freq) + 10), breaks=seq(0,max(df_no_of_ASV_assigned_to_phyla$phylum.Freq), 500))
ggsave("Fig.2A.tiff", plot = g_no_of_ASV_assigned_to_phyla, width=32, height=24, dpi=300, compression="lzw")

############################################
#    Relative abundance at phylum level    #
############################################
phylum_uniq <- unique(my_table_2["Phylum"])[[1]]
myL2 <- data.frame(matrix(nrow = length(phylum_uniq), ncol=91)) # phylum-sample table (not otu table)
for (i in 1:length(phylum_uniq)){ # merge OTUs w/ the same phyla
  index <- which(my_table_2["Phylum"]==as.character(phylum_uniq[i]))
  myL2[i,1:N] <- apply(my_table_2[index,1:N],2,sum)
  myL2[i,91] <- as.character(phylum_uniq[i])
  names(myL2) <- c(names(my_table_2[1:N]),"Phylum")
}

# No. of samples containing each of the 20 phyla
myL2_temp <- as.data.frame(ifelse(myL2[,1:90] == 0, 0, 1))
df_fig_temp <- data.frame(value=apply(myL2_temp,1,sum),
                        phylum=myL2$Phylum)
df_fig_temp$phylum <- factor(df_fig_temp$phylum, levels = levels(df_no_of_ASV_assigned_to_phyla$phylum.Var1))

g_fig_temp <- ggplot() +
  geom_bar(data=df_fig_temp, aes(x=phylum, y=value), fill="white", col="black", stat="identity") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 10), colour = "black"),
    axis.title.x = element_text(size = 28, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 28, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Phyla") + 
  ylab("No. of samples containing the phylum") +
  annotate("text", x= 1:nrow(df_fig_temp), 
           y=max(df_fig_temp$value) * 0.445, 
           label=levels(df_fig_temp$phylum), angle=90, size=10, fontface="bold") +
  scale_y_continuous(limits=c(0, max(df_fig_temp$value) + 5), breaks=seq(0,max(df_fig_temp$value), 10))
# ggsave("Fig.999_phylum_in_no_of_samples.tiff", plot = g_fig_temp, width=16, height=12, dpi=300, compression="lzw")

# Relative abundance
index_shared_by_all_L2 <- which(apply(myL2[,1:90], 1, function(x) all(x != 0))  == TRUE) # which phyla shared by all sample?
myL2[index_shared_by_all_L2, "Phylum"] # "Patescibacteria" "Firmicutes" "Bacteroidetes" "Proteobacteria" "[k]Bacteria" 
sum(apply(myL2[index_shared_by_all_L2, 1:90], 1, sum))/sum(apply(myL2[,1:90], 1, sum)) # ratio of reads in these shared phyla to reads in all phyla = 0.9452208
sum_L2 <- apply(myL2[1:N],2,sum) # overall count each sample: should equal 12637

myL2_percentage <- data.frame(matrix(nrow = length(phylum_uniq), ncol=N)) # taxa percentage for each sample 
for (i in 1: length(sum_L2)){
  myL2_percentage[,i] <- myL2[,i]/sum_L2[i]
}

table(apply(myL2_percentage[,1:90],2,function(x) which(x == max(x))))

L2_average <- apply(myL2_percentage,1,mean)    
myTable_L2 <- data.frame(L2_average,myL2["Phylum"]) # overall percentage for each phylum
 
myTable_L2_by_stage <- data.frame(matrix(nrow=nrow(myL2_percentage), # stage-wise percentage
                                        ncol = length(col_index_list)))
for (i in 1:length(col_index_list)){
  myTable_L2_by_stage[,i] <- apply(myL2_percentage[,col_index_list[[i]]], 1, mean)
}

myTable_L2_by_stage[,13] <- myL2["Phylum"]    
apply(myTable_L2_by_stage[,1:12],2,function(x) which(x == max(x))) # firmicutes dominat in 11 units except HC (Bacteroidetes most dominant)

# Merge minor phyla for visualization (phylum having < 2% in all stages have been merged into one entry)
index_to_sum_L2 <- which(apply(myTable_L2_by_stage[1:12],1,max) < 0.02)
new_sum_L2 <- apply(myTable_L2_by_stage[index_to_sum_L2, 1:12],2,sum)
plot_data_L2 <- rbind(myTable_L2_by_stage[-(index_to_sum_L2),1:12], new_sum_L2)    
stage_name_order <- paste0("stage", 1:12)
names(plot_data_L2) <- stage_name_order
plot_data_L2$Phylum <-c(myTable_L2_by_stage[-(index_to_sum_L2),"Phylum"], "Other/Unassigned")

plot_data_L2$percentage <- apply(plot_data_L2[,1:12], 1, mean)
plot_data_L2 <- plot_data_L2[order(plot_data_L2$percentage),] # ascending order
plot_data_L2 <- plot_data_L2[,1:13]

# change format for plotting
myTable_L2_by_stage_longF <- melt(plot_data_L2, id.vars = "Phylum")
myTable_L2_by_stage_longF$Phylum <- factor(myTable_L2_by_stage_longF$Phylum, levels=plot_data_L2$Phylum)

p_L2 <- 
  ggplot() + geom_bar(data=myTable_L2_by_stage_longF, aes(y=value, x=variable, fill=Phylum), colour="black", stat="identity") + 
  scale_fill_manual(values=col_vector) +
  labs(x="Management Unit ID",
       y="Relative abundance (%)",
       main="Composition at Phylum") +
  scale_x_discrete(labels=1:12) +
  # scale_y_continuous(labels = scales::percent) +
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
    legend.direction = "vertical",
    legend.key.size =  unit(1, "in"),
    legend.key.height=unit(6,"line"),    # legend space
    legend.box.margin=margin(0,0,0,50),
    plot.margin=unit(c(0.8,0.5,0.8,0.5), "cm"))    # t,r,b,l
ggsave("Fig.2B_composition_at_phylum.tiff", plot= p_L2, width=32, height=24, compression="lzw")
