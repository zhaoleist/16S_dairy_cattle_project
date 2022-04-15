######################################################################
#    Plot no. of ASVs in each sample & each ASV in no. of samples    #
######################################################################
myL0 <- as.data.frame(ifelse(my_table_2[,1:90] == 0, 0, 1))

# Fig. not in paper
no_of_samples_per_ASV <- apply(myL0,1,sum)
# hist(no_of_samples_per_ASV, xlab="ASVs in the number of samples", breaks=9, xlim=c(0,90), col="blue", border="white", cex.axis=1.2)
df_no_of_samples_per_ASV<- data.frame(x=no_of_samples_per_ASV)
g_no_of_samples_per_ASV <- ggplot(df_no_of_samples_per_ASV) + 
  geom_histogram(aes(x=x), binwidth = 2, color="black", fill="lightblue", size=1.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=32, face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0), color="black"),
    axis.text.y = element_text(size = 32, face = "plain", margin = margin(t = 0, r = 10, b = 0, l = 0), colour = "black"),
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 32, face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
    plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Frequency") + 
  xlab("ASVs in the number of samples") +
  scale_x_continuous(breaks=seq(0,90,10),limits=c(1,91)) +
  scale_y_continuous(breaks=seq(0,500,100),limits=c(0,500))
# ggsave(filename = "Fig.999_ASV_in_no_of_samples.tiff", plot=g_no_of_samples_per_ASV, width=16, height=6, compression="lzw")

# Fig. not in paper
no_of_ASV_per_sample <- apply(myL0,2,sum)
df_ASV_per_sample <- data.frame(sample=names(no_of_ASV_per_sample),
                                value=no_of_ASV_per_sample,
                                unit=stage)
df_ASV_per_sample$unit <- sub("-[A-Z]+[a-z]*","",df_ASV_per_sample$unit)
df_ASV_per_sample$sample <- factor(df_ASV_per_sample$sample, levels=production_tracks) # correct the level order of sample names
g_ASV_per_sample <- ggplot() + 
  geom_bar(data=df_ASV_per_sample, aes(x=sample, y=value, fill=unit), colour="black", stat = "identity") +
  scale_fill_manual(values=my_colors) +
  labs(x="Samples",
       y="No. of ASVs",
       main="No. of ASVs in all the samples") +
  scale_y_continuous(limits=c(0, max(df_ASV_per_sample$value) + 20), breaks=seq(0,1800,250)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
  theme(
    axis.line = element_line(colour = 'black'),
    # axis.text.x = element_text(size = 36, colour='black', face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(size = 54, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 54, face = "bold"),
    legend.title = element_text(size = 42, face="bold"),
    legend.text = element_text(size = 60, face="bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.direction = "vertical",
    legend.key.size =  unit(0.8, "in"),
    legend.key.height=unit(2,"line"),    # legend space
    legend.box.margin=margin(0,50,0,50), 
    plot.margin=unit(c(0.8,0.5,0.8,0.5), "cm")) + # t,r,b,l
  guides(fill=guide_legend(title="Mgmt.\nUnit ID"))
# ggsave("Fig.999_ASV_per_sample.tiff", plot=g_ASV_per_sample, width=48, height=18, compression="lzw")


##########################################
#    ASV Relative Abundance by sample    #
##########################################
my_otu_ra <- my_table_2
thrhld <- 0.0075 # 0.75%
index_my_otu_ra <- which(apply(my_otu_ra[,1:90],1,sum) > (MIN_DEPTH * 90 * thrhld)) # only 1 ASV has > 2% ra
other_otu_ra <- apply(my_otu_ra[-index_my_otu_ra,1:90],2,sum)
df_my_otu_ra <- as.data.frame(rbind(my_otu_ra[index_my_otu_ra,1:90], other_otu_ra))
df_my_otu_ra$genus <- c(unlist(my_otu_ra[index_my_otu_ra,96]), "Other")
df_my_otu_ra$genus <- paste0("ASV_", row.names(df_my_otu_ra),": ", df_my_otu_ra$genus)
df_my_otu_ra[length(index_my_otu_ra)+1,91] <- "Other"
df_my_otu_ra$genus <- factor(df_my_otu_ra$genus, levels=df_my_otu_ra[order(apply(df_my_otu_ra[,1:90],1,sum)),91])
df_my_otu_ra[,1:90] <- df_my_otu_ra[,1:90]/MIN_DEPTH
df_my_otu_ra_longF <- melt(df_my_otu_ra, id.vars="genus")
df_my_otu_ra_longF$stage <- rep(stage,each = 13)
thrhld <- NULL

# Fig. not in paper
g_otu_ra <- ggplot(data=df_my_otu_ra_longF) + 
  geom_bar(aes(x=variable, y=value, fill=genus), colour="black", stat="identity") +
  scale_fill_manual(values=col_vector[5:17]) +
  labs(x="Samples",
       y="Relative Abundance (%)",
       main="Composition of ASVs") +
  # scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(labels = function(x) x*100, expand=c(0,0)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",fill = NA, size = 2)) +
  theme(
    axis.line = element_line(colour = 'black'),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, colour='black', face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(size = 54, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 54, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 60, face="bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.direction = "vertical",
    legend.key.size =  unit(0.4, "in"),
    legend.key.height=unit(5.5,"line"),    # legend space
    plot.margin=unit(c(0.8,0.5,0.8,0.8), "cm")) +   # t,r,b,l
    annotate("segment", x=as.vector(unlist(lapply(col_index_list,head, n=1))) - 0.5, # "rect" cant cahnge color?
         xend=as.vector(unlist(lapply(col_index_list,tail, n=1))) + 0.5,
         y=-0.05, yend=-0.05, col=my_colors, size=50)
# ggsave("Fig.999_ASV_relative_abundance_by_sample.tiff", plot=g_otu_ra, width=48, height=18, compression="lzw")
