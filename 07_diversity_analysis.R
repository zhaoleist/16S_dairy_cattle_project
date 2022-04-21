###################################
#    Alpha diversity: use         #
#       --> chao 1: richness      #
#       --> shannon: diversity    #
###################################

# shannon diversity index
my_table_2_run <- my_table_2
H <- diversity(t(my_table_2_run[,1:90]))    # from package vegan

shapiro.test(H)    # normality test
stage_for_H <- stage
# allPair <- do.call(c, apply(combn(1:12, 2), 2, list))
# shannon(my_table_2_run[,1])    from OTUtable package

# chao1 diversity index
myChao1 <- function(x){
  Sobs <- length(which(x != 0))
  f1 <- length(which(x == 1))
  f2 <- length(which(x == 2))
  Schao1 <- Sobs + f1 * (f1 - 1)/(2 * (f2 + 1))
  return(Schao1)
}
chao1_values <- apply(my_table_2_run[,1:90], 2, myChao1)

shapiro.test(chao1_values) 
stage_for_chao1_values <- stage

shapiro.test(no_of_ASV_per_sample)
meta <- data.frame(sample=names(chao1_values), 
                   stage=stage_for_chao1_values,
                   chao1=chao1_values,
                   shannon=H,
                   OTU=no_of_ASV_per_sample) # stage=gsub("\\.[0-9]{3}","",names(H)))
meta$stage <- factor(meta$stage, levels=unique(meta$stage))    # factor here, so split can preserve natural stage order

Z <- split(meta, meta$stage)
for(i in Z) print(shapiro.test(i$chao1)$p.value)
for(i in Z) print(shapiro.test(i$shannon)$p.value)

leveneTest(meta$chao1,meta$stage) # or leveneTest(chao1 ~ stage, meta); to test if variance equal (p-value > 0.05 : equal variance)
leveneTest(meta$shannon, meta$stage)  
bartlett.test(meta$chao1,meta$stage)
bartlett.test(meta$shannon,meta$stage)
# difference between the two tests
# check http://atomic.phys.uni-sofia.bg/local/nist-e-handbook/e-handbook/eda/section3/eda35a.htm


# equal variance but not from normally distributed population --> not ANOVA

# aov_chao1 <- aov(chao1 ~ stage, data=meta)
# summary(aov_chao1)
# TukeyHSD(aov_chao1)
# 
# text_for_chao <- c("3~12","3,7,\n11,12","1,2","1","1","1","1,2,8","1,7","1","1","1,2","1,2")
# par(mfrow=c(2,1), mai=c(1.02,0.82,0.82,0.42))
# plot(chao1 ~ stage, data=meta, xaxt='n',  xlab="Management Unit ID", ylab="Chao 1 richness index", cex.lab=1.5, col=myColors)
# stripchart(chao1 ~ stage, data = meta, vertical = TRUE, pch = 19, add=TRUE, col="red")
# axis(side=1, at=1:12, labels=1:12)
# # text(x=1:12,y=2700,label=text_for_chao, col=myColors )
# mtext(text_for_chao, font=2, cex=1.4, col=myColors, side=3, at=1:12, outer = FALSE)

# gChao1 <- ggplot(meta, aes(x=stage, y=chao1)) + 
#   geom_boxplot() +
#   stat_boxplot(geom ='errorbar') +
#   geom_dotplot(binaxis = "y",
#                stackdir = "center",
#                dotsize = 20,
#                fill = "red",
#                binwidth = 1) +
#   # geom_signif(comparisons = list(),
#   #             map_signif_level = TRUE, test = "t.test") + 
#   # annotate("text", label="3-12", x=1, y=500) +
#   scale_x_discrete(labels=myLevels) + 
#   labs(x="Stage", y="Chao Richness Index") +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black", size=0.6))
# gChao1

# plot alpha diversity
tiff(filename = "Fig.3_alpha_diversity.tiff", res=300, width=16, height=12, units = "in", compression="lzw")
par(mfrow=c(2,1), mar=c(7,10,4,2) + 0.1)    # b,l,t,r
kruskal.test(chao1 ~ stage, data=meta)
pairwise.wilcox.test(meta$chao1, meta$stage, p.adjust.method="BH")
text_for_chao <- c("3~12","3,4,\n9,11","1,2,7","1,2,7","1,10","1,10",
                   "1,3,\n4,10","1","1,2","1,5,\n6,7","1,2","1")

plot(chao1 ~ stage, data=meta, xaxt='n', yaxt='n', xlab="", ylab="",  
     cex.axis=2, cex.lab=2, col=my_colors, ylim=c(350,2750))
stripchart(chao1 ~ stage, data = meta, vertical = TRUE, pch = 19, add=TRUE, col="violetred")
axis(side=1, at=1:12, labels=1:12, cex.axis=2.4)
axis(side=2, cex.axis=2.5, las=2)
mtext(side=2, line=6, "Chao 1 Richness Index", col="black", cex=2.4) # adjust y-axis
mtext(text_for_chao, font=2.4, cex=2, col=my_colors, side=3, at=1:12, outer = FALSE)

kruskal.test(shannon ~ stage, data=meta)
pairwise.wilcox.test(meta$shannon, meta$stage, p.adjust.method="BH")
text_for_shannon <- c("3~12","3,4,9","1,2,\n5,7","1,2,5,\n6,7", "1,3,\n4,10","1,4,10",
                      "1,3,\n4,10", "1", "1,2,10", "1,5,6,\n7,9", "1", "1")
plot(shannon ~ stage, data=meta, xaxt='n', yaxt='n', xlab="", ylab="", 
     cex.axis=3.0, cex.lab=3.0, col=my_colors)
stripchart(shannon ~ stage, data = meta, vertical = TRUE, pch = 19, add=TRUE, col="violetred")
axis(side=1, at=1:12, labels=1:12, cex.axis=2.4)
axis(side=2, cex.axis=2.4, las=2)
mtext(side=1, line=4, "Manegement Unit ID", cex=2.4)
mtext(side=2, line=6, "Shannon Diversity Index", col="black", cex=2.4) # adjust y-axis
mtext(text_for_shannon, font=2.4, cex=2,col=my_colors, side=3, at=1:12, outer = FALSE)
# save as 1200 x 800 pixels
dev.off() # reset par() to default by calling dev.off()

kruskal.test(OTU ~ stage, data=meta)
pairwise.wilcox.test(meta$OTU, meta$stage, p.adjust.method = "BH")


########################
#    Beta diversity    #
########################
my_table_2_beta <- as.data.frame(t(my_table_2))
names(my_table_2_beta) <- paste0("ASV_",1:ncol(my_table_2_beta))
my_table_2_beta <- my_table_2_beta[-(91:96),]
my_table_2_beta <- sapply(my_table_2_beta, function(x) as.numeric(as.character(x)))

# Bray Curtis
BC.nmds <- metaMDS(my_table_2_beta, distance = "bray", k=2, trymax=1000)
# par(mfrow=c(1,1),mar=c(5,5,5,5))
# plot(BC.nmds, type="n", main="Bray-Curtis")
# points(BC.nmds, display="sites", pch=20, col=myColors[meta$stage], cex=3.5)
# legend(-2.5,-1.5, ncol=2, legend= factor(unique(stage_for_H)), col=myColors, 
#        pch=20, cex = 0.5, pt.cex = 1.5)

plot(scores(BC.nmds, display="site")[,1],scores(BC.nmds, display="site")[,2])
df_nmds <- data.frame(NMDS1=scores(BC.nmds)[,1], NMDS2=scores(BC.nmds)[,2], Stage=stage)
df_nmds$Stage <- sub("-[A-Z]+[a-z]*","",df_nmds$Stage)
# df_nmds$Stage <- factor(df_nmds$Stage, levels=unique(df_nmds$Stage))
df_nmds$Stage <- as.numeric(df_nmds$Stage)
df_nmds$Stage <- factor(df_nmds$Stage, levels=1:12)
df_nmds$sample <- factor(meta$sample)

gNMDS <- 
  ggplot(df_nmds) +
  geom_point(aes(x=NMDS1, y=NMDS2,  col=Stage, fill=Stage, shape=Stage), size=8, stroke=1.6) +
  # geom_label_repel(aes(x=NMDS1, y=NMDS2, label=sample)) +
  scale_color_manual(labels=unique(df_nmds$Stage),values=c(rep("black",11),"red")) +
  scale_fill_manual(labels=unique(df_nmds$Stage),values=my_colors) +
  scale_shape_manual(labels=unique(df_nmds$Stage),values=c(rep(22, each=4),
                                                   rep(24, each=2),
                                                   rep(21, each=5),
                                                   rep(8, each=1))) +
  # labs(x="sample",
  #      y="relative abundance",
  #      main= text) +
  coord_cartesian(xlim = c(-4, 2)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, size = 1.2)) +
  theme(
    # axis.line = element_line(colour = 'black', size=1),
    axis.text.x = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0), colour="black"),
    axis.text.y = element_text(size = 24, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0), colour="black"),
    axis.title.x = element_text(size = 24, face = "bold", colour="black"),
    axis.title.y = element_text(size = 24, face = "bold", colour="black"),
    legend.title = element_text(size= 24, face="bold"),
    legend.text = element_text(size = 24, face="bold"),
    # legend.justification=c(1,0),
    # legend.position= c(0.90,0.82),
    # legend.position= c(0.20,0.15),
    legend.position= "right",
    # legend.position="None",
    legend.box.background = element_rect(colour="white"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.direction = "vertical",
    legend.key.size =  unit(0.4, "in")) +
  # annotate("text", x=1.73, y=0.9, label="Management Unit ID", size=5.5, fontface="bold") +
  guides(fill=guide_legend(nrow=12,byrow=TRUE, title="Mgmt. \nUnit ID"),
         col=guide_legend(nrow=12,byrow=TRUE, title="Mgmt. \nUnit ID"),
         shape=guide_legend(nrow=12,byrow=TRUE, title="Mgmt. \nUnit ID"))
gNMDS

# lower claster
sampleID_lower_cluster <- c("BH.008", "CUC.004","CUC.005","Hos.001","MLLH.003",
                            "FC.006","FH.006","FC.005","FC.003","MLC.007",
                            "FC.004","MLC.005","Hos.005", "FC.002")
index_lower_cluster <- c()
for (iter in 1:length(sampleID_lower_cluster)){
  index_lower_cluster[iter] <-  which(df_nmds$sample == sampleID_lower_cluster[iter]) 
}
df_nmds_lower <- df_nmds[index_lower_cluster,]

# left cluster
sampleID_left_cluster <- c("PWH.003","HC.002","HC.009","HC.004","HC.001","HC.003",
                           "HC.006","HC.007", "HC.008", "PWH.001","HC.005")
index_left_cluster <- c()
for (iter in 1:length(sampleID_left_cluster)){
  index_left_cluster[iter] <-  which(df_nmds$sample == sampleID_left_cluster[iter]) 
}
df_nmds_left <- df_nmds[index_left_cluster,]

# df_nmds_upper <- df_nmds[-c(index_lower_cluster, index_left_cluster),]
# upper cluster
index_upper_cluster <- setdiff(1:90, c(index_lower_cluster,index_left_cluster)) 
sampleID_upper_cluster <- setdiff(names(myL2[1:90]), c(sampleID_lower_cluster,sampleID_left_cluster))
df_nmds_upper <- df_nmds[index_upper_cluster,]

gNMDS_2 <- gNMDS + 
  geom_encircle(aes(x=NMDS1, y=NMDS2), 
                data=df_nmds_lower, 
                color="red", 
                size=6, 
                expand=0.045) +
  geom_encircle(aes(x=NMDS1, y=NMDS2), 
                data=df_nmds_upper, 
                color="black", 
                size=6, 
                expand=0.04) +
  geom_encircle(aes(x=NMDS1, y=NMDS2), 
                data=df_nmds_left, 
                color="blue", 
                size=6, 
                expand=0.04) +
  annotate(geom="text", x=-3, y=0.35, label=as.roman(1), color="blue", size=16) +
  annotate(geom="text", x=1.05, y=0.5, label=as.roman(2), color="black", size=16) +
  annotate(geom="text", x=1.35, y=-1, label=as.roman(3), color="red", size=16)
gNMDS_2
ggsave("Fig.4.tiff", plot=gNMDS_2, width=20, height=12, compression="lzw")

# encircle based on stage
gNMDS + geom_encircle(aes(x=NMDS1, y=NMDS2, col=Stage), 
                data=df_nmds, 
                size=3, 
                expand=0.04) 
  
# # Jaccard
# J.nmds <- metaMDS(my_table_2_beta, distance = "jaccard", k=2, trymax=1000)
# par(mar=c(2,2,2,4), new=TRUE)
# plot(J.nmds, type="n", main="Jarcard")
# points(J.nmds, display="sites", pch=20, col=my_colors, cex=2.5)
# legend("topleft", legend= factor(unique(stage)), col=my_colors, 
#        pch=20, cex = 0.6, pt.cex = 2,text.width = 2, inset=c(1,0),xpd=TRUE)


# PERMANOVA
# BC.dist <- vegdist(my_table_2_beta, distance="bray")
# adonis(BC.dist ~ stage, data=meta, permutations = 1000)
# beta dispersion
# disp.stage <- betadisper(BC.dist, meta$stage)
# permutest(disp.stage, pairwise=TRUE, permutations=1000)
# simper(my_table_2_beta,meta$stage,  permutations=5)



# Patescibacteria
# df_Patescibacteria <- data.frame(count=unlist(c(myL2[1,sampleID_lower_cluster], myL2[1,sampleID_upper_cluster])),
#                                  cluster=c(rep("lower", length(sampleID_lower_cluster)), rep("upper", length(sampleID_upper_cluster))))


##############################################
#    Explore two clusters at phylum level    #
##############################################
df_lower <- myL2[,c(sampleID_lower_cluster,"Phylum")]
# df_upper <- myL2[,setdiff(production_tracks[c(1:90,92)],c(sampleID_lower_cluster,sampleID_left_cluster))]
df_upper <- myL2[,c(sampleID_upper_cluster,"Phylum")]

df_lower_longF <- melt(df_lower, id.vars="Phylum")
ggplot(df_lower_longF) + geom_bar(aes(x=variable,y=value, fill=Phylum), stat="identity") +
  scale_fill_manual(values=col_vector) 

df_upper_longF <- melt(df_upper, id.vars="Phylum")
ggplot(df_upper_longF) + geom_bar(aes(x=variable,y=value, fill=Phylum), stat="identity") +
  scale_fill_manual(values=col_vector)
 
# Fig.4B & S3
for (i in 1:nrow(myL2)){
  df_explore_two_clusters <- data.frame(t((myL2[i,1:90]/MIN_DEPTH) *100))
  df_explore_two_clusters$sampleID <- row.names(df_explore_two_clusters)
  df_upper_pates <- df_explore_two_clusters[df_explore_two_clusters$sampleID %in% sampleID_upper_cluster,]
  df_lower_pates <- df_explore_two_clusters[df_explore_two_clusters$sampleID %in% sampleID_lower_cluster,]
  # median(df_lower_pates$X1)
  # mean(df_lower_pates$X1)
  # median(df_upper_pates$X1)
  # mean(df_upper_pates$X1)
  df_upper_pates$cluster <- rep(paste0("cluster ",as.roman(2)), nrow(df_upper_pates))
  df_lower_pates$cluster <- rep(paste0("cluster ",as.roman(3)), nrow(df_lower_pates))
  df_for_pates <- as.data.frame(rbind(df_upper_pates,df_lower_pates))
  names(df_for_pates)[1] <- "X1"
  my_plot <- ggboxplot(df_for_pates, x = "cluster", y = "X1",
               color = "cluster",
               add = "jitter", lwd=1.35)
# Add p-value
# my_plot <- my_plot + stat_compare_means()
# Change method
  my_plot <- my_plot + stat_compare_means(method = "t.test", size=12, label.x=0.75, label.y=max(df_for_pates$X1)) + 
    labs(x="cluster in Fig.4A", y="Relative Abundance (%)", color=paste0(myL2[i,"Phylum"],"      ","\n\n\n")) +
    scale_colour_manual(values=c("black", "red")) + theme_classic(base_size=52)
  ggsave(filename = paste0("Fig.S3_",myL2[i,"Phylum"],".png"),plot = my_plot, width = 16, height = 12, dpi = 300)
}
