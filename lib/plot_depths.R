# Analysis of Sequencing Depth impact on the algorithm
library(ggplot2)

##### Loading Data #####
# Bring in varying depths counts tables (depths 500k, 1M, 2M, 5M)
depths_1 <- read.csv("data/counts_tables/counts_depth_1.00.txt", header=T, sep=',')
depths_99 <- read.csv("data/counts_tables/counts_depth_0.99.txt", header=T, sep=',')
depths_98 <- read.csv("data/counts_tables/counts_depth_0.98.txt", header=T, sep=',')
depths_1$Identity <- "1.00 identity"
depths_99$Identity <- "0.99 identity"
depths_98$Identity <- "0.98 identity"

# 20M depth results
tab_100 <- read.csv('data/counts_tables/closedlabel_interv_1.00.txt', header=T, sep=',')
tab_100$Depth <- "20M"; tab_100$Identity <- "1.00 identity"
tab_99 <- read.csv('data/counts_tables/closedlabel_interv_0.99.txt', header=T, sep=',')
tab_99$Depth <- "20M"; tab_99$Identity <- "0.99 identity"
tab_98 <- read.csv('data/counts_tables/closedlabel_interv_0.98.txt', header=T, sep=',')
tab_98$Depth <- "20M"; tab_98$Identity <- "0.98 identity"
plc_100 <- read.csv('data/counts_tables/placebo_1.00.txt', header=T, sep=",")
plc_100$Depth <- "20M"; plc_100$Identity <- "1.00 identity"
plc_99 <- read.csv('data/counts_tables/placebo_0.99.txt', header=T, sep=",")
plc_99$Depth <- "20M"; plc_99$Identity <- "0.99 identity"
plc_98 <- read.csv('data/counts_tables/placebo_0.98.txt', header=T, sep=",")
plc_98$Depth <- "20M"; plc_98$Identity <- "0.98 identity"
wk0_100 <- read.csv('data/counts_tables/counts_wk0_1.00.txt', header=T, sep=",")
wk0_100$Depth <- "20M"; wk0_100$Identity <- "1.00 identity"
wk0_99 <- read.csv('data/counts_tables/counts_wk0_0.99.txt', header=T, sep=",")
wk0_99$Depth <- "20M"; wk0_99$Identity <- "0.99 identity"
wk0_98 <- read.csv('data/counts_tables/counts_wk0_0.98.txt', header=T, sep=",")
wk0_98$Depth <- "20M"; wk0_98$Identity <- "0.98 identity"
## add them to the other depths results
depths_1 <- rbind(depths_1, tab_100, wk0_100, plc_100)
depths_99 <- rbind(depths_99, tab_99, wk0_99, plc_99)
depths_98 <- rbind(depths_98, tab_98, wk0_98, plc_98)
## omit UCFMT009BL - problems with this sample in depths
depths_1 <- depths_1[!(depths_1$Sample == "UCFMT009_BL"),]
depths_99 <- depths_99[!(depths_99$Sample == "UCFMT009_BL"),]
depths_98 <- depths_98[!(depths_98$Sample == "UCFMT009_BL"),]

# Get expected rates from self alignments
exp_100 <- read.csv('data/self_estimates/identity_1.00_estimates_self.txt', sep='\t', header=T)
exp_100$Rate <- exp_100$Count_Align / exp_100$Total_Reads
exp_99 <- read.csv('data/self_estimates/identity_0.99_estimates_self.txt', sep='\t', header=T)
exp_99$Rate <- exp_99$Count_Align / exp_99$Total_Reads
exp_98 <- read.csv('data/self_estimates/identity_0.98_estimates_self.txt', sep='\t', header=T)
exp_98$Rate <- exp_98$Count_Align / exp_98$Total_Reads
depths_1$Exp_Rate <- mean(exp_100$Rate)
depths_99$Exp_Rate <- mean(exp_99$Rate)
depths_98$Exp_Rate <- mean(exp_98$Rate)



##### Organize Tables #####
# Combine tables & create variables/factors as needed
depths <- rbind(depths_1, depths_99, depths_98)
depths$Sample <- gsub("_", "", depths$Sample)
depths$ID <- gsub("[W|B].*$", "", depths$Sample)
depths$Depth <- factor(depths$Depth, levels=c("500k", "1M", "2M", "5M", "20M"))
depths$Depth_num <- as.numeric(gsub("k|M","",depths$Depth))
depths$Depth_num[depths$Depth_num == 500] <- 0.5
depths$Total_Sum <- depths$Total_Donor + depths$Total_Pre + depths$Unmapped
depths$Corr_Engraftment <- depths$Total_Donor / (depths$Total_Sum)
depths$Week <- factor(paste0("Week ", gsub(".*BL$", "0", gsub(".*W","", depths$Sample))),
                      levels=c("Week 0", "Week 4", "Week 8", "Week 12"))
depths$Identity <- factor(depths$Identity, levels=c("1.00 identity", "0.99 identity", "0.98 identity"))


# Plot "Error" distance from 20M results (box plots per timepoint)
## note : excluding case UCFMT011 since it was only used for the case study in this depths analysis
depths <- depths[order(depths$Week),]
depths <- depths[order(depths$ID),]
depths <- depths[order(depths$Depth),]
depths <- depths[order(depths$Identity),]
tmp_truth <- c( rep(depths$Corr_Engraftment[293:365],5),
                rep(depths$Corr_Engraftment[658:730],5),
                rep(depths$Corr_Engraftment[1023:1095],5) )
depths$Dist_to_20M <- depths$Corr_Engraftment - tmp_truth
depths$Jitter <- jitter(depths$Depth_num, factor=1.2)
## Scale the X-axis for readability
depths$Scaled_Depth <- log(depths$Depth_num)
depths$Scaled_Jitter <- jitter(depths$Scaled_Depth, factor=1.2)
tiff('figures/supplemental_depths_error.tiff', units="in", width=3, height=4, res=600)
ggplot(depths[depths$ID != "UCFMT011",], aes(x=Scaled_Depth, y=Dist_to_20M, color=Identity,
                                             fill=Identity, group=interaction(Scaled_Depth,Identity))) +
  geom_hline(aes(yintercept=0), col="grey", lwd=0.6) +
  geom_path(aes(x=Scaled_Jitter, group=Sample), alpha=0.1) +
  geom_boxplot(outlier.shape=NA, alpha=0.2, lwd=0.4) +
  geom_point(aes(x=Scaled_Jitter), alpha=0.2, cex=0.4) +
  scale_color_manual(values=rev(c("purple","green3","blue"))) +
  scale_fill_manual(values=rev(c("purple","green3","blue"))) +
  labs(title="Sequencing Depth Sensitivity",
       x="Sequencing Depth (post samples)",
       y="Distance to Engraftment at full depth") +
  scale_x_continuous(breaks=unique(depths$Scaled_Depth), labels=c("500k","1M","2M","5M","20M")) +
  facet_wrap(~Identity, nrow=3) +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        strip.text.x = element_text(size=7, family="Helvetica"))
dev.off()


# Varying Depths of Databases too : UCFMTDepth_num# Varying Depths of Databases too : UCFMT011 case study
# for some reason the UCFMT011 samples have the greatest variance in the depths at weeks 4-12
tmp_plt <- depths[depths$ID == "UCFMT011",]
ggplot(tmp_plt, aes(x=Depth_num, y=Dist_to_20M, color=Identity,
                    group=interaction(Depth_num,Identity), shape=Week)) +
  geom_hline(aes(yintercept=0), col="grey", lwd=0.75) +
  geom_path(aes(x=Jitter, group=Sample), alpha=0.6) +
  geom_point(aes(x=Jitter), alpha=0.6, cex=2) +
  scale_color_manual(values=rev(c("purple","green3","blue"))) +
  scale_shape_manual(values=c(16, 17, 18, 15)) +
  labs(title="Case Study : UCFMT011",
       x="Sequencing Depth : post samples (millions)",
       y="Distance to Engraftment at full depth") +
  ylim(-0.007, 0.001) +
  facet_wrap(~Identity, nrow=3) +
  theme_bw() +
  theme(title = element_text(face="bold", hjust=0.5))
# lets look at the varying depths by all samples for UCFMT011
uc011 <- read.csv('data/counts_tables/depth011_0.98.txt', header=T, sep=",")
uc011 <- rbind(uc011, read.csv('data/counts_tables/depth011_0.99.txt', header=T, sep=","))
uc011 <- rbind(uc011, read.csv('data/counts_tables/depth011_1.00.txt', header=T, sep=","))
uc011$Identity <- factor(rep(c("0.98 identity","0.99 identity","1.00 identity"), each=15),
                         levels=c("1.00 identity","0.99 identity","0.98 identity"))
uc011$Depth_num <- as.numeric(gsub("[k|M]","",uc011$Depth))
uc011$Depth_num[uc011$Depth_num == 500] <- 0.5
uc011$Week <- factor(paste0("Week ", gsub("UCFMT.*[W|BL]","",uc011$Sample)),
                     levels=c("Week 4", "Week 8", "Week 12"))
uc011 <- uc011[order(uc011$Week),]
uc011 <- uc011[order(uc011$Depth_num),]
uc011 <- uc011[order(uc011$Identity),]
rownames(uc011) <- 1:nrow(uc011)
uc011$Total_Sum <- uc011$Total_Donor + uc011$Total_Pre + uc011$Unmapped
uc011$Perc_Unmapped <- uc011$Unmapped / uc011$Total_Sum
uc011$Corr_Engraftment <- uc011$Total_Donor / (uc011$Total_Sum * 0.9212651)
tmp_truth011 <- c(rep(uc011$Corr_Engraftment[13:15], 5),
                  rep(uc011$Corr_Engraftment[28:30], 5),
                  rep(uc011$Corr_Engraftment[43:45], 5))
uc011$Dist_to_20M <- uc011$Corr_Engraftment - tmp_truth011
uc011$Scaled_Depth <- log(uc011$Depth_num)
tiff('figures/supplemental_case_depths.tiff', units="in", width=3, height=4, res=600)
ggplot(uc011[rev(1:nrow(uc011)),], aes(x=Scaled_Depth, y=Dist_to_20M, color=Identity,
                  group=interaction(Scaled_Depth,Identity), shape=Week)) +
  geom_hline(aes(yintercept=0), col="grey", lwd=0.6) +
  geom_path(aes(group=interaction(Sample, Identity)), alpha=0.5, show.legend=F) +
  geom_point(alpha=0.6, cex=2) +
  scale_color_manual(values=rev(c("purple","green3","blue"))) +
  scale_shape_manual(values=c(17, 18, 15)) +
  scale_x_continuous(breaks=unique(uc011$Scaled_Depth), labels=c("500k","1M","2M","5M","20M")) +
  labs(title="Case Study: Database & Sample Depths",
       x="Sequencing Depth",
       y="Distance to Engraftment at full depth", shape="") +
  guides(color="none") +
  facet_wrap(~Identity, nrow=3) +
  theme_bw() +
  theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        strip.text.x = element_text(size=7, family="Helvetica"),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,-10,-10),
        legend.text = element_text(size=6, family="Helvetica"))
dev.off()


