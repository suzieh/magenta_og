# Looking at the output from the Ulcerative colitis dataset,
# let's compare the placebo and intervention groups.
library(ggplot2)
library(ggpubr)
library(dplyr)

# Loading closed label data
tab_100 <- read.csv('data/counts_tables/closedlabel_interv_1.00.txt', header=T, sep=',')
tab_99 <- read.csv('data/counts_tables/closedlabel_interv_0.99.txt', header=T, sep=',')
tab_98 <- read.csv('data/counts_tables/closedlabel_interv_0.98.txt', header=T, sep=',')

# Laoding placebo data
plc_100 <- read.csv('data/counts_tables/placebo_1.00.txt', header=T, sep=",")
plc_99 <- read.csv('data/counts_tables/placebo_0.99.txt', header=T, sep=",")
plc_98 <- read.csv('data/counts_tables/placebo_0.98.txt', header=T, sep=",")

# Loading wk0 data from both placebo & donor
wk0_100 <- read.csv('data/counts_tables/counts_wk0_1.00.txt', header=T, sep=",")
wk0_99 <- read.csv('data/counts_tables/counts_wk0_0.99.txt', header=T, sep=",")
wk0_98 <- read.csv('data/counts_tables/counts_wk0_0.98.txt', header=T, sep=",")

# Loading expected rates from self alignments
exp_100 <- read.csv('data/self_estimates/identity_1.00_estimates_self.txt', sep='\t', header=T)
exp_100$Rate <- exp_100$Count_Align / exp_100$Total_Reads
exp_99 <- read.csv('data/self_estimates/identity_0.99_estimates_self.txt', sep='\t', header=T)
exp_99$Rate <- exp_99$Count_Align / exp_99$Total_Reads
exp_98 <- read.csv('data/self_estimates/identity_0.98_estimates_self.txt', sep='\t', header=T)
exp_98$Rate <- exp_98$Count_Align / exp_98$Total_Reads

# ----- Choose preferred identity here -----
df <- rbind(tab_100, plc_100, wk0_100)
df$Identity <- "1.00 identity"
exp_rate <- max(exp_100$Rate)

# Set up variables for plotting
df$Type <- c( rep("MTT", nrow(tab_98)), rep("placebo", nrow(plc_98)),
              rep("MTT", 10), rep("placebo", 12))
df <- df[!grepl("OL", df$Sample),]
df$Sample <- gsub("W","_W", df$Sample)
df$Sample <- gsub("BL","_BL", df$Sample)
df$Sample <- gsub("__","_", df$Sample)
df$ID <- gsub("_.*", "", df$Sample)
df$Week <- as.numeric(gsub("W","",gsub("UCFMT[0-9][0-9][0-9]_","",df$Sample)))
df$Week[is.na(df$Week)] <- 0
df <- df[order(df$Week),]
df$Total_Sum <- df$Total_Donor + df$Total_Pre + df$Unmapped
df$Expected_Rate <- exp_rate
df$Raw_Engraftment <- df$Total_Donor / df$Total_Sum
df$Corr_Engraftment <- df$Total_Donor / (df$Total_Sum * df$Expected_Rate)

# Filter out incomplete courses (UCFMT001)
df <- df[df$ID != "UCFMT001",]


# Get Donor ID per sample (optional for plot)
cl_d <- read.csv("data/comparison_tables/post_table.txt", sep="\t", header=T)
pl_d <- read.csv("data/comparison_tables/placebo_post_table.txt", sep="\t", header=T)
df$Donor_ID <- c(cl_d$Donor_Set, pl_d$Donor_Set)[match(df$ID, c(cl_d$Post_ID, pl_d$Post_ID))]
df$Label_ID <- df$ID
df$Label_ID[which(df$Week != 8)] <- NA
df$Label_X <- NA
df$Label_X[which(df$Week == 8)] <- 8 + runif(sum(df$Week == 8), min=-2, max=2)


# Compare by true donor & placebo:donor pairing
translate <- data.frame(Original=c("71","70a","70b","70c","70z","102a","102b","102c","102d","102e","102f","102g"),
                        New=c("Donor B","Donor A1","Donor A2","Donor A3","Donor A1 ",
                              "Donor C3","Donor C1","Donor C-","Donor C4","Donor C-", "Donor C2", "Donor C5"))
df$Donor_ID_Uniq <- df$Donor_ID
df$Donor_ID_Uniq[df$ID %in% c("UCFMT009","UCFMT014")] <- "70z"
df$Donor_ID_Uniq <- translate$New[match(df$Donor_ID_Uniq, translate$Original)]
plot_df <- df[!(df$Donor_ID %in% c("102c","102e")),]
rownames(plot_df) <- 1:nrow(plot_df)
plot_df$Corr_Engraftment <- plot_df$Corr_Engraftment * 100
p <- ggplot(plot_df[nrow(plot_df):1,], aes(x=Week, y=Corr_Engraftment)) +
  geom_line(aes(group=ID, col=Type), alpha=0.7, lwd=1) +
  geom_point(aes(col=Type), alpha=0.7, size=2) +
  scale_color_manual(values=c("blue","pink2")) + 
  scale_x_continuous(breaks=c(0,4,8,12), limits=c(0,15)) +
  facet_wrap(~Donor_ID_Uniq, nrow=1) +
  labs(title="MTT vs Placebo Engraftment per Donor", x="Week", y="% Engraftment", col="") +
  theme_bw()
p + theme(title=element_text(face="bold"))
p + geom_text(aes(x=Label_X, label=Label_ID, col=Type), vjust=-2, show.legend=F)
## save figure
tiff("figures/trial_engraft_lines.tiff", units="in", width=6, height=2.5, res=600)
p + theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
          axis.title = element_text(face="bold", size=7, family="Helvetica"),
          axis.text = element_text(size=6, family="Helvetica"),
          strip.text.x = element_text(size=6, family="Helvetica"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.text = element_text(size=7, family="Helvetica"))
dev.off()


# Plot all lines together, include means
## calculate group means
grp_means <- data.frame(Week=rep(c(0,4,8,12),2),
                        Type=rep(c("MTT", "placebo"), each=4))
grp_means$Corr_Engraftment <- NA
for (i in 1:nrow(grp_means)) {
  grp_means$Corr_Engraftment[i] <- mean(df$Corr_Engraftment[df$Week == grp_means$Week[i] &
                                                              df$Type == grp_means$Type[i]])
}
## plot all engraftment lines & mean lines
tiff("figures/supplemental_engraft_mean_lines.tiff", units="in", width=3, height=3, res=600)
ggplot(df, aes(x=Week, y=Corr_Engraftment)) + 
  geom_line(aes(group=ID, col=Type), alpha=0.3, lwd=0.6, show.legend=F) +
  geom_line(data=grp_means, mapping=aes(x=Week, y=Corr_Engraftment, col=Type, group=Type), lwd=1, alpha=0.7, lty="dashed") +
  scale_color_manual(values=c("blue","deeppink")) + 
  scale_x_continuous(breaks=c(0,4,8,12), limits=c(0,12)) +
  labs(title="MTT vs. Placebo Engraftment", x="Week", y="% Engraftment", col="Group") +
  theme_bw() +
  theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        legend.text = element_text(size=6, family="Helvetica"),
        legend.title = element_text(size=7, face="bold", family="Helvetica"),
        #legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(2,2,2,2),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0.4, "cm"),
        legend.position = c(0.16,0.85),
        legend.box.background = element_rect(colour = "grey"))
dev.off()


# Intervention vs. Placebo Boxplots
# df$Dist_Grp_Mean <- NA
# for (i in 1:nrow(df)) {
#   df$Dist_Grp_Mean[i] <- abs(df$Corr_Engraftment[i] - 
#                                grp_means$Corr_Engraftment[grp_means$Week == df$Week[i] & 
#                                                             grp_means$Type == df$Type[i]])
# }
df$JitteredWk <- df$Week + ifelse(df$Type == "MTT", -0.75, 0.75)
typep <- df %>% group_by(Week) %>%
  summarise(wilcox.p = wilcox.test(Corr_Engraftment ~ Type, exact = FALSE)$p.value)
typep$Type="MTT"
df$Corr_Engraftment <- df$Corr_Engraftment * 100
tiff("figures/engraftment_boxplots.tiff", units="in", width=3.5, height=2.5, res=600)
ggplot(df, aes(x=Week, y=Corr_Engraftment, col=Type, fill=Type,
               group=interaction(Type, Week))) +
  geom_boxplot(alpha=0.4, outlier.shape=NA) + 
  geom_jitter(aes(x=JitteredWk), size=1, width=0.25, alpha=0.7) +
  scale_color_manual(values=c("blue","pink2")) +
  scale_fill_manual(values=c("blue","pink2")) +
  scale_x_continuous(breaks=c(0,4,8,12), limits=c(-2,14)) +
  geom_text(data=typep, aes(x=Week, y=65, label=paste0("p=",round(wilcox.p,3),ifelse(wilcox.p<0.05,"*",""))),
            col="grey8", size=6/.pt, family="Helvetica") +
  labs(title="Engraftment per Week", x="Week", y="% Engraftment", col="", fill="") +
  theme_bw() +
  theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        legend.text = element_text(size=6, family="Helvetica"))
dev.off()




