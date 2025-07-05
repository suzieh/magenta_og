# Given an outcome table of the following format:
## Simulation | Uniquely Match Donor | Uniquely Match Pre | Ambiguous | Total Donor | Total Pre | Unknown/Unique |
# We want to make a graph showing % donor simulated vs. % donor found
## note that "Simulation" is % Pre sample
library(scales)

# To determine corrected engraftment values, we adjust by the amount of alignment
#   between donor & pre samples with themselves (their own MAGs databases) at
#   each identity.
self_1 <- read.csv('data/self_estimates/simengraft_1.00.txt', sep='\t', header=T)
self_1$Rate <- self_1$Count_Align / self_1$Total_Reads
self_1 <- max(self_1$Rate)     # 0.751
self_99 <- read.csv('data/self_estimates/simengraft_0.99.txt', sep='\t', header=T)
self_99$Rate <- self_99$Count_Align / self_99$Total_Reads
self_99 <- max(self_99$Rate)   # 0.807
self_98 <- read.csv('data/self_estimates/simengraft_0.98.txt', sep='\t', header=T)
self_98$Rate <- self_98$Count_Align / self_98$Total_Reads
self_98 <- max(self_98$Rate)   # 0.827

# Collect and organize tables
## Loading tables
tab_100 = read.csv('data/counts_tables/simengraft_identity_1.00.txt', header=T, sep=',')
tab_99 = read.csv('data/counts_tables/simengraft_identity_0.99.txt', header=T, sep=',')
tab_98 = read.csv('data/counts_tables/simengraft_identity_0.98.txt', header=T, sep=',')
## Combine tables, add identity value
tab_100$Identity <- "1.00"
tab_99$Identity <- "0.99"
tab_98$Identity <- "0.98"
sim_tab <- rbind(tab_100, tab_99, tab_98)
## Clean Up
sim_tab$Simulation_Donor <- 1.0 - sim_tab$Simulation               # Add Simulation_Donor
sim_tab$Total_Sum <- sim_tab$Total_Donor + sim_tab$Total_Pre + sim_tab$Unmapped  # Get total reads
sim_tab$Engrafted_Raw <- sim_tab$Total_Donor / sim_tab$Total_Sum   # Raw Engraftment (including all unmapped reads)
sim_tab$Expected_Align_Rate <- rep(c(self_1, self_99, self_98), each=nrow(tab_100)) # Best expected alignment rate from self alignments
sim_tab$Engrafted_Corr <- (sim_tab$Total_Donor) / (sim_tab$Total_Sum * sim_tab$Expected_Align_Rate) # Corrected Engraftment by expected alignment rate
sim_tab <- sim_tab[rev(1:nrow(sim_tab)),] # reverse order prior to plotting

# Additional simulations (to confirm results)
## loading additional simulations
add_100 <- read.csv('data/counts_tables/counts_simulations_1.00.txt', header=T, sep=',')
add_99 <- read.csv('data/counts_tables/counts_simulations_0.99.txt', header=T, sep=',')
add_98 <- read.csv('data/counts_tables/counts_simulations_0.98.txt', header=T, sep=',')
add_100$Identity <- "1.00"
add_99$Identity <- "0.99"
add_98$Identity <- "0.98"
sim_add <- rbind(add_100, add_99, add_98)
sim_add$Simulation_Donor <- sim_add$Simulation_Donor / 100
sim_add$Total_Sum <- 8750000
### get other self-alignments
sim_add$Expected_Align_Rate <- 1.0 # default value of 100% alignment
for (sub in unique(sim_add$Subject)) {
  # get pre/donor self alignment percentages - results are per identity
  pre <- sim_add[sim_add$Subject == sub & sim_add$Simulation_Donor == 0, "Total_Pre"] / 8750000
  donor <- sim_add[sim_add$Subject == sub & sim_add$Simulation_Donor == 1, "Total_Donor"] / 8750000
  # set average of the donor/pre self alignments as expected alignment rate
  sim_add$Expected_Align_Rate[sim_add$Subject == sub] <- c(rep(mean(c(pre[1], donor[1])),11),
                                                           rep(mean(c(pre[2], donor[2])),11),
                                                           rep(mean(c(pre[3], donor[3])),11) )
}
sim_add$Expected_Align_Rate <- round(sim_add$Expected_Align_Rate, 5)
table(sim_add$Expected_Align_Rate)
sim_add$Engrafted_Corr <- (sim_add$Total_Donor) / (sim_add$Total_Sum * sim_add$Expected_Align_Rate)


# Plot results for all identities, single simulation:
plot(seq(0,1,0.1), seq(0,1,0.1), col=alpha("black", 0.9), pch=18, cex=2, type="b", lwd=2,
     main="Simulated Engraftment", xlab="% Donor Simulated", ylab="% Donor Engraftment")
lines(seq(0,1,0.1), sim_tab$Engrafted_Raw[sim_tab$Identity == "1.00"], col=alpha("blue",0.5), pch=17, cex=1, type="b", lty=2, lwd=2)
lines(seq(0,1,0.1), sim_tab$Engrafted_Corr[sim_tab$Identity == "1.00"]-0.02, col=alpha("blue",0.7), pch=19, cex=1.5, type="b", lty=2, lwd=3)
lines(seq(0,1,0.1), sim_tab$Engrafted_Raw[sim_tab$Identity == "0.99"], col=alpha("green3",0.5), pch=17, cex=1, type="b", lty=2, lwd=2)
lines(seq(0,1,0.1), sim_tab$Engrafted_Corr[sim_tab$Identity == "0.99"]-0.01, col=alpha("green3",0.7), pch=19, cex=1.5, type="b", lty=2, lwd=3)
lines(seq(0,1,0.1), sim_tab$Engrafted_Raw[sim_tab$Identity == "0.98"], col=alpha("purple",0.5), pch=17, cex=1, type="b", lty=2, lwd=2)
lines(seq(0,1,0.1), sim_tab$Engrafted_Corr[sim_tab$Identity == "0.98"], col=alpha("purple",0.7), pch=19, cex=1.5, type="b", lty=2, lwd=3)
legend(0.05, 1.0, legend=c("Expected", "Observed", "Corrected", "1.00", "0.99", "0.98"),
       col=c("black","darkgrey","grey","blue","green3","purple"), lty = c(1,2,2,1,1,1), pch=c(18,17,19,19,19,19), cex=0.8, lwd=2)


# Plot results for all identities (corrected values only, single simulation)
plot(seq(0,1,0.1), seq(0,1,0.1), col=alpha("black", 0.9), pch=17, cex=2, type="b", lwd=2,
     main="Simulated Engraftment", xlab="% Donor Simulated", ylab="% Donor Engraftment")
lines(seq(0,1,0.1), sim_tab$Engrafted_Corr[sim_tab$Identity == "1.00"]-0.02, col=alpha("blue",0.7), pch=19, cex=1.5, type="b", lty=2, lwd=3)
lines(seq(0,1,0.1), sim_tab$Engrafted_Corr[sim_tab$Identity == "0.99"]-0.01, col=alpha("green3",0.7), pch=19, cex=1.5, type="b", lty=2, lwd=3)
lines(seq(0,1,0.1), sim_tab$Engrafted_Corr[sim_tab$Identity == "0.98"], col=alpha("purple",0.7), pch=19, cex=1.5, type="b", lty=2, lwd=3)
legend(0.05, 1.0, legend=c("Expected", "1.00", "0.99", "0.98"),
       col=c("black","blue","green3","purple"), lty = c(1,1,1,1), pch=c(17,19,19,19), cex=0.8, lwd=2)


# Error plot per identity (corrected values)
sim_tab$Error <- sim_tab$Engrafted_Corr - sim_tab$Simulation_Donor
sim_tab$Error_Raw <- sim_tab$Engrafted_Raw - sim_tab$Simulation_Donor
plot(x=seq(0,1,0.1), y=sim_tab$Error_Raw[sim_tab$Identity == "1.00"], col=alpha("blue", 0.5), pch=17, type='b', cex=1, lwd=2, lty=2,
     main="Simulated Engraftment Error", xlab="% Donor Simulated", ylab="Error (Simulation - Observed)")
lines(seq(0,1,0.1), sim_tab$Error[sim_tab$Identity == "1.00"], col=alpha("blue", 0.7), pch=19, type='b', cex=1, lwd=3, lty=2)
lines(seq(0,1,0.1), sim_tab$Error_Raw[sim_tab$Identity == "0.99"], col=alpha("green3", 0.5), pch=17, type='b', cex=1, lwd=2, lty=2)
lines(seq(0,1,0.1), sim_tab$Error[sim_tab$Identity == "0.99"]+0.002, col=alpha("green3", 0.7), pch=19, type='b', cex=1, lwd=3, lty=2)
lines(seq(0,1,0.1), sim_tab$Error_Raw[sim_tab$Identity == "0.98"], col=alpha("purple", 0.5), pch=17, type='b', cex=1, lwd=2, lty=2)
lines(seq(0,1,0.1), sim_tab$Error[sim_tab$Identity == "0.98"], col=alpha("purple", 0.7), pch=19, type='b', cex=1, lwd=3, lty=2)
legend(0.05, -0.1, legend=c("Observed", "Corrected", "1.00", "0.99", "0.98"),
       col=c("darkgrey","grey","blue","green3","purple"), lty = c(2,2,1,1,1), pch=c(17,19,19,19,19), cex=0.8, lwd=2)



# Combine all simulated results
columns <- c("Subject", "Simulation_Donor","Engrafted_Corr","Total_Sum","Expected_Align_Rate","Identity")
sim_tab$Subject <- "MGH06E"
df_all <- rbind(sim_tab[,columns], sim_add[,columns])


# False Negative Rate per Identity - all samples
df_all$Error <- df_all$Engrafted_Corr - df_all$Simulation_Donor
df_all$Identity <- gsub(" identity","", df_all$Identity)
tiff("figures/recall_error_sims.tiff", units="in", width=2.4, height=2.4, res=600)
ggplot(df_all[df_all$Simulation_Donor == 0.5,], aes(x=Identity, y=Error, col=Identity, fill=Identity)) + 
  geom_hline(yintercept=0, col="darkgrey", lwd=1, alpha=0.5) +
  geom_boxplot(lwd=0.4, alpha=0.4, show.legend=F, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.7, size=2, show.legend=F) +
  scale_color_manual(values=rev(c("blue","green3","purple"))) +
  scale_fill_manual(values=rev(c("blue","green3","purple"))) +
  labs(x="Alignment Identity", y="False Negative Rate", title="Recall Error (50% Donor)") +
  theme_bw() +
  theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text.y = element_text(size=6, family="Helvetica"),
        axis.text.x = element_text(size=7, family="Helvetica"))
dev.off()


# Boxplots - Multiple simulation results
df_all$Label <- ""
df_all$Label[df_all$Subject == "MGH07R" & df_all$Simulation_Donor == 0.3] <- "MGH07R"
df_all$Label[df_all$Subject == "MGH11R" & df_all$Simulation_Donor == 0.5] <- "MGH11R"
df_all$Label[df_all$Subject == "MGH01R" & df_all$Simulation_Donor == 0.5] <- "MGH01R"
df_all$Label[df_all$Subject == "MGH17R" & df_all$Simulation_Donor == 0.7] <- "MGH17R"
df_all$Identity <- paste0(df_all$Identity, " identity")
df_all$Engrafted_Corr <- df_all$Engrafted_Corr*100
df_all$Simulation_Donor <- df_all$Simulation_Donor*100
tiff("figures/simulations.tiff", units="in", width=6, height=2.5, res=600)
ggplot(df_all, aes(x=Simulation_Donor, y=Engrafted_Corr, col=Identity,
                   group=Simulation_Donor)) +
  geom_abline(slope=1, intercept=0, col="darkgrey", lwd=1, alpha=0.5) +
  geom_boxplot(alpha=0.2, lwd=0.4, outlier.shape=NA) +
  geom_line(aes(group=interaction(Subject,Identity)), lwd=0.4, alpha=0.4) +
  geom_point(alpha=0.4, size=0.5) +
  scale_color_manual(values=rev(c("blue","green3","purple"))) +
  labs(x="Simulated Donor (%)", y="Engraftment (%)", title="Simulated Engraftment") +
  facet_grid(~Identity) +
  scale_x_continuous(limits=c(-10,110), breaks=seq(0,100,20)) + 
  scale_y_continuous(limits=c(-10,110), breaks=seq(0,100,20)) + 
  #geom_text(aes(label=Label), vjust=1, size=7/.pt) +
  theme_bw() +
  theme(aspect.ratio=1,
        legend.position = "none",
        plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        strip.text = element_text(size=7, family="Helvetica"))
dev.off()


# Quick look at outliers
outliers <- c("MGH07R","MGH11R","MGH01R","MGH17R") # low outliers (2), high outliers (2)
nonoutliers <- unique(df_all$Subject[!(df_all$Subject %in% outliers)])
## any metadata factors explain this?
smillie <- read.delim("data/smillie_metadata/samples.meta.txt", header=T, sep="\t")
smillie[match(outliers, smillie$subject),] # high outliers (01, 17) have lower health scores 1-2 (maximum 10 in study) from questionnaires
smillie[match(nonoutliers, smillie$subject),]
table(smillie$health[smillie$subject %in% df_all$Subject])
triads <- read.delim('data/smillie_metadata/triads.txt', sep="\t", header=T)
table(triads$success[triads$subject %in% outliers], triads$subject[triads$subject %in% outliers])  # low outliers (07, 11) have failures
table(triads$success[triads$subject %in% nonoutliers], triads$subject[triads$subject %in% nonoutliers])
## what about overlap between donor & pre?
out_df <- sim_add[sim_add$Subject %in% outliers,]
out_df[out_df$Simulation_Donor == 0, c("Subject","Total_Donor","Identity")] # low outliers (07,11 have high overlap Donor:Pre)
out_df[out_df$Simulation_Donor == 1, c("Subject","Total_Pre","Identity")]
mean(sim_add$Total_Donor[sim_add$Subject %in% nonoutliers & sim_add$Simulation_Donor == 0])
ovl_df <- sim_add[sim_add$Simulation_Donor == 0,]
ovl_df$Outlier <- ifelse(ovl_df$Subject %in% outliers, "yes", "no")
ggplot(ovl_df, aes(x=Outlier, y=Total_Donor, col=Outlier)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(width=0.25) +
  geom_text(aes(label=Subject)) +
  scale_color_manual(values=c("orange","cyan")) +
  facet_grid(~Identity) +
  theme_bw()
# any failed samples, including low outliers (07,11) and sample 09 all have high overlap with donors

