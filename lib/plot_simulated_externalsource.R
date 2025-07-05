# Given our output tables from the external source, we want
# to evaluate the error introduced by an unknown extra source.
# We expect the error to resemble the overlap between the donor
# and extra source.
library(ggplot2)

# Bring in third source data
third_1 <- read.csv("data/counts_tables/simextra_counts_1.00.txt", header=T, sep=',')
third_99 <- read.csv("data/counts_tables/simextra_counts_0.99.txt", header=T, sep=',')
third_98 <- read.csv("data/counts_tables/simextra_counts_0.98.txt", header=T, sep=',')
third_1$Identity <- "1.00 identity"
third_99$Identity <- "0.99 identity"
third_98$Identity <- "0.98 identity"


# Combine & determine variables for Plotting
third <- rbind(third_1, third_99, third_98)
third$Extra_Source <- as.numeric(gsub("sim.*extra_","",third$Sample))
third$Simulation_Pre <- as.numeric(gsub("sim_","", gsub("_extra_.*","",third$Sample)))
third$Simulation_Donor <- 1 - third$Simulation_Pre
third$True_Sim_Donor <- (1 - third$Extra_Source) * third$Simulation_Donor
third$True_Sim_Pre <- (1 - third$Extra_Source) * third$Simulation_Pre
third$Total_Sum <- third$Total_Donor + third$Total_Pre + third$Unmapped # note: all 8750000
third$Expect_Rate <- rep(c(0.75, 0.81, 0.83), each=121) # % donor reads matching itself per identity
third$Corr_Engraftment <- third$Total_Donor / (third$Total_Sum * third$Expect_Rate)
third$Diff_Engraftment <- third$Corr_Engraftment - third$True_Sim_Donor


# Plot % Donor error vs % Third Source
third$num_id_source <- as.numeric(third$Extra_Source) + rep(c(0.025,0,-0.025), each=121)
third$num_id_source <- jitter(third$num_id_source, factor=1.5)
third$id_samp <- paste0(third$Simulation_Donor, " donor; ", third$Identity)
p <- ggplot(data=third, aes(x=Extra_Source, y=Diff_Engraftment, color=Identity,
                  group=interaction(Extra_Source,Identity))) +
  geom_hline(aes(yintercept=0), col="grey", lwd=0.75) +
  geom_boxplot() +
  geom_path(aes(x=num_id_source, group=id_samp), alpha=0.2) +
  geom_point(aes(x=num_id_source), alpha=0.3, cex=0.5) +
  scale_color_manual(values=c("purple","green3","blue")) +
  labs(title="Sensitivity : External Source Simulations",
       x="External Source %", y="Engraftment Error (Observed - Truth)") +
  ylim(-0.1,0.7) +
  theme_bw() +
  theme(legend.position = c(0.3,0.75),
        title = element_text(face="bold", hjust=0.5))
p


# Add in other external source simulations
## loading additional simulations
con_100 <- read.csv('data/counts_tables/counts_contam_sim_1.00.txt', header=T, sep=',')
con_99 <- read.csv('data/counts_tables/counts_contam_sim_0.99.txt', header=T, sep=',')
con_98 <- read.csv('data/counts_tables/counts_contam_sim_0.98.txt', header=T, sep=',')
con_100$Identity <- "1.00 identity"
con_99$Identity <- "0.99 identity"
con_98$Identity <- "0.98 identity"
con_all <- rbind(con_100, con_99, con_98)
con_all$Simulation_Donor <- con_all$Simulation_Donor / 100
con_all$Extra_Source <- 0.2
con_all$True_Sim_Donor <- (1 - con_all$Extra_Source) * con_all$Simulation_Donor
con_all$Total_Sum <- con_all$Total_Donor + con_all$Total_Pre + con_all$Unmapped
con_all$Expect_Rate <- sim_add$Expected_Align_Rate
con_all$Corr_Engraftment <- con_all$Total_Donor / (con_all$Total_Sum * con_all$Expect_Rate)
con_all$Diff_Engraftment <- con_all$Corr_Engraftment - con_all$True_Sim_Donor


# 0.2 Extra Source Accuracy
columns <- c("Subject","True_Sim_Donor","Simulation_Donor","Extra_Source","Corr_Engraftment","Diff_Engraftment","Identity")
third$Subject <- "MGH06E"
twenty <- rbind(third[third$Extra_Source == 0.2, columns], con_all[,columns])
twenty$Identity <- gsub(" identity","",twenty$Identity)
tiff("figures/accuracy_error_extrasim.tiff", units="in", width=2.4, height=2.4, res=600)
ggplot(twenty[twenty$Simulation_Donor == 0.5,], aes(x=Identity, y=Diff_Engraftment, col=Identity, fill=Identity)) + 
  geom_hline(yintercept=0, col="darkgrey", lwd=1, alpha=0.5) +
  geom_boxplot(lwd=0.6, alpha=0.4, outlier.shape=NA, show.legend=F) + 
  geom_jitter(width=0.25, alpha=0.7, size=2, show.legend=F) +
  scale_color_manual(values=rev(c("blue","green3","purple"))) +
  scale_fill_manual(values=rev(c("blue","green3","purple"))) +
  scale_y_continuous(breaks=seq(-0.1, 0.21, 0.05)) +
  labs(x="Alignment Identity", y="False Positive Rate", title="Accuracy (20% Contaminate, 40% Donor)") +
  #geom_text(aes(label=Subject), size=7/.pt, show.legend=F) + # note that sample 17 has outlying error above, sample 11 outlying error under
  theme_bw() +
  theme(plot.title=element_text(face="bold", size=7, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text.y = element_text(size=6, family="Helvetica"),
        axis.text.x = element_text(size=7, family="Helvetica"))
dev.off()










# Additional (exploratory) Expected over-estimation of engraftment : 
#     (% overlap between donor & extra) * (% extra in sample)
# overlap <- rep(c(mean(third$Corr_Engraftment[third$Extra_Source == 1.0 & third$Identity == "1.00 identity"]),
#                  mean(third$Corr_Engraftment[third$Extra_Source == 1.0 & third$Identity == "0.99 identity"]),
#                  mean(third$Corr_Engraftment[third$Extra_Source == 1.0 & third$Identity == "0.98 identity"])), each=nrow(third_1))
# third$Error_expected <- third$Extra_Source * overlap
# third$Err_Obs_Exp <- third$Diff_Engraftment - third$Error_expected
# p <- ggplot(data=third, aes(x=Extra_Source, y=Diff_Engraftment, color=Identity,
#                             group=interaction(Extra_Source,Identity))) +
#   geom_hline(aes(yintercept=0), col="grey", lwd=0.75) +
#   geom_boxplot() +
#   geom_path(aes(x=num_id_source, group=id_samp), alpha=0.2) +
#   geom_point(aes(x=num_id_source), alpha=0.3, cex=0.5) +
#   scale_color_manual(values=c("purple","green3","blue")) +
#   labs(title="Sensitivity : External Source Simulations",
#        x="External Source %", y="Engraftment Error (Observed - Truth)") +
#   ylim(-0.1,0.7) +
#   theme_bw()
# p + geom_boxplot(aes(x=Extra_Source, y=Err_Obs_Exp, group=Extra_Source), col="grey", alpha=0.3)
# 
# 


