# Plotting Output from Closed Label
library(readxl)
library(ggplot2)

##### Get expected rate from self alignments #####
# Self alignment results
exp_100 <- read.csv('data/self_estimates/identity_1.00_estimates_self.txt', sep='\t', header=T)
exp_100$Rate <- exp_100$Count_Align / exp_100$Total_Reads
exp_99 <- read.csv('data/self_estimates/identity_0.99_estimates_self.txt', sep='\t', header=T)
exp_99$Rate <- exp_99$Count_Align / exp_99$Total_Reads
exp_98 <- read.csv('data/self_estimates/identity_0.98_estimates_self.txt', sep='\t', header=T)
exp_98$Rate <- exp_98$Count_Align / exp_98$Total_Reads
# Supplemental : plot the confidence per identity group
df <- data.frame(value=c(exp_100$Rate, exp_99$Rate, exp_98$Rate),
                 id=c(rep("1.00",nrow(exp_100)), rep("0.99", nrow(exp_99)), rep("0.98", nrow(exp_98))))
ggplot(df, aes(x=id, y=value, fill=id, color=id)) +
  geom_boxplot(alpha=0.4, lwd=0.8) +
  geom_jitter(width=0.25, alpha=0.6, size=2) +
  scale_color_manual(values=c("purple","green3","blue")) +
  scale_fill_manual(values=c("purple","green3","blue")) +
  labs(title="Self-Mapping Reads (Donor & Pre Samples)", y="Proportion Reads Aligned",
       x="Alignment Identity", color="Identity", fill="Identity") +
  ylim(0.75, 1.0) +
  theme_bw()


##### Loading Counts Results #####
# Get intervention results per identity
tab_100 <- read.csv('data/counts_tables/closedlabel_interv_1.00.txt', header=T, sep=',')
tab_99 <- read.csv('data/counts_tables/closedlabel_interv_0.99.txt', header=T, sep=',')
tab_98 <- read.csv('data/counts_tables/closedlabel_interv_0.98.txt', header=T, sep=',')
# Get placebo results per identity
plc_100 <- read.csv('data/counts_tables/placebo_1.00.txt', header=T, sep=",")
plc_99 <- read.csv('data/counts_tables/placebo_0.99.txt', header=T, sep=",")
plc_98 <- read.csv('data/counts_tables/placebo_0.98.txt', header=T, sep=",")
# Get wk0 results per identity
wk0_100 <- read.csv('data/counts_tables/counts_wk0_1.00.txt', header=T, sep=",")
wk0_99 <- read.csv('data/counts_tables/counts_wk0_0.99.txt', header=T, sep=",")
wk0_98 <- read.csv('data/counts_tables/counts_wk0_0.98.txt', header=T, sep=",")


##### Combine Table & Add Relevant Variables #####
combine <- rbind(tab_100, tab_99, tab_98, plc_100, plc_99, plc_98, wk0_100, wk0_99, wk0_98)
combine$Identity <- c(rep("1.00 identity", nrow(tab_100)), rep("0.99 identity", nrow(tab_99)), rep("0.98 identity", nrow(tab_98)),
                      rep("1.00 identity", nrow(plc_100)), rep("0.99 identity", nrow(plc_99)), rep("0.98 identity", nrow(plc_98)),
                      rep("1.00 identity", nrow(wk0_100)), rep("0.99 identity", nrow(wk0_99)), rep("0.98 identity", nrow(wk0_98)))
combine$Type <- c(rep("intervention", nrow(tab_100) + nrow(tab_99) + nrow(tab_98)),
                  rep("placebo", nrow(plc_100) + nrow(plc_99) + nrow(plc_98)),
                  rep(c(rep("intervention",10), rep("placebo",12)), 3))
# --- reduce to closed label
combine <- combine[!(grepl("OL", combine$Sample)),]
# --- add general ID
combine$Sample <- gsub("W","_W", combine$Sample)
combine$Sample <- gsub("BL","_BL", combine$Sample)
combine$Sample <- gsub("__","_", combine$Sample)
combine$ID = gsub("_.*", "", combine$Sample)
# --- calculate week from Sample ID
combine$Week <- as.numeric(gsub("UCFMT[0-9][0-9][0-9]_W","", combine$Sample))
combine$Week[is.na(combine$Week)] <- 0
combine <- combine[order(combine$Week),]
combine <- combine[order(combine$ID),]
# --- calculate sum per simulation
combine$Total_Sum <- combine$Total_Donor + combine$Total_Pre + combine$Unmapped
combine$Raw_Perc_Donor <- combine$Total_Donor / combine$Total_Sum
combine$Raw_Perc_Pre <- combine$Total_Pre / combine$Total_Sum
combine$Raw_Perc_Unk <- combine$Unmapped / combine$Total_Sum
combine$Expected_Rate <- c(mean(exp_98$Rate), mean(exp_99$Rate), mean(exp_100$Rate))[factor(combine$Identity)]
combine$Corr_Perc_Donor <- combine$Total_Donor / (combine$Total_Sum * combine$Expected_Rate)
combine$Corr_Perc_Pre <- combine$Total_Pre / (combine$Total_Sum * combine$Expected_Rate)
cl_d <- read.csv("data/comparison_tables/post_table.txt", sep="\t", header=T)
pl_d <- read.csv("data/comparison_tables/placebo_post_table.txt", sep="\t", header=T)
combine$Donor_ID <- c(cl_d$Donor_Set, pl_d$Donor_Set)[match(combine$ID, c(cl_d$Post_ID, pl_d$Post_ID))]


##### Visualize Closed Label Samples #####
# Engraftment per donor per sample
rownames(combine) <- 1:nrow(combine)
combine$Identity_Type <- factor(paste(combine$Identity, combine$Type, combine$ID))
combine <- combine[!(combine$Donor_ID %in% c("102c","102e")),]
ggplot(combine, aes(x=Week, y=Corr_Perc_Donor, color=Type, shape=Identity, linetype=Identity)) +
  geom_line(aes(group=Identity_Type), alpha=0.7, lwd=1) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=c("blue","pink")) +
  scale_shape_manual(values=c(16,17,18)) +
  facet_wrap(~Donor_ID, nrow=1) +
  labs(title="Engraftment per sample (Corrected)", x="Week", y="Engraftment (Mapped Reads)",
       col="Randomization", shape="Identity", linetype="Identity") +
  xlim(0,15) +
  theme_bw()


##### Compare Source Tracker & Our Results #####
# Get the Source Tracker results
st <- read_xlsx('data/sample_labels_moutsoglou/all pre post maps.xlsx', sheet=1,
                range="A29:D139", trim_ws=T, col_names=c("Sample","Donor","Pre","Unknown"))
st <- st[!is.na(st$Sample),]      # omit 10 blank lines
st <- st[!grepl("#", st$Sample),] # omit comment lines
st$Type <- "intervention"
st_plc <- read.csv('data/sample_labels_moutsoglou/st16s_placebo_results.csv', header=T,
                   col.names = c("Sample","Donor","Pre","Unknown"))
st_plc$Type <- "placebo"
st <- rbind(st, st_plc)
st$ID <- gsub("_.*", "", st$Sample)
st <- st[st$ID %in% combine$ID,]  # keep 154 samples
st <- st[st$Sample != "UCFMT017_W12a",] # drop extra sample
st$Sample <- gsub("b","", st$Sample)
st$Type <- paste0(st$Type, " (16S)")
st$Week <- as.numeric(gsub(".*_W","", st$Sample))
st$Identity <- "16S"
st$Donor_ID <- combine$Donor_ID[match(st$ID,combine$ID)]
st <- as.data.frame(st[!is.na(st$Week),])
st$Corr_Perc_Donor <- round(as.numeric(st$Donor),3)

# Combine datasets for comparison (keep only 0.98 identity results for now)
combine$Type <- paste0(combine$Type, " (WGS)")
comp <- rbind(st[,c("Sample","ID","Donor_ID","Week","Corr_Perc_Donor","Type","Identity")],
              combine[,c("Sample","ID","Donor_ID","Week","Corr_Perc_Donor","Type","Identity")])
comp <- comp[!(comp$ID %in% c("UCFMT001")),]
comp <- comp[comp$Identity %in% c("0.98 identity","16S"),]

# Translate Donor IDs to better names
comp$Donor_ID_Uniq <- comp$Donor_ID
comp$Donor_ID_Uniq[comp$ID %in% c("UCFMT009","UCFMT014")] <- "70z"
translate <- data.frame(Original=c("71","70a","70b","70c","70z","102a","102b","102c","102d","102e","102f","102g"),
                        New=c("Donor B","Donor A1","Donor A2","Donor A3","Donor A1 ",
                              "Donor C3","Donor C1","Donor C-","Donor C4","Donor C-", "Donor C2", "Donor C5"))
comp$Donor_ID_Uniq <- translate$New[match(comp$Donor_ID_Uniq, translate$Original)]

# Plot Together - Engraftment Lines
comp$Identity_Type <- paste0(comp$Identity, comp$Type, comp$ID)
ggplot(comp, aes(x=Week, y=Corr_Perc_Donor, color=Type)) +
  geom_line(aes(group=Identity_Type), alpha=0.7, lwd=1) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=c("#708090","blue","#e5cece","pink2")) +
  scale_x_continuous(breaks=c(0,4,8,12), limits=c(0,15)) +
  facet_wrap(~Donor_ID_Uniq) +
  labs(title="Engraftment per Participant : Comparison to SourceTracker",
       x="Week", y="% Engraftment (Corrected)") +
  theme_bw() +
  theme(title=element_text(face="bold", hjust = 0.5))

# Plot Together But Separate Placebo & intervention - Engraftment Lines
comp$Rand_Type <- gsub(" [(].*$","", comp$Type)
comp$Identity_Type <- paste0(comp$Identity, comp$Type, comp$ID)
comp$Corr_Perc_Donor <- comp$Corr_Perc_Donor * 100
tiff("figures/supplemental_st_engraft_lines.tiff", units="in", width=6.5, height=4, res=600)
ggplot(comp, aes(x=Week, y=Corr_Perc_Donor, color=Type)) +
  geom_line(aes(group=Identity_Type), alpha=0.7, lwd=1) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=c("#708090","blue","#e5cece","pink2")) +
  scale_x_continuous(breaks=c(0,4,8,12), limits=c(0,15)) +
  facet_grid(Rand_Type ~ Donor_ID_Uniq) +
  labs(title="Engraftment per Participant : Comparison to SourceTracker",
       x="Week", y="% Engraftment", col="Group (Data Type)") +
  theme_bw() +
  theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        strip.text.x = element_text(size=7, family="Helvetica"),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,0,-10),
        legend.text = element_text(size=7, family="Helvetica"),
        legend.title = element_text(size=7, family="Helvetica"))
dev.off()

# Plot Together - Engraftment Lines (only weeks 0,4,8,12)
ggplot(comp[comp$Week %in% c(4,8,12),], aes(x=Week, y=Corr_Perc_Donor, color=Type)) +
  geom_line(aes(group=Identity_Type), alpha=0.7, lwd=1) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=c("#708090","blue","#e5cece","pink2")) +
  scale_x_continuous(breaks=c(0,4,8,12), limits=c(0,15)) +
  facet_wrap(~Donor_ID_Uniq, nrow=1) +
  labs(title="Engraftment per Participant: Comparison to SourceTracker", x="Week", y="% Engraftment (Corrected)") +
  theme_bw() +
  theme(title=element_text(face="bold", hjust = 0.5))

# Stat compare same samples : differences in donor values
wgs_ovlp <- combine[(combine$Sample %in% st$Sample) & (combine$Identity == "0.98 identity"),]
st_ovlp <- st[match(wgs_ovlp$Sample, st$Sample),]
reduce <- data.frame(result_WGS = wgs_ovlp$Corr_Perc_Donor,
                     result_16S = st_ovlp$Corr_Perc_Donor,
                     Type = wgs_ovlp$Type)
reduce$Type <- gsub(" [(]WGS[)]", "", as.character(reduce$Type))
reduce$result_WGS <- reduce$result_WGS * 100
reduce$result_16S <- reduce$result_16S * 100
p <- ggplot(reduce, aes(x=result_WGS, y=result_16S, col=Type, fill=Type)) +
  geom_abline(slope=1, intercept=0, col="grey", lwd=1) +
  stat_smooth(method=lm, fullrange=T, se=F) +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values=c("blue","pink3")) +
  scale_fill_manual(values=c("lightblue","pink")) +
  xlim(0,100) + ylim(0,100) +
  labs(title="Agreement with SourceTracker", x="WGS result (%)", y="16S rRNA gene result (%)", col="", fill="") +
  theme_bw()
p + theme(title=element_text(face="bold", hjust = 0.5))
cor.test(reduce$result_16S, reduce$result_WGS, method="pearson") # p < 0.001***
cor.test(reduce$result_16S[reduce$Type=="intervention"],
         reduce$result_WGS[reduce$Type=="intervention"], method="pearson") # p < 0.001***
cor.test(reduce$result_16S[reduce$Type=="placebo"],
         reduce$result_WGS[reduce$Type=="placebo"], method="pearson") # p = 0.445
tiff('figures/st_agreement.tiff', units='in', width=3.5, height=2.5, res=600)
p +
  annotate("text", x=0.87, y=0.64, label="p < 0.001***", col="blue", size=7/.pt, family="Helvetica")+
  annotate("text", x=0.83, y=0.2, label="p = 0.445", col="pink2", size=7/.pt, family="Helvetica") +
  annotate("text", x=0.25, y=0.85, label="together\np < 0.001***", col="darkgrey", size=7/.pt, family="Helvetica") +
  theme(plot.title=element_text(face="bold", size=8, family="Helvetica"),
        axis.title = element_text(face="bold", size=7, family="Helvetica"),
        axis.text = element_text(size=6, family="Helvetica"),
        legend.text = element_text(size=6, family="Helvetica"),
        legend.title = element_text(size=7, family="Helvetica"))
dev.off()


