# MAG abundance plots
library(ggplot2)
library(gridExtra)
library(readr)
library(stringr)
library(tidyverse)

##### Loading data #####
# Read in MAG table (only reading in 1.00 identity for now) - using readr and tibble for faster reading
mag_tab <- readr::read_delim('data/mags_tables/binned_mags_table_1.00.txt', col_names=T, delim='\t')
filenames <- c(mag_tab[,1])$FileName

# Read in MAG Bin table
bin_tab <- read.delim('data/metabinner/donor_70a_metabinner_result.tsv', header=F, sep='\t', col.names=c("MAGID","Bin"))
bin_tab$MAGID <- paste0("donor_70a_", bin_tab$MAGID) # 15367 MAGs
tmp <- read.delim('data/metabinner/donor_70b_metabinner_result.tsv', header=F, sep='\t', col.names=c("MAGID","Bin"))
tmp$MAGID <- paste0("donor_70b_", tmp$MAGID) # 20305 MAGs
bin_tab <- rbind(bin_tab, tmp)
tmp <- read.delim('data/metabinner/donor_70c_metabinner_result.tsv', header=F, sep='\t', col.names=c("MAGID","Bin"))
tmp$MAGID <- paste0("donor_70c_", tmp$MAGID) # 22249 MAGs
bin_tab <- rbind(bin_tab, tmp)

# Read in total matching reads from counts table
ct_tab <- read.csv('data/counts_tables/closedlabel_interv_1.00.txt', header=T, sep=',')


##### Cleaning Tables #####
# Reduce MAGs to those in bins
mag_bin <- select(mag_tab, which(colnames(mag_tab) %in% bin_tab$MAGID))
bin_tab <- bin_tab[match(colnames(mag_bin), bin_tab$MAGID),]

# Clean up tables prior to plotting
## get IDs from filenames in MAGs table
samples <- gsub("_donor.*.txt","", filenames)
mag_bin <- mag_bin[samples %in% ct_tab$Sample,] # drop 3 samples
filenames <- filenames[samples %in% ct_tab$Sample]
samples <- samples[samples %in% ct_tab$Sample]
# get donor IDs from MAGIDs
donors <- colnames(mag_bin)
donors <- gsub("_NODE.*","",donors)
table(donors)
donor_simp <- gsub("[a|b|c|d|e|f|g]$", "", donors)
table(donor_simp)
# get table for % engrafted MAGs
ct_tab <- ct_tab[match(samples, ct_tab$Sample),]  # make sure same order as mag_bin
ct_tab$Donor_Group <- gsub(".*_donor_","",gsub("_alignment.*","",filenames))
ct_tab$Donor_ID <- paste0("donor_", gsub("[a|b|c|d|e|f|g]","",ct_tab$Donor_Group))
ct_tab$Week <- as.numeric(gsub(".*W","",ct_tab$Sample))
prc_tab <- sweep(mag_bin, 1, ct_tab$Total_Donor, FUN="/")  # get % of engrafted strains
# simplify to combined total %s per bin
tmp <- t(prc_tab)
tmp <- aggregate(tmp, by=list(bin_tab$Bin), FUN=sum)
rownames(tmp) <- tmp$Group.1
tmp$Group.1 <- NULL
colnames(tmp) <- samples
prc_bin <- t(tmp)

##### Analysis & Plotting #####
# Translation for sample names
translate <- data.frame(Subject=c("UCFMT004","UCFMT009","UCFMT011","UCFMT012","UCFMT016",
                                  "UCFMT017","UCFMT023","UCFMT028","UCFMT032","UCFMT033"),
                        Donor=c("Donor B","Donor A1 ", "Donor A1","Donor A2","Donor A3",
                                "Donor C3","Donor C1","Donor C4","Donor C2","Donor C5") )
# Colors for plotting
ten_colors <- c("#f95700","#f89212","#08653c","#4abb88","#59d43a",
                "#9d02d7","#ec64ff","#ffc5fc","#0000ff","#34c8fd")
color_sets <- list(c(10,3,5,6,2,8,9,7,4,1),  # cool to warm (ish)
                   c(9,4,7,5,8,2,1,10,6,3),  # mixing warm/cool
                   c(2,6,7,4,8,10,3,1,5,9))  # warm to cool (ish)
# Split by donor (specific to donor group)
d_list <- c("70a","70b","70c")
title_prefixes <- c("MAG Bin Abundances: ","MAG Bins: ","MAG Bins: ")
axis_titles <- c("% of Engrafted Reads","","")
bin_ids <- c("A1", "A2","A3")
plot_list <- lapply(d_list,
                    function (curr_donor) {
                      i <- which(d_list == curr_donor)
                      # Restrict to current donor
                      donor_tab <- as.data.frame(prc_bin[ct_tab$Donor_Group == curr_donor,])
                      # Get top 10 bin totals per donor
                      donor_tab <- donor_tab[,order(colSums(donor_tab), decreasing=T)]
                      donor_tab <- donor_tab[,c(1:10)]
                      colnames(donor_tab) <- paste0("Bin ", colnames(donor_tab))
                      donor_tab$Other <- 1 - rowSums(donor_tab)
                      donor_tab$Sample <- rownames(donor_tab)
                      # Melt prior to plotting
                      melted <- reshape2::melt(donor_tab, id.vars="Sample", variable.name="Bin_ID", value.name="Abundance", )
                      melted$Week <- as.numeric(gsub(".*_W","",melted$Sample))
                      melted$Donor <- translate$Donor[match(gsub("_W.*","",melted$Sample), translate$Subject)]
                      melted$Abundance <- melted$Abundance * 100
                      # Plot results
                      p <- ggplot(melted, aes(x=Week, y=Abundance, fill=Bin_ID)) +
                        geom_bar(position="stack", stat="identity") +
                        scale_fill_manual(values=c(ten_colors[color_sets[[i]]], "lightgrey")) +
                        labs(x="Week", y=axis_titles[i],
                             title=paste0(title_prefixes[i],trimws(melted$Donor[1])), 
                             fill="") +
                        scale_x_continuous(limits=c(2,14), breaks=c(4,8,12)) +
                        facet_wrap(~Donor) +
                        guides(fill=guide_legend(nrow=4)) +
                        theme_bw() +
                        theme(plot.title=element_text(face="bold", size=7, family="Helvetica"),
                              axis.title = element_text(face="bold", size=6, family="Helvetica"),
                              axis.text = element_text(size=5, family="Helvetica"),
                              strip.text = element_text(size=6, family="Helvetica"),
                              legend.position = "bottom",
                              legend.box.margin = margin(-5,-10,0,-10),
                              legend.margin = margin(0,0,0,0),
                              legend.key.size = unit(0.3, "lines"),
                              legend.text = element_text(size=5, family="Helvetica"),
                              legend.title = element_text(face="bold", size=6, family="Helvetica"))
                      # Return the plot
                      p
                    })
# Plot all Samples & Donors
tiff("figures/mag_bins_donor70.tiff", units="in", width=6.5, height=3, res=600)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
             widths=c(1.8,1,1), nrow=1)
dev.off()



