#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# GLM
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

glm.matrix.neg<-read.csv(args[1], header=TRUE)

glm.matrix.neg.subset1<-subset(glm.matrix.neg, length <= 40)
glm.matrix.neg.subset2<-subset(glm.matrix.neg, length >= 50)

glm.matrix.neg.subset2$length[glm.matrix.neg.subset2$length == 50] <- 41
glm.matrix.neg.subset2$length[glm.matrix.neg.subset2$length == 51] <- 42
glm.matrix.neg.subset2$length[glm.matrix.neg.subset2$length == 52] <- 43

glm.merge<-rbind(glm.matrix.neg.subset1,glm.matrix.neg.subset2)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# RF
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

rf.matrix<-read.csv(args[2], header=TRUE)

rf.matrix.subset1<-subset(rf.matrix, length <= 40)
rf.matrix.subset2<-subset(rf.matrix, length >= 50)

rf.matrix.subset2$length[rf.matrix.subset2$length == 50] <- 41
rf.matrix.subset2$length[rf.matrix.subset2$length == 51] <- 42
rf.matrix.subset2$length[rf.matrix.subset2$length == 52] <- 43

rf.merge<-rbind(rf.matrix.subset1,rf.matrix.subset2)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Plot heatmaps - Draw boxes around the different sets fo predictors
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

box_glm <-ggplot(data=glm.merge, aes(x=position, y=length, fill=standardized_coefficients)) + 
  geom_tile() + 
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=42.5,ymax=43.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=41.5,ymax=42.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=-0.5,ymin=40.5,ymax=41.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-0.5,xmax=0.5,ymin=40.5,ymax=41.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=0.5,xmax=20.5,ymin=40.5,ymax=41.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=19.5,ymax=40.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=18.5,ymax=19.5), color="black", size=0.5, fill=NA) +
  scale_x_continuous(limits = c(-21,21), expand = c(0, 0), breaks=c(-20,-15,-10,-5,0,5,10,15,20)) +
  scale_y_continuous(limits = c(18,44), expand = c(0, 0), breaks=c(19,20,25,30,35,40,41,42,43), labels=c("Nucleotide_sequence","20","25","30","35","40","Proportion_of_reads","ORF_FPKM","Number_of_upstream_TIS")) + 
  scale_fill_gradientn(colours=c("#5e3c99", "#b2abd2", "white", "#fee0b6", "#f1a340", "#b35806"), values=c(-6,-0.0001,0,0.0001,0.1,6), rescaler = function(x,...) x, oob = identity, name="Score") +
  xlab("Position Relative to Start Codon") + 
  ylab("Fragment Length") +
  ggtitle("Scaled coefficients from the GLM") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

box_rf <-ggplot(data=rf.merge, aes(x=position, y=length, fill=value)) + 
  geom_tile() + 
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=42.5,ymax=43.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=41.5,ymax=42.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=-0.5,ymin=40.5,ymax=41.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-0.5,xmax=0.5,ymin=40.5,ymax=41.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=0.5,xmax=20.5,ymin=40.5,ymax=41.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=19.5,ymax=40.5), color="black", size=0.5, fill=NA) +
  geom_rect(aes(xmin=-20.5,xmax=20.5,ymin=18.5,ymax=19.5), color="black", size=0.5, fill=NA) +
  scale_x_continuous(limits = c(-21,21), expand = c(0, 0), breaks=c(-20,-15,-10,-5,0,5,10,15,20)) +
  scale_y_continuous(limits = c(18,44), expand = c(0, 0), breaks=c(19,20,25,30,35,40,41,42,43), labels=c("Nucleotide_sequence","20","25","30","35","40","Proportion_of_reads","ORF_FPKM","Number_of_upstream_TIS")) + 
  scale_fill_gradientn(colours=c("white", "#ece2f0", "#a6bddb", "#1c9099"), values=c(0,0.05,0.1,1), rescaler = function(x,...) x, oob = identity, name="Score") +
  xlab("Position Relative to Start Codon") + 
  ylab("Fragment Length") +
  ggtitle("Scaled variable importance from the random forest") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# output
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

both.heatmaps<-grid.arrange(box_glm, box_rf , ncol=2)
ggsave(file=args[3], both.heatmaps, width = 400, height = 120, unit="mm", dpi = 1200)
