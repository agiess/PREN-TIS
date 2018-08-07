#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)

##
# 1 unshifted heatmaps
##

start_m5<-read.table(args[1])
start_m3<-read.table(args[2])

mono.start5.raw<-read.csv(args[3], header=FALSE, check.names=FALSE, sep=",")
mono.start3.raw<-read.csv(args[4], header=FALSE, check.names=FALSE, sep=",")

start_m5_50 <- subset(start_m5, V1<=30 & V1 >= -30);
start_m5_50 <- start_m5_50 %>% group_by(V2) %>% mutate(LenSum=sum(V3), LenSD=sd(V3), LenMean=mean(V3)) %>% mutate(LenZ=(V3-LenMean)/LenSD)
start_m5_50 <- start_m5_50 %>% group_by(V1) %>% mutate(PosSum=sum(V3), PosSD=sd(V3), PosMean=mean(V3)) %>% mutate(PosZ=(V3-PosMean)/PosSD)

start_m3_50 <- subset(start_m3, V1<=50 & V1 >= -10);
start_m3_50 <- start_m3_50 %>% group_by(V2) %>% mutate(LenSum=sum(V3), LenSD=sd(V3), LenMean=mean(V3)) %>% mutate(LenZ=(V3-LenMean)/LenSD)
start_m3_50 <- start_m3_50 %>% group_by(V1) %>% mutate(PosSum=sum(V3), PosSD=sd(V3), PosMean=mean(V3)) %>% mutate(PosZ=(V3-PosMean)/PosSD)

mLen5<-ggplot(subset(start_m5_50, V2 <= 40 & V2 >= 20 & V1 <= 30 & V1 >= -30) , aes(x=V1, y=V2, fill=LenZ)) + geom_tile()  +
               scale_fill_gradientn(colours=c("grey", "grey", "white", "white", "red1", "red2", "red3", "red3"), values=c(-4,-2,0,2,4,6,8,10), rescaler = function(x,...) x, oob = identity, name="Z-score by fragment length") +
               xlab("Position relative to start codon") + ylab("Fragment length") +
               scale_x_continuous(limits = c(-31,31), expand = c(0, 0), breaks=c(-30, -20, -10, 0, 10, 20, 30)) +
               scale_y_continuous(limits = c(19,41), expand = c(0, 0), breaks=c(20, 25, 30, 35, 40)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               theme(text = element_text(size = 12)) +
               theme(panel.background=element_rect(fill="black"))

mLen3<-ggplot(subset(start_m3_50, V2 <= 40 & V2 >= 20 & V1 <= 50 & V1 >= -10) , aes(x=V1, y=V2, fill=LenZ)) + geom_tile()  +
               scale_fill_gradientn(colours=c("grey", "grey", "white", "white", "red1", "red2", "red3", "red3"), values=c(-4,-2,0,2,4,6,8,10), rescaler = function(x,...) x, oob = identity, name="Z-score by fragment length") +
               xlab("Position relative to start codon") + ylab("Fragment length") +
               scale_x_continuous(limits = c(-11,51), expand = c(0, 0), breaks=c(-10, 0, 10, 20, 30, 40, 50)) +
               scale_y_continuous(limits = c(19,41), expand = c(0, 0), breaks=c(20, 25, 30, 35, 40)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               theme(text = element_text(size = 12)) +
               theme(panel.background=element_rect(fill="black"))

mPos5<-ggplot(subset(start_m5_50, V2 <= 40 & V2 >= 20 & V1 <= 30 & V1 >= -30) , aes(x=V1, y=V2, fill=PosZ)) + geom_tile()  +
               scale_fill_gradientn(na.value = "black", colours=c("grey", "grey", "white", "white", "red1", "red2", "red3", "red3"), values=c(-3,-2,-1,0,1,2,3,4,5,6,7), rescaler = function(x,...) x, oob = identity, name="Z-score by position") +
               xlab("Position relative to start codon") + ylab("Fragment length") +
               scale_x_continuous(limits = c(-31,31), expand = c(0, 0), breaks=c(-30, -20, -10, 0, 10, 20, 30)) +
               scale_y_continuous(limits = c(19,41), expand = c(0, 0), breaks=c(20, 25, 30, 35, 40)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               theme(text = element_text(size = 12)) +
               theme(panel.background=element_rect(fill="black"))

mPos3<-ggplot(subset(start_m3_50, V2 <= 40 & V2 >= 20 & V1 <= 50 & V1 >= -10) , aes(x=V1, y=V2, fill=PosZ)) + geom_tile()  +
               scale_fill_gradientn(na.value = "black", colours=c("grey", "grey", "white", "white", "red1", "red2", "red3", "red3"), values=c(-3,-2,-1,0,1,2,3,4,5,6,7), rescaler = function(x,...) x, oob = identity, name="Z-score by position") +
               xlab("Position relative to start codon") + ylab("Fragment length") +
               scale_x_continuous(limits = c(-11,51), expand = c(0, 0), breaks=c(-10, 0, 10, 20, 30, 40, 50)) +
               scale_y_continuous(limits = c(19,41), expand = c(0, 0), breaks=c(20, 25, 30, 35, 40)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               theme(text = element_text(size = 12)) +
               theme(panel.background=element_rect(fill="black"))

##
# 2 unshifted barplots
##

#scale the data by the sum of the window
subset_mono5.raw<-filter(mono.start5.raw, V1 >=-30, V1 <= 30)
subset_mono5.raw <- mutate(subset_mono5.raw, proportion = V2/sum(V2))

subset_mono3.raw<-filter(mono.start3.raw, V1 >=-11, V1 <= 51)
subset_mono3.raw <- mutate(subset_mono3.raw, proportion = V2/sum(V2))

max_y <- max(subset_mono3.raw$proportion,subset_mono5.raw$proportion)

mb5 <- ggplot(data=subset_mono5.raw, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
             geom_bar(stat="identity",width=1, colour="black", size=0.6) +
             theme_bw() +
             xlab("Position relative to start codon") +
             ylab("Proportion of reads") +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
             scale_x_continuous(limits = c(-31, 31), expand = c(0, 0), breaks=c(-30, -20, -10, 0, 10, 20, 30)) +
             scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
             scale_fill_manual(values=c("#cc9968", "#cd7130", "#665a50"), name = "Frame") +
             theme(text = element_text(size = 12))

mb3 <- ggplot(data=subset_mono3.raw, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
             geom_bar(stat="identity",width=1, colour="black", size=0.6) +
             theme_bw() +
             xlab("Position relative to start codon") +
             ylab("Proportion of reads") +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
             scale_x_continuous(limits = c(-11, 51), expand = c(0, 0), breaks=c(-10, 0, 10, 20, 30, 40, 50)) +
             scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
             scale_fill_manual(values=c("#cc9868", "#cd7130", "#665a40"), name = "Frame") +
             theme(text = element_text(size = 12))

###
# Pretify plots
###

# Function to save legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Save the legend
#+++++++++++++++++++++++
legend_heatmapP5<- get_legend(mPos5)    
legend_heatmapP3<- get_legend(mPos3)
legend_heatmapL5<- get_legend(mLen5)
legend_heatmapL3<- get_legend(mLen3)
legend_barchart5 <- get_legend(mb5)
legend_barchart3 <- get_legend(mb3)

mPos5 <- mPos5 + theme(legend.position="none")
mPos3 <- mPos3 + theme(legend.position="none")
mLen5 <- mLen5 + theme(legend.position="none")
mLen3 <- mLen3 + theme(legend.position="none")
mb5 <- mb5 + theme(legend.position="none")
mb3 <- mb3 + theme(legend.position="none")

##
# Output
##

lay <- rbind(c(1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,3,3,4,4,4),
             c(5,5,5,5,5,5,5,5,6,6,6,7,7,7,7,7,7,7,7,8,8,8),
             c(9,9,9,9,9,9,9,9,10,10,10,11,11,11,11,11,11,11,11,12,12,12))

grid.arrange(mb5, legend_barchart5, mb3, legend_barchart3, mLen5, legend_heatmapL5, mLen3, legend_heatmapL3, mPos5, legend_heatmapP5, mPos3, legend_heatmapP3, layout_matrix = lay)

out1<-arrangeGrob(mb5, legend_barchart5, mb3, legend_barchart3, mLen5, legend_heatmapL5, mLen3, legend_heatmapL3, mPos5, legend_heatmapP5, mPos3, legend_heatmapP3, layout_matrix = lay)
ggsave(file=,args[5], out1, width=450, height=200, unit="mm", dpi=600) 
