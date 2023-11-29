library(ggplot2)
library(plyr)

# expected command line args
args = commandArgs(trailingOnly=TRUE)
srr <- args[1]
ref <- args[2]
srr_ref <- paste(srr, ref, sep="_")

ref_inserts = args[3]
human_inserts = args[4]


data1 <- read.tabledata <- read.table(human_inserts, header = F)
data2 <- read.tabledata <- read.table(ref_inserts, header = F)
length_cutoff <- 500
d1 <- data1$V1[data1$V1 <= length_cutoff]
d2 <- data2$V1[data2$V1 <= length_cutoff]

ref_legend <- paste(ref, paste("(", length(d2), ")", sep=""), sep=" ")
human_legend <- paste("Human", paste("(", length(d1), ")", sep="")) #ifelse(length(d1)<=100000, paste("(", length(d1), ")", sep=""), "(100000)"), sep=" ")

df1 <- data.frame(PLOT_DATA=human_legend, V2=d1 )
df1 <- df1[sample(nrow(df1), 100000), ]
df2 <- data.frame(PLOT_DATA=ref_legend, V2=d2)

df <- rbind(df1, df2)
colnames(df)[1] <- srr
PLOT_COL <- sym(srr)

# plot with median
mu <- ddply(df, srr, summarise, grp.median=median(V2))
p <- ggplot(df, aes(x=V2, fill=!!PLOT_COL))+geom_density(alpha=.25)+scale_fill_manual(values=c("blue","red"))+
  scale_x_continuous(name = 'Fragment Size, bp',limits=c(50,500),breaks = seq(50,500,50),position ='bottom')+
  scale_y_continuous(name = 'Fraction of Reads',limits=c(0,0.035))+
  theme_bw()+theme(panel.grid=element_blank())+
  geom_vline(data=mu, aes(xintercept=grp.median), color=c("blue","red"), linetype="dashed") + theme(aspect.ratio=1)


median_png <- paste(srr_ref, "dedup_inserts", "median.png", sep="_")
ggsave(
  median_png,
  plot = p
)

# plot with mode
getmode <- function(v) { 
  uniqv <- unique(v) 
  uniqv[which.max(tabulate(match(v, uniqv)))] 
}
mo1 <- getmode(d1)
mo2 <- getmode(d2)
p <- ggplot(df, aes(x=V2, fill=!!PLOT_COL))+geom_density(alpha=.25)+scale_fill_manual(values=c("blue","red"))+
  scale_x_continuous(name = 'Fragment Size, bp',limits=c(50,500),breaks = seq(50,500,50),position ='bottom')+
  scale_y_continuous(name = 'Fraction of Reads',limits=c(0,0.035))+
  theme_bw()+theme(panel.grid=element_blank())+
  geom_vline(aes(xintercept=mo1), color="blue", linetype="dashed")+
  geom_vline(aes(xintercept=mo2), color="red", linetype="dashed") + theme(aspect.ratio=1)

mode_png <- paste(srr_ref, "dedup_inserts", "mode.png", sep="_")
ggsave(
  mode_png,
  plot = p
)
