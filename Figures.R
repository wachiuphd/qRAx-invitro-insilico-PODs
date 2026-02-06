library(tidyverse)
library(GGally)
library(scales)

res.wide.df <- read.csv("PODs.bychem.source-withUFD.csv")
res.wide.df$Chemical <- factor(res.wide.df$Chemical,
                            levels=res.wide.df$Chemical[order(res.wide.df$Chemical)])
res.wide.df$CASRN <- factor(res.wide.df$CASRN,
                       levels=res.wide.df$CASRN[order(res.wide.df$Chemical)])

res.df <- read.csv("PODs.bychem-withUFD.csv")
res.df$Chemical <- factor(res.df$Chemical,
                               levels=res.wide.df$Chemical[order(res.wide.df$Chemical)])
res.df$CASRN <- factor(res.df$CASRN,
                       levels=res.wide.df$CASRN[order(res.wide.df$Chemical)])

res.df$name[res.df$name=="PPRTV"] <- "qRAx (PPRTV)"
res.df$name[res.df$name=="Aurisano"] <- "ToxValDB (Aurisano)"
res.df$name[res.df$name=="ToxCast"] <- "IVIVE (ToxCast)"
res.df$name[res.df$name=="CTV"] <- "QSAR (Wignall)"
res.df$name[res.df$name=="Kvasnicka"] <- "QSAR (Kvasnicka)"
res.df$name[res.df$name=="von Borries"] <- "QSAR (von Borries)"
res.df$name <- factor(res.df$name,
                           levels=c("QSAR (von Borries)",
                                    "QSAR (Kvasnicka)",
                                    "QSAR (Wignall)",
                                    "IVIVE (ToxCast)",
                                    "ToxValDB (Aurisano)",
                                    "qRAx (PPRTV)"))

res.success <- data.frame(name=names(res.wide.df[,-(1:4)]),
  n=apply(res.wide.df[,-(1:4)],2,function(x) { sum(!is.na(x))}))[1:7,]
res.success$name <- gsub("pod.med_","",res.success$name)
res.success$name[res.success$name=="PPRTV"] <- "qRAx (PPRTV)"
res.success$name[res.success$name=="Aurisano"] <- "ToxValDB (Aurisano)"
res.success$name[res.success$name=="ToxCast"] <- "IVIVE (ToxCast)"
res.success$name[res.success$name=="CTV"] <- "QSAR (Wignall)"
res.success$name[res.success$name=="Kvasnicka"] <- "QSAR (Kvasnicka)"
res.success$name[res.success$name=="von.Borries"] <- "QSAR (von Borries)"
res.success$name <- factor(res.success$name,
                           levels=c("QSAR (von Borries)",
                                    "QSAR (Kvasnicka)",
                                    "QSAR (Wignall)",
                                    "IVIVE (ToxCast)",
                                    "ToxValDB (Aurisano)",
                                    "qRAx (PPRTV)"))

psuccess <-
  ggplot(subset(res.success,!is.na(name)))+
    geom_col(aes(x=name,y=n,fill=name))+
    geom_label(aes(x=name,y=n,label=paste0(round(100*n/41),"%")),hjust=1)+
    coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
    theme_bw()+theme(legend.position = "none")+
    ggtitle("Successful PODs out of 41 attempted substances")+
    xlab("")+ylab("n success (%)")+geom_hline(yintercept = 41,linetype="dashed")
print(psuccess)
ggsave("Fig2-success.pdf",psuccess,height=4,width=6,scale=0.75)

chem.names <- as.character(res.df$Chemical)
names(chem.names)<-res.df$CASRN
ppods <-
  ggplot(subset(res.df,!is.na(name)))+geom_point(aes(x=pod.med,y=name,color=name))+
    geom_errorbarh(aes(xmin=pod.med.lb,xmax=pod.med,y=name,color=name))+
    theme_bw()+theme(legend.position = "none")+
    scale_color_viridis_d(direction = -1, end=0.9)+
    facet_wrap(~CASRN,labeller=labeller(CASRN=chem.names),ncol=3)+
    scale_x_log10("Human Equivalent Dose POD (mg/kg-d)")
print(ppods)
ggsave("Fig3-PODs.bychem.pdf",ppods,height=10,width=6,scale=1.5)

## Most sensitive

min.pods <- data.frame(
  name=rev(c("QSAR (von Borries)",
             "QSAR (Kvasnicka)",
             "QSAR (Wignall)",
             "IVIVE (ToxCast)",
             "ToxValDB (Aurisano)",
             "qRAx (PPRTV)")),
  value=table(apply(res.wide.df[,c(5,8,11,7,9,10)],1,which.min))[paste0(1:6)]
) 
min.pods$name <- factor(min.pods$name,
                      levels=c("QSAR (von Borries)",
                               "QSAR (Kvasnicka)",
                               "QSAR (Wignall)",
                               "IVIVE (ToxCast)",
                               "ToxValDB (Aurisano)",
                               "qRAx (PPRTV)"))
pmin <-
  ggplot(min.pods)+
  geom_col(aes(x=name,y=value.Freq,fill=name))+
  geom_label(aes(x=name,y=value.Freq,label=paste0(round(100*value.Freq/41),"%")))+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Minimum PODs")+
  xlab("")+ylab("Frequency out of 41 (%)")
print(pmin)
ggsave("Fig3B-mostsens.med.pdf",pmin,height=4,width=6,scale=0.75)

min.pods.lb <- data.frame(
  name=rev(c("QSAR (von Borries)",
             "QSAR (Kvasnicka)",
             "QSAR (Wignall)",
             "IVIVE (ToxCast)",
             "ToxValDB (Aurisano)",
             "qRAx (PPRTV)")),
  value=table(apply(res.wide.df[,7+c(5,8,11,7,9,10)],1,which.min))[paste0(1:6)]
) 
min.pods.lb$value.Freq[is.na(min.pods.lb$value.Freq)]<-0
min.pods.lb$name <- factor(min.pods.lb$name,
                        levels=c("QSAR (von Borries)",
                                 "QSAR (Kvasnicka)",
                                 "QSAR (Wignall)",
                                 "IVIVE (ToxCast)",
                                 "ToxValDB (Aurisano)",
                                 "qRAx (PPRTV)"))
pmin.lb <-
  ggplot(min.pods.lb)+
  geom_col(aes(x=name,y=value.Freq,fill=name))+
  geom_label(aes(x=name,y=value.Freq,label=paste0(round(100*value.Freq/41),"%")))+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Minimum PODs (lower bound)")+
  xlab("")+ylab("Frequency out of 41 (%)")
print(pmin.lb)
ggsave("Fig3B-mostsens.lb.pdf",pmin.lb,height=4,width=6,scale=0.75)


min.pods.lb.med <- data.frame(
  name=rev(c("QSAR (von Borries)",
             "QSAR (Kvasnicka)",
             "QSAR (Wignall)",
             "IVIVE (ToxCast)",
             "ToxValDB (Aurisano)",
             "qRAx (PPRTV)")),
  value=table(apply(res.wide.df[,c(12,8,11,7,9,10)],1,which.min))[paste0(1:6)]
) 
min.pods.lb.med$value.Freq[is.na(min.pods.lb$value.Freq)]<-0
min.pods.lb.med$name <- factor(min.pods.lb.med$name,
                           levels=c("QSAR (von Borries)",
                                    "QSAR (Kvasnicka)",
                                    "QSAR (Wignall)",
                                    "IVIVE (ToxCast)",
                                    "ToxValDB (Aurisano)",
                                    "qRAx (PPRTV)"))

pmin.lb.med <-
  ggplot(min.pods.lb.med)+
  geom_col(aes(x=name,y=value.Freq,fill=name))+
  geom_label(aes(x=name,y=value.Freq,label=paste0(round(100*value.Freq/41),"%")))+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Minimum PODs (PPRTV lb)")+
  xlab("")+ylab("Frequency out of 41 (%)")
print(pmin.lb.med)
ggsave("Fig3B-mostsens.lb.med.pdf",pmin.lb.med,height=4,width=6,scale=0.75)


## Confidence intervals

res.df$name.rev <- factor(res.df$name,
                          levels=rev(levels(res.df$name)))
pci<-
  ggplot(subset(res.df,!is.na(name))) + 
  geom_histogram(aes(x=log10(pod.med/pod.med.lb),fill=name),
                 breaks=seq(-0.25,3.25,0.5)) +
  scale_fill_viridis_d(direction = -1, end=0.9)+
  facet_wrap(~name.rev,ncol=1)+
  scale_x_continuous(breaks=0:3,labels=c(1,10,100,1000))+
  xlab("Fold uncertainty")+
  annotation_logticks(side="b",outside=TRUE,
                      short=unit(0.02,"cm"),
                      mid=unit(0.05,"cm"),
                      long=unit(0.1,"cm"))+
  coord_cartesian(clip = "off")+
  theme_bw()+theme(legend.position = "none")

pcibox<-
  ggplot(subset(res.df,!is.na(name))) + 
  geom_boxplot(aes(y=name,x=log10(pod.med/pod.med.lb),fill=name)) +
  scale_fill_viridis_d(direction = -1, end=0.9,alpha = 0.75)+
  ggtitle("Precision of PODs")+
  scale_x_continuous(breaks=0:3,labels=c(1,10,100,1000))+
  xlab("Fold uncertainty")+
  annotation_logticks(side="b",outside=TRUE,
                      short=unit(0.02,"cm"),
                      mid=unit(0.05,"cm"),
                      long=unit(0.1,"cm"))+
  coord_cartesian(clip = "off")+
  ylab("")+
  theme_bw()+theme(legend.position = "none")

ggsave("Fig4-CI.pdf",pcibox,height=4,width=6,scale=0.75)


## point vs. median

res.wide.mat <- res.wide.df[,c(5,8,11,7,9,10)]
names(res.wide.mat) <- rev(c("QSAR (von Borries)",
                                    "QSAR (Kvasnicka)",
                                    "QSAR (Wignall)",
                                    "IVIVE (ToxCast)",
                                    "ToxValDB (Aurisano)",
                                    "qRAx (PPRTV)"))
res.wide.mat$UF <- res.wide.df[,5]/res.wide.df[,12]
res.wide.mat$UF.bin <- round(10^(round(2*log10((res.wide.df[,5]/res.wide.df[,12])))/2),-1)
res.wide.mat$UF.bin[res.wide.mat$UF.bin==0] <- 1
res.wide.mat$UF.bin <- factor(res.wide.mat$UF.bin)

lowerfun <- function(data,mapping) {
  ggplot(data = data, mapping = mapping)+
    geom_point()+annotation_logticks(linewidth=0.25)+
    scale_x_log10(limits=c(1e-5,1e3))+scale_y_log10(limits=c(1e-5,1e3))+
    geom_abline(slope=1,intercept = c(-1,0,1),color="grey",
                linetype=c("dotted","dashed","dotted"))
}
pscat.med <- 
  ggpairs(res.wide.mat,columns = 1:6,aes(color=UF.bin),
          lower = list(continuous = wrap(lowerfun)),
          diag = list(continuous = ggally_blankDiag)
          )+
  scale_color_viridis_d(option="magma",end=0.8)+
    ggtitle("Human Equivalent Dose POD (mg/kg-d) [point estimate]")+
    xlab("mg/kg-d")+ylab("mg/kg-d")+
    theme_bw()+theme(panel.grid = element_blank())
ggsave("Fig4-A-point-median.pdf",pscat.med,height=6,width=6,scale=1.4)

## lb vs. lb

res.wide.lb.mat <- res.wide.df[,7+c(5,8,11,7,9,10)]
names(res.wide.lb.mat) <- rev(c("QSAR (von Borries)",
                             "QSAR (Kvasnicka)",
                             "QSAR (Wignall)",
                             "IVIVE (ToxCast)",
                             "ToxValDB (Aurisano)",
                             "qRAx (PPRTV)"))
res.wide.lb.mat$UF <- res.wide.df[,5]/res.wide.df[,12]
res.wide.lb.mat$UF.bin <- round(10^(round(2*log10((res.wide.df[,5]/res.wide.df[,12])))/2),-1)
res.wide.lb.mat$UF.bin[res.wide.lb.mat$UF.bin==0] <- 1
res.wide.lb.mat$UF.bin <- factor(res.wide.lb.mat$UF.bin)

lowerfun.lb <- function(data,mapping) {
  ggplot(data = data, mapping = mapping)+
    geom_point()+annotation_logticks(linewidth=0.25)+
    scale_x_log10(limits=c(1e-6,1e2))+scale_y_log10(limits=c(1e-6,1e2))+
    geom_abline(slope=1,intercept = c(-1,0,1),color="grey",
                linetype=c("dotted","dashed","dotted"))
}
pscat.lb <-
  ggpairs(res.wide.lb.mat,columns = 1:6,aes(color=UF.bin),
          lower = list(continuous = wrap(lowerfun.lb)),
          diag = list(continuous = ggally_blankDiag)
          )+
  scale_color_viridis_d(option="magma",end=0.8)+
  ggtitle("Human Equivalent Dose POD (mg/kg-d) [lower confidence bound]")+
    xlab("mg/kg-d")+ylab("mg/kg-d")+
    theme_bw()+theme(panel.grid = element_blank())
ggsave("Fig4-B-lb-lb.pdf",pscat.lb,height=6,width=6,scale=1.4)


## PPRTV lb vs.med

res.wide.lb.med.mat <- res.wide.df[c(12,8,11,7,9,10)]
names(res.wide.lb.med.mat) <- rev(c("QSAR (von Borries)",
                                "QSAR (Kvasnicka)",
                                "QSAR (Wignall)",
                                "IVIVE (ToxCast)",
                                "ToxValDB (Aurisano)",
                                "qRAx (PPRTV)"))
res.wide.lb.med.mat$UF <- res.wide.df[,5]/res.wide.df[,12]
res.wide.lb.med.mat$UF.bin <- round(10^(round(2*log10((res.wide.df[,5]/res.wide.df[,12])))/2),-1)
res.wide.lb.med.mat$UF.bin[res.wide.lb.med.mat$UF.bin==0] <- 1
res.wide.lb.med.mat$UF.bin <- factor(res.wide.lb.med.mat$UF.bin)

lowerfun.lb.med <- function(data,mapping) {
  ggplot(data = data, mapping = mapping)+
    geom_point()+annotation_logticks(linewidth=0.25)+
    scale_x_log10(limits=c(1e-6,1e3))+scale_y_log10(limits=c(1e-6,1e3))+
    geom_abline(slope=1,intercept = c(-1,0,1),color="grey",
                linetype=c("dotted","dashed","dotted"))
}
pscat.lb.med<-
  ggpairs(res.wide.lb.med.mat,columns = 1:6,aes(color=UF.bin),
          lower = list(continuous = wrap(lowerfun.lb.med)),
          diag = list(continuous = ggally_blankDiag)
          )+
  scale_color_viridis_d(option="magma",end=0.8)+
  ggtitle("Human Equivalent Dose POD (mg/kg-d) [PPRTV lb vs. point]")+
    xlab("mg/kg-d")+ylab("mg/kg-d")+
    theme_bw()+theme(panel.grid = element_blank())
ggsave("Fig4-C-lb-med.pdf",pscat.lb.med,height=6,width=6,scale=1.4)

#### Within 10-fold of qRAx

res.fold <- pivot_longer(log10(res.wide.mat[,2:6]/res.wide.mat[,1]),cols=1:5)
res.fold$name <- factor(res.fold$name,
                        levels=levels(res.df$name)[-6])
res.fold.lb <- pivot_longer(log10(res.wide.lb.mat[,2:6]/res.wide.lb.mat[,1]),cols=1:5)
res.fold.lb$name <- factor(res.fold.lb$name,
                        levels=levels(res.df$name)[-6])

pfold<-ggplot(res.fold)+
  annotate("rect", ymin = -Inf, ymax = Inf, 
           xmin = -1, xmax = 1, fill = 'grey',alpha=0.33)+
  geom_boxplot(aes(y=name,x=value,fill=name)) +
    scale_fill_viridis_d(direction = -1, end=0.9,alpha = 0.75)+
    ggtitle("NAM PODs vs. qRAx (PPRTV)")+
    scale_x_continuous(limits=c(-2,5),breaks=-2:5,
                       labels=math_format(10^.x))+
    xlab("Fold-difference")+
    annotation_logticks(side="b",outside=TRUE,
                        short=unit(0.02,"cm"),
                        mid=unit(0.05,"cm"),
                        long=unit(0.1,"cm"))+
    geom_vline(xintercept = c(-1,0,1),linetype=c("dotted","dashed","dotted"))+
    coord_cartesian(clip = "off")+
    ylab("")+
    theme_bw()+theme(legend.position = "none")

pfold.lb<-ggplot(res.fold.lb)+
  annotate("rect", ymin = -Inf, ymax = Inf, 
           xmin = -1, xmax = 1, fill = 'grey',alpha=0.33)+
    geom_boxplot(aes(y=name,x=value,fill=name)) +
    scale_fill_viridis_d(direction = -1, end=0.9,alpha = 0.75)+
    ggtitle("NAM PODs vs. qRAx [lower bounds]")+
    scale_x_continuous(limits=c(-2,5),breaks=-2:5,
                       labels=math_format(10^.x))+
    xlab("Fold-difference")+
    annotation_logticks(side="b",outside=TRUE,
                        short=unit(0.02,"cm"),
                        mid=unit(0.05,"cm"),
                        long=unit(0.1,"cm"))+
  geom_vline(xintercept = c(-1,0,1),linetype=c("dotted","dashed","dotted"))+
  coord_cartesian(clip = "off")+
    ylab("")+
    theme_bw()+theme(legend.position = "none")

ggsave("Fig4-fold.pdf",pfold,height=4,width=6,scale=0.75)
ggsave("Fig4-fold.lb.pdf",pfold.lb,height=4,width=6,scale=0.75)


#### Combinations

## coverage

res.suc.mat <- !is.na(res.wide.mat[,-(7:8)])

# pairs
combinations.2 <- as.data.frame(t(combn(6,2)))
combinations.2$NAM.1 <- colnames(res.suc.mat)[combinations.2$V1]
combinations.2$NAM.2 <- colnames(res.suc.mat)[combinations.2$V2]
combinations.2$NAM.comb <- 
  paste0(combinations.2$NAM.1,
         "+\n",
         combinations.2$NAM.2)
combinations.2$NAM.comb <- factor(combinations.2$NAM.comb,
                                  levels=rev(combinations.2$NAM.comb))
combinations.2$nsuccess <- 0
for (j in 1:nrow(combinations.2)) {
  combinations.2$nsuccess[j]<-
    sum(apply(res.suc.mat[,as.numeric(combinations.2[j,1:2])],1,sum)>0)
}
combinations.2$pct <- combinations.2$nsuccess/41
psuccess.comb.2 <-
  ggplot(combinations.2)+
  geom_col(aes(x=NAM.comb,y=nsuccess,fill=NAM.comb))+
  geom_label(aes(x=NAM.comb,y=nsuccess,label=paste0(round(100*nsuccess/41,1),"%")),hjust=1)+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Coverage from pairs")+
  xlab("")+ylab("n success (%)")+geom_hline(yintercept = 41,linetype="dashed")
print(psuccess.comb.2)

ggsave("Fig5A-pairs.pdf",psuccess.comb.2,height=4,width=6,scale=1.25)

# trios
combinations.3 <- as.data.frame(t(combn(6,3)))
combinations.3$NAM.1 <- colnames(res.suc.mat)[combinations.3$V1]
combinations.3$NAM.2 <- colnames(res.suc.mat)[combinations.3$V2]
combinations.3$NAM.3 <- colnames(res.suc.mat)[combinations.3$V3]
combinations.3$NAM.comb <- 
  paste0(combinations.3$NAM.1,
         "+\n",
         combinations.3$NAM.2,
         "+\n",
         combinations.3$NAM.3)
combinations.3$NAM.comb <- factor(combinations.3$NAM.comb,
                                  levels=rev(combinations.3$NAM.comb))
combinations.3$nsuccess <- 0
for (j in 1:nrow(combinations.3)) {
  combinations.3$nsuccess[j]<-
    sum(apply(res.suc.mat[,as.numeric(combinations.3[j,1:3])],1,sum)>0)
}
combinations.3$pct <- combinations.3$nsuccess/41
psuccess.comb.3 <-
  ggplot(combinations.3)+
  geom_col(aes(x=NAM.comb,y=nsuccess,fill=NAM.comb))+
  geom_label(aes(x=NAM.comb,y=nsuccess,label=paste0(round(100*nsuccess/41,1),"%")),hjust=1)+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Coverage from three")+
  xlab("")+ylab("n success (%)")+geom_hline(yintercept = 41,linetype="dashed")
print(psuccess.comb.3)

ggsave("Fig5B-three.pdf",psuccess.comb.3,height=6,width=6,scale=1.5)

#### Comparison to PPRTV

res.ratio.df <- pivot_longer(cbind(res.wide.df[,1:4],
                                   res.wide.df[,c(8,11,7,9,10)]/
                                     res.wide.df[,5]),cols=5:9)
res.ratio.df$name <- gsub("pod.med_","",res.ratio.df$name)
res.ratio.df$name[res.ratio.df$name=="PPRTV"] <- "qRAx (PPRTV)"
res.ratio.df$name[res.ratio.df$name=="Aurisano"] <- "ToxValDB (Aurisano)"
res.ratio.df$name[res.ratio.df$name=="ToxCast"] <- "IVIVE (ToxCast)"
res.ratio.df$name[res.ratio.df$name=="CTV"] <- "QSAR (Wignall)"
res.ratio.df$name[res.ratio.df$name=="Kvasnicka"] <- "QSAR (Kvasnicka)"
res.ratio.df$name[res.ratio.df$name=="von.Borries"] <- "QSAR (von Borries)"
res.ratio.df$name <- factor(res.ratio.df$name,
                           levels=c("QSAR (von Borries)",
                                    "QSAR (Kvasnicka)",
                                    "QSAR (Wignall)",
                                    "IVIVE (ToxCast)",
                                    "ToxValDB (Aurisano)",
                                    "qRAx (PPRTV)"))
ggplot(res.ratio.df)+geom_boxplot(aes(x=name,y=value))+
  geom_hline(yintercept=1,linetype="dashed")+
  scale_y_log10()+coord_flip()

res.ratio.lb.df <- pivot_longer(cbind(res.wide.df[,1:4],
                                   res.wide.df[,7+c(8,11,7,9,10)]/
                                     res.wide.df[,7+5]),cols=5:9)
res.ratio.lb.df$name <- gsub("pod.med.lb_","",res.ratio.lb.df$name)
res.ratio.lb.df$name[res.ratio.lb.df$name=="PPRTV"] <- "qRAx (PPRTV)"
res.ratio.lb.df$name[res.ratio.lb.df$name=="Aurisano"] <- "ToxValDB (Aurisano)"
res.ratio.lb.df$name[res.ratio.lb.df$name=="ToxCast"] <- "IVIVE (ToxCast)"
res.ratio.lb.df$name[res.ratio.lb.df$name=="CTV"] <- "QSAR (Wignall)"
res.ratio.lb.df$name[res.ratio.lb.df$name=="Kvasnicka"] <- "QSAR (Kvasnicka)"
res.ratio.lb.df$name[res.ratio.lb.df$name=="von.Borries"] <- "QSAR (von Borries)"
res.ratio.lb.df$name <- factor(res.ratio.lb.df$name,
                            levels=c("QSAR (von Borries)",
                                     "QSAR (Kvasnicka)",
                                     "QSAR (Wignall)",
                                     "IVIVE (ToxCast)",
                                     "ToxValDB (Aurisano)",
                                     "qRAx (PPRTV)"))
ggplot(res.ratio.lb.df)+geom_boxplot(aes(x=name,y=value))+
  geom_hline(yintercept=1,linetype="dashed")+
  scale_y_log10()+coord_flip()

## min POD combinations...

res.val.mat <- res.wide.mat[,-(7:8)]
apply(res.val.mat[,as.numeric(combinations.2[j,1:2])],1,min,na.rm=T)
