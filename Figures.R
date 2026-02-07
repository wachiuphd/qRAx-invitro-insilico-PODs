library(tidyverse)
library(GGally)
library(scales)
resultsfolder <- "results"
res.wide.df <- read.csv(file.path(resultsfolder,"PODs.bychem.source-withUFD.csv"))
res.wide.df$Chemical <- factor(res.wide.df$Chemical,
                            levels=res.wide.df$Chemical[order(res.wide.df$Chemical)])
res.wide.df$CASRN <- factor(res.wide.df$CASRN,
                       levels=res.wide.df$CASRN[order(res.wide.df$Chemical)])

res.df <- read.csv(file.path(resultsfolder,"PODs.bychem-withUFD.csv"))
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
ggsave(file.path(resultsfolder,"Fig3A-success.pdf"),psuccess,height=4,width=6,scale=0.75)

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
ggsave(file.path(resultsfolder,"Fig2-PODs.bychem.pdf"),ppods,height=10,width=6,scale=1.5)

## Most sensitive - median estimates
varnames <- c(pod.med_von.Borries="QSAR (von Borries)",
              pod.med_Kvasnicka="QSAR (Kvasnicka)",
              pod.med_CTV="QSAR (Wignall)",
              pod.med_ToxCast="IVIVE (ToxCast)",
              pod.med_Aurisano="ToxValDB (Aurisano)",
              pod.med_PPRTV="qRAx (PPRTV)")
min.pods <- data.frame(
  var=names(res.wide.df[,5:10]),
  value=as.numeric(table(apply(res.wide.df[,5:10],1,which.min)))
) 
min.pods$name <- factor(as.character(varnames[min.pods$var]),
                      levels=as.character(varnames))
pmin <-
  ggplot(min.pods)+
  geom_col(aes(x=name,y=value,fill=name))+
  geom_label(aes(x=name,y=value,label=paste0(round(100*value/41),"%")))+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Minimum PODs")+
  xlab("")+ylab("Frequency out of 41 (%)")
print(pmin)
ggsave(file.path(resultsfolder,"Fig3B-mostsens.med.pdf"),pmin,height=4,width=6,scale=0.75)

## Lower bounds
varnames.lb <- c(pod.med.lb_von.Borries="QSAR (von Borries)",
              pod.med.lb_Kvasnicka="QSAR (Kvasnicka)",
              pod.med.lb_CTV="QSAR (Wignall)",
              pod.med.lb_ToxCast="IVIVE (ToxCast)",
              pod.med.lb_Aurisano="ToxValDB (Aurisano)",
              pod.med.lb_PPRTV="qRAx (PPRTV)")
min.pods.lb <- data.frame(
  var=names(res.wide.df[,11:16]),
  value=as.numeric(table(apply(res.wide.df[,11:16],1,which.min))[paste(1:6)])
) 
min.pods.lb$value[is.na(min.pods.lb$value)]<-0
min.pods.lb$name <- factor(as.character(varnames.lb[min.pods.lb$var]),
                        levels=as.character(varnames.lb))
pmin.lb <-
  ggplot(min.pods.lb)+
  geom_col(aes(x=name,y=value,fill=name))+
  geom_label(aes(x=name,y=value,label=paste0(round(100*value/41),"%")))+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Minimum PODs (lower bound)")+
  xlab("")+ylab("Frequency out of 41 (%)")
print(pmin.lb)
ggsave(file.path(resultsfolder,"FigS3A-mostsens.lb.pdf"),pmin.lb,height=4,width=6,scale=0.75)

## Lower bound PPRTV with median for others
varnames.lb.med <- c(pod.med_von.Borries="QSAR (von Borries)",
                 pod.med_Kvasnicka="QSAR (Kvasnicka)",
                 pod.med_CTV="QSAR (Wignall)",
                 pod.med_ToxCast="IVIVE (ToxCast)",
                 pod.med_Aurisano="ToxValDB (Aurisano)",
                 pod.med.lb_PPRTV="qRAx (PPRTV)")
min.pods.lb.med <- data.frame(
  var=names(res.wide.df[,c(6:11)]),
  value=as.numeric(table(apply(res.wide.df[,c(6:11)],1,which.min))[paste(1:6)])
) 
min.pods.lb.med$value[is.na(min.pods.lb.med$value)]<-0
min.pods.lb.med$name <- factor(as.character(varnames.lb.med[min.pods.lb.med$var]),
                           levels=as.character(varnames.lb.med))
pmin.lb.med <-
  ggplot(min.pods.lb.med)+
  geom_col(aes(x=name,y=value,fill=name))+
  geom_label(aes(x=name,y=value,label=paste0(round(100*value/41),"%")))+
  coord_flip()+scale_fill_viridis_d(direction = -1,end=0.9)+
  theme_bw()+theme(legend.position = "none")+
  ggtitle("Minimum PODs (PPRTV lower bound)")+
  xlab("")+ylab("Frequency out of 41 (%)")
print(pmin.lb.med)
ggsave(file.path(resultsfolder,"FigS3B-mostsens.lb.med.pdf"),pmin.lb.med,height=4,width=6,scale=0.75)


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

ggsave(file.path(resultsfolder,"Fig3C-CI.pdf"),pcibox,height=4,width=6,scale=0.75)

## point vs. median

res.wide.mat <- res.wide.df[,rev(names(varnames))]
names(res.wide.mat) <- varnames[rev(names(varnames))]
res.wide.mat$UF <- res.wide.df$pod.med_PPRTV/res.wide.df$pod.med.lb_PPRTV
res.wide.mat$UF.bin <- round(10^(round(2*log10(res.wide.mat$UF))/2),-1)
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
          upper = list(continuous = wrap("cor", method = "spearman")),
          lower = list(continuous = wrap(lowerfun)),
          diag = list(continuous = ggally_blankDiag)
          )+
  scale_color_viridis_d(option="magma",end=0.8)+
    ggtitle("Human Equivalent Dose POD (mg/kg-d) [point estimate]")+
    xlab("mg/kg-d")+ylab("mg/kg-d")+
    theme_bw()+theme(panel.grid = element_blank())
ggsave(file.path(resultsfolder,"Fig5-point-median.pdf"),pscat.med,height=6,width=6,scale=1.4)

## lb vs. lb
res.wide.lb.mat <- res.wide.df[,rev(names(varnames.lb))]
names(res.wide.lb.mat) <- varnames.lb[rev(names(varnames.lb))]
res.wide.lb.mat$UF <- res.wide.df$pod.med_PPRTV/res.wide.df$pod.med.lb_PPRTV
res.wide.lb.mat$UF.bin <- round(10^(round(2*log10(res.wide.lb.mat$UF))/2),-1)
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
          upper = list(continuous = wrap("cor", method = "spearman")),
          lower = list(continuous = wrap(lowerfun.lb)),
          diag = list(continuous = ggally_blankDiag)
          )+
  scale_color_viridis_d(option="magma",end=0.8)+
  ggtitle("Human Equivalent Dose POD (mg/kg-d) [lower confidence bound]")+
    xlab("mg/kg-d")+ylab("mg/kg-d")+
    theme_bw()+theme(panel.grid = element_blank())
ggsave(file.path(resultsfolder,"FigS4-lb-lb.pdf"),pscat.lb,height=6,width=6,scale=1.4)


## PPRTV lb vs.med
res.wide.lb.med.mat <- res.wide.df[,rev(names(varnames.lb.med))]
names(res.wide.lb.med.mat) <- varnames.lb.med[rev(names(varnames.lb.med))]
res.wide.lb.med.mat$UF <- res.wide.df$pod.med_PPRTV/res.wide.df$pod.med.lb_PPRTV
res.wide.lb.med.mat$UF.bin <- round(10^(round(2*log10(res.wide.lb.med.mat$UF))/2),-1)
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
          upper = list(continuous = wrap("cor", method = "spearman")),
          lower = list(continuous = wrap(lowerfun.lb.med)),
          diag = list(continuous = ggally_blankDiag)
          )+
  scale_color_viridis_d(option="magma",end=0.8)+
  ggtitle("Human Equivalent Dose POD (mg/kg-d) [PPRTV lb vs. point]")+
    xlab("mg/kg-d")+ylab("mg/kg-d")+
    theme_bw()+theme(panel.grid = element_blank())
ggsave(file.path(resultsfolder,"FigS5-lb-med.pdf"),pscat.lb.med,height=6,width=6,scale=1.4)

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

ggsave(file.path(resultsfolder,"Fig4A-fold.pdf"),pfold,height=4,width=6,scale=0.75)
ggsave(file.path(resultsfolder,"Fig4B-fold.lb.pdf"),pfold.lb,height=4,width=6,scale=0.75)


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

ggsave(file.path(resultsfolder,"FigS1-pairs.pdf"),psuccess.comb.2,height=4,width=6,scale=1.25)

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

ggsave(file.path(resultsfolder,"FigS2-three.pdf"),psuccess.comb.3,height=6,width=6,scale=1.5)

