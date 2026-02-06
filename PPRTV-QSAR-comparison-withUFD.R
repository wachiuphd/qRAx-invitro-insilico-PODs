library(tidyverse)
library(readxl)
library(ctxR)
datafolder <- "data"
resultsfolder <- "results"

## Get DTXSID from Comptox Dashboard for each PPRTV
my_key <- '01cbab76-904f-11ee-92b7-325096b39f47'
register_ctx_api_key(key=my_key)
pprtv <- read_excel(file.path(datafolder,"03282025 PPRTV with Read-across-Oral-PODs.xlsx"),sheet=2)
pprtv.CASRN.search <- chemical_equal_batch(word_list = (pprtv$CASRN))
pprtv.CASRN.search.details<-get_chemical_details_batch(DTXSID = pprtv.CASRN.search$dtxsid)
pprtv.CASRN.search.valid <- left_join(pprtv.CASRN.search,
                                      pprtv.CASRN.search.details)
names(pprtv.CASRN.search.valid)[names(pprtv.CASRN.search.valid)=="casrn"]<-"CASRN"
names(pprtv.CASRN.search.valid)[names(pprtv.CASRN.search.valid)=="dtxsid"]<-"DTXSID"
names(pprtv.CASRN.search.valid)[names(pprtv.CASRN.search.valid)=="inchikey"]<-"INCHIKEY"
pprtv.CASRN.search.valid <- pprtv.CASRN.search.valid[,-"descriptorStringTsv"]
pprtv <- left_join(pprtv,pprtv.CASRN.search.valid)
#write.csv(pprtv,"PPRTV-all-substances.csv",row.names=FALSE)
#pprtv <- subset(pprtv,!is.na(DTXSID) & !is.na(POD))
pprtv$`POD-HED.lb` <- pprtv$`POD-HED`/(pprtv$UF_L*pprtv$UF_S*pprtv$UF_D)
# write.csv(pprtv,"PPRTV-for-comparison-withUFD.csv",row.names = FALSE)

## ToxCast assays
toxcast.dat <- pprtv
toxcast.dat$AC50.05 <- NA
toxcast.dat$AC50.50 <- NA
toxcast.dat$AC50.active.05 <- NA
toxcast.dat$AC50.active.50 <- NA
toxcast.dat$HED.05 <- NA
toxcast.dat$HED.50 <- NA
toxcast.dat$HED.active.05 <- NA
toxcast.dat$HED.active.50 <- NA

for (j in 1:nrow(toxcast.dat)) {
  tmp <- get_bioactivity_details(toxcast.dat$DTXSID[j])
  tmp.httk <- get_httk_data(toxcast.dat$DTXSID[j])
  if (length(tmp)>0 & length(tmp.httk>0)) {
    toxcast.dat$AC50.05[j] <- quantile(tmp$ac50,0.05,na.rm=T)
    toxcast.dat$AC50.50[j] <- quantile(tmp$ac50,0.5,na.rm=T)
    toxcast.dat$AC50.active.05[j] <- quantile(subset(tmp,hitcall>0.95)$ac50,0.05,na.rm=T)
    toxcast.dat$AC50.active.50[j] <- quantile(subset(tmp,hitcall>0.95)$ac50,0.5,na.rm=T)
    css <- subset(tmp.httk,parameter=="Css" & model=="3compartmentss" & percentile=="95%" & species=="Human")$predicted
    css.uM <- 1000 * css / toxcast.dat$monoisotopicMass
    toxcast.dat$HED.05[j] <- toxcast.dat$AC50.05[j] / css.uM
    toxcast.dat$HED.50[j] <- toxcast.dat$AC50.50[j] / css.uM
    toxcast.dat$HED.active.05[j] <- toxcast.dat$AC50.active.05[j] / css.uM
    toxcast.dat$HED.active.50[j] <- toxcast.dat$AC50.active.50[j] / css.uM
  }
}

cat("ToxCast missing",sum(is.na(toxcast.dat$HED.active.05)),"out of",nrow(pprtv),"\n")
ggplot(toxcast.dat)+geom_point(aes(x=`POD-HED`,y=HED.active.50))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=HED.active.50))+
  geom_errorbar(aes(x=`POD-HED`,ymax=HED.active.50,ymin=HED.active.05))+
  scale_x_log10(limits=c(1e-5,100))+scale_y_log10(limits=c(1e-5,100))+
  ggtitle("PPRTV vs. ToxCast")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
# ggsave("PPRTV.ToxCast-withUFD.pdf",height=6,width=6)
summary(lm(log10(HED.active.50)~log10(`POD-HED`),data=toxcast.dat))
summary(lm(log10(HED.active.50)-log10(`POD-HED`)~0,data=toxcast.dat))
summary(lm(log10(HED.active.50)-log10(`POD-HED`)~1,data=toxcast.dat))
summary(lm(log10(HED.active.05)~log10(`POD-HED`),data=toxcast.dat))
summary(lm(log10(HED.active.05)-log10(`POD-HED`)~0,data=toxcast.dat))
summary(lm(log10(HED.active.05)-log10(`POD-HED`)~1,data=toxcast.dat))
summary(lm(log10(HED.active.05)~log10(`POD-HED.lb`),data=toxcast.dat))
summary(lm(log10(HED.active.05)-log10(`POD-HED.lb`)~0,data=toxcast.dat))
summary(lm(log10(HED.active.05)-log10(`POD-HED.lb`)~1,data=toxcast.dat))

## CTV
ctv.input <- data.frame(CAS=pprtv$CASRN,
                        SMILES=pprtv$qsarReadySmiles,
                        Name.Check=pprtv$Chemical,
                        AVERAGE_MASS=pprtv$monoisotopicMass)
write.csv(ctv.input,file.path(datafolder,"PPRTV-CTV-input.csv",row.names=FALSE))
ctv.dat <- as_tibble(read.csv(file.path(datafolder,"PPRTV-CTV_ToxValuePredictions.csv")))
ctv.dat <- rename(ctv.dat,qsarReadySmiles = SMILES)
ctv.dat <- rename(ctv.dat,CASRN = CAS)
ctv.dat <- subset(ctv.dat,tv %in% c("RfDNOAEL","RfDBMDL"))
ctv.dat <- subset(ctv.dat,!is.na(appl.domain)) # has appl domain
ctv.dat <- subset(ctv.dat,appl.domain < 3)
ctv.dat <- pivot_wider(ctv.dat,names_from=tv,values_from = 6:14)
# Tox Values - assume rat to convert to HED
ctv.dat$HED_RfDBMDL <- ctv.dat$prediction_RfDBMDL * 0.25
ctv.dat$HED_RfDNOAEL <- ctv.dat$prediction_RfDNOAEL * 0.25
ctv.dat$HED_RfDBMDL.lb <- ctv.dat$lower95_RfDBMDL * 0.25
ctv.dat$HED_RfDNOAEL.lb <- ctv.dat$lower95_RfDNOAEL * 0.25
ctv.dat$HED_RfDBMDL.ub <- ctv.dat$upper95_RfDBMDL * 0.25
ctv.dat$HED_RfDNOAEL.ub <- ctv.dat$upper95_RfDNOAEL * 0.25
ctv.dat$HED_min <- pmin(ctv.dat$HED_RfDNOAEL,ctv.dat$HED_RfDBMDL,na.rm=T)
ctv.dat$HED_min.lb <- pmin(ctv.dat$HED_RfDNOAEL.lb,ctv.dat$HED_RfDBMDL.lb,na.rm=T)

ctv.dat <- left_join(pprtv,ctv.dat)

cat("CTV missing",sum(is.na(ctv.dat$HED_min)),"out of",nrow(pprtv),"\n")
ggplot(ctv.dat)+geom_point(aes(x=`POD-HED`,y=HED_min))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=HED_min))+
  geom_errorbar(aes(x=`POD-HED`,ymax=HED_min,ymin=HED_min.lb))+
  scale_x_log10(limits=c(1e-5,100))+scale_y_log10(limits=c(1e-5,100))+
  ggtitle("PPRTV vs. CTV")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave("PPRTV.CTV-withUFD.pdf",height=6,width=6)
summary(lm(log10(HED_min)~log10(`POD-HED`),data=ctv.dat))
summary(lm(log10(HED_min)-log10(`POD-HED`)~0,data=ctv.dat))
summary(lm(log10(HED_min)-log10(`POD-HED`)~1,data=ctv.dat))
summary(lm(log10(HED_min.lb)~log10(`POD-HED`),data=ctv.dat))
summary(lm(log10(HED_min.lb)-log10(`POD-HED`)~0,data=ctv.dat))
summary(lm(log10(HED_min.lb)-log10(`POD-HED`)~1,data=ctv.dat))
summary(lm(log10(HED_min.lb)~log10(`POD-HED.lb`),data=ctv.dat))
summary(lm(log10(HED_min.lb)-log10(`POD-HED.lb`)~0,data=ctv.dat))
summary(lm(log10(HED_min.lb)-log10(`POD-HED.lb`)~1,data=ctv.dat))

## Aurisano
aurisano.dat <- as_tibble(read.csv("Aurisano_ExcelTableS5.csv"))
aurisano.dat <- rename(aurisano.dat,DTXSID = dtxsid)
aurisano.dat <- rename(aurisano.dat,HEDgen = PODgen)
# Tox Values
aurisano.dat$HEDgen.lb <- aurisano.dat$HEDgen*aurisano.dat$PODgen.GSD2final^(-1.65/2)
aurisano.dat$HEDgen.ub <- aurisano.dat$HEDgen*aurisano.dat$PODgen.GSD2final^(1.65/2)
aurisano.dat <- rename(aurisano.dat,HEDrd = PODrd)
aurisano.dat$HEDrd.lb <- aurisano.dat$HEDrd*aurisano.dat$PODrd.GSD2final^(-1.65/2)
aurisano.dat$HEDrd.ub <- aurisano.dat$HEDrd*aurisano.dat$PODrd.GSD2final^(1.65/2)
aurisano.dat$HEDmin <- pmin(aurisano.dat$HEDgen,aurisano.dat$HEDrd,na.rm=T)
aurisano.dat$HEDmin.lb <- pmin(aurisano.dat$HEDgen.lb,aurisano.dat$HEDrd.lb,na.rm=T)

aurisano.dat <- left_join(pprtv,aurisano.dat)
cat("Aurisano missing",sum(is.na(aurisano.dat$HEDmin)),"out of",nrow(pprtv),"\n")
ggplot(aurisano.dat)+geom_point(aes(x=`POD-HED`,y=HEDmin))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=HEDmin))+
  geom_errorbar(aes(x=`POD-HED`,ymax=HEDmin,ymin=HEDmin.lb))+
  scale_x_log10(limits=c(1e-5,100))+scale_y_log10(limits=c(1e-5,100))+
  ggtitle("PPRTV vs. Aurisano")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave("PPRTV.Aurisano-withUFD.pdf",height=6,width=6)
summary(lm(log10(HEDmin)~log10(`POD-HED`),data=aurisano.dat))
summary(lm(log10(HEDmin)-log10(`POD-HED`)~0,data=aurisano.dat))
summary(lm(log10(HEDmin)-log10(`POD-HED`)~1,data=aurisano.dat))
summary(lm(log10(HEDmin.lb)~log10(`POD-HED`),data=aurisano.dat))
summary(lm(log10(HEDmin.lb)-log10(`POD-HED`)~0,data=aurisano.dat))
summary(lm(log10(HEDmin.lb)-log10(`POD-HED`)~1,data=aurisano.dat))
summary(lm(log10(HEDmin.lb)~log10(`POD-HED.lb`),data=aurisano.dat))
summary(lm(log10(HEDmin.lb)-log10(`POD-HED.lb`)~0,data=aurisano.dat))
summary(lm(log10(HEDmin.lb)-log10(`POD-HED.lb`)~1,data=aurisano.dat))

# load kvasnicka
kvasnicka.dat <- as_tibble(read.csv(file.path("Two-Stage-ML-Results-Browser",
                                              "PODs800k.csv")))
kvasnicka.dat <- subset(kvasnicka.dat, DTXSID %in% pprtv$DTXSID)
kvasnicka.dat <- rename(kvasnicka.dat, CASRN = CAS_RN_Dashboard)
kvasnicka.dat$pod.min <- pmin(kvasnicka.dat$pod.gen,kvasnicka.dat$pod.repdev,na.rm=T)
kvasnicka.dat$pod.min.lb <- pmin(kvasnicka.dat$lb.gen,kvasnicka.dat$lb.repdev,na.rm=T)
kvasnicka.dat <- left_join(pprtv,kvasnicka.dat)
cat("Kvasnicka missing",sum(is.na(kvasnicka.dat$pod.min.lb)),"out of",nrow(pprtv),"\n")
ggplot(kvasnicka.dat)+geom_point(aes(x=`POD-HED`,y=pod.min))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=pod.min))+
  geom_errorbar(aes(x=`POD-HED`,ymax=pod.min,ymin=pod.min.lb))+
  scale_x_log10(limits=c(1e-5,100))+scale_y_log10(limits=c(1e-5,100))+
  ggtitle("PPRTV vs. Kvasnicka")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave("PPRTV.Kvasnicka-withUFD.pdf",height=6,width=6)
summary(lm(log10(pod.min)~log10(`POD-HED`),data=kvasnicka.dat))
summary(lm(log10(pod.min)-log10(`POD-HED`)~0,data=kvasnicka.dat))
summary(lm(log10(pod.min)-log10(`POD-HED`)~1,data=kvasnicka.dat))
summary(lm(log10(pod.min.lb)~log10(`POD-HED`),data=kvasnicka.dat))
summary(lm(log10(pod.min.lb)-log10(`POD-HED`)~0,data=kvasnicka.dat))
summary(lm(log10(pod.min.lb)-log10(`POD-HED`)~1,data=kvasnicka.dat))
summary(lm(log10(pod.min.lb)~log10(`POD-HED.lb`),data=kvasnicka.dat))
summary(lm(log10(pod.min.lb)-log10(`POD-HED.lb`)~0,data=kvasnicka.dat))
summary(lm(log10(pod.min.lb)-log10(`POD-HED.lb`)~1,data=kvasnicka.dat))

## vonBorries
vonBorries.dat <- as_tibble(read.csv("vonBorries.POD.all.csv"))
vonBorries.dat <- left_join(pprtv,vonBorries.dat,by="INCHIKEY")
vonBorries.dat <- vonBorries.dat[,-which("CASRN.y"==names(vonBorries.dat))]
vonBorries.dat <- rename(vonBorries.dat,CASRN = CASRN.x)
# load externally predicted values
vonBorries.dat.2 <- as_tibble(read.csv("PODUAM_results-2.csv"))
vonBorries.dat.2 <- rename(vonBorries.dat.2,qsarReadySmiles = Original.SMILES)
indx <- which(vonBorries.dat$qsarReadySmiles %in% vonBorries.dat.2$qsarReadySmiles |
                vonBorries.dat$msReadySmiles %in% vonBorries.dat.2$qsarReadySmiles)
indx.2 <- which(vonBorries.dat.2$qsarReadySmiles %in% vonBorries.dat$qsarReadySmiles |
                  vonBorries.dat.2$qsarReadySmiles %in% vonBorries.dat$msReadySmiles)
vonBorries.dat$POD.nc.pred[indx] <- vonBorries.dat.2$Predicted.PODnc..log10.mg.kg.d.[indx.2]
vonBorries.dat$POD.nc.pred.lb[indx] <- vonBorries.dat.2$PODnc..2.5..[indx.2]
vonBorries.dat$POD.nc.pred.ub[indx] <- vonBorries.dat.2$PODnc..97.5..[indx.2]
vonBorries.dat$POD.rd.pred[indx] <- vonBorries.dat.2$Predicted.PODrd..log10.mg.kg.d.[indx.2]
vonBorries.dat$POD.rd.pred.lb[indx] <- vonBorries.dat.2$PODrd..2.5..[indx.2]
vonBorries.dat$POD.rd.pred.ub[indx] <- vonBorries.dat.2$PODrd..97.5..[indx.2]

vonBorries.dat$HED.nc <- 10^vonBorries.dat$POD.nc.pred
vonBorries.dat$HED.nc.lb <- 10^vonBorries.dat$POD.nc.pred.lb
vonBorries.dat$HED.nc.ub <- 10^vonBorries.dat$POD.nc.pred.ub
vonBorries.dat$HED.rd<- 10^vonBorries.dat$POD.rd.pred
vonBorries.dat$HED.rd.lb <- 10^vonBorries.dat$POD.rd.pred.lb
vonBorries.dat$HED.rd.ub <- 10^vonBorries.dat$POD.rd.pred.ub
vonBorries.dat$HED.min <- pmin(vonBorries.dat$HED.nc,vonBorries.dat$HED.rd,na.rm=T)
vonBorries.dat$HED.min.lb <- pmin(vonBorries.dat$HED.nc.lb,vonBorries.dat$HED.rd.lb,na.rm=T)

cat("von Borries missing",sum(is.na(vonBorries.dat$HED.min.lb)),"out of",nrow(pprtv),"\n")
ggplot(vonBorries.dat)+geom_point(aes(x=`POD-HED`,y=HED.min))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=HED.min))+
  geom_errorbar(aes(x=`POD-HED`,ymax=HED.min,ymin=HED.min.lb))+
  scale_x_log10(limits=c(1e-5,100))+scale_y_log10(limits=c(1e-5,100))+
  ggtitle("PPRTV vs. VonBorries")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave("PPRTV.VonBorries-withUFD.pdf",height=6,width=6)
summary(lm(log10(HED.min)~log10(`POD-HED`),data=vonBorries.dat))
summary(lm(log10(HED.min)-log10(`POD-HED`)~0,data=vonBorries.dat))
summary(lm(log10(HED.min)-log10(`POD-HED`)~1,data=vonBorries.dat))
summary(lm(log10(HED.min.lb)~log10(`POD-HED`),data=vonBorries.dat))
summary(lm(log10(HED.min.lb)-log10(`POD-HED`)~0,data=vonBorries.dat))
summary(lm(log10(HED.min.lb)-log10(`POD-HED`)~1,data=vonBorries.dat))
summary(lm(log10(HED.min.lb)~log10(`POD-HED.lb`),data=vonBorries.dat))
summary(lm(log10(HED.min.lb)-log10(`POD-HED.lb`)~0,data=vonBorries.dat))
summary(lm(log10(HED.min.lb)-log10(`POD-HED.lb`)~1,data=vonBorries.dat))

## All together - PODs

res.df <- pprtv[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
res.df$name <- "PPRTV"
res.df$pod.med <- pprtv$`POD-HED`
res.df$pod.med.lb <- pprtv$`POD-HED.lb`
# TTC
tmp.df <- ttc.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "TTC"
tmp.df$pod.med <- ttc.dat$TTC_HEDmin
tmp.df$pod.med.lb <- NA
res.df <- rbind(res.df,tmp.df)
# CTV
tmp.df <- ctv.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "CTV"
tmp.df$pod.med <- ctv.dat$HED_min
tmp.df$pod.med.lb <- ctv.dat$HED_min.lb
res.df <- rbind(res.df,tmp.df)
# Aurisano
tmp.df <- aurisano.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "Aurisano"
tmp.df$pod.med <- aurisano.dat$HEDmin
tmp.df$pod.med.lb <- aurisano.dat$HEDmin.lb
res.df <- rbind(res.df,tmp.df)
# Kvasnicka
tmp.df <- aurisano.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "Kvasnicka"
tmp.df$pod.med <- kvasnicka.dat$pod.min
tmp.df$pod.med.lb <- kvasnicka.dat$pod.min.lb
res.df <- rbind(res.df,tmp.df)
# von Borries
tmp.df <- vonBorries.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "von Borries"
tmp.df$pod.med <- vonBorries.dat$HED.min
tmp.df$pod.med.lb <- vonBorries.dat$HED.min.lb
res.df <- rbind(res.df,tmp.df)
# ToxCast
tmp.df <- toxcast.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "ToxCast"
tmp.df$pod.med <- toxcast.dat$HED.active.50
tmp.df$pod.med.lb <- toxcast.dat$HED.active.05
res.df <- rbind(res.df,tmp.df)

res.df$name <- factor(res.df$name,levels=c("PPRTV","ToxCast","von Borries","Kvasnicka","CTV","Aurisano","TTC"))

ggplot(res.df)+geom_point(aes(x=pod.med,y=name,color=name))+
  geom_errorbarh(aes(xmin=pod.med.lb,xmax=pod.med,y=name,color=name))+
  facet_wrap(~CASRN)+scale_x_log10()
ggsave("PODs.bychem-withUFD.pdf",height=6,width=10)
write.csv(res.df,"PODs.bychem-withUFD.csv",row.names = FALSE)

res.wide.df <- pivot_wider(res.df,names_from=name,values_from = 6:7)

write.csv(res.wide.df,"PODs.bychem.source-withUFD.csv",row.names = FALSE)

## All together - ratios

res.ratio.df <- pprtv[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
res.ratio.df$name <- "PPRTV"
res.ratio.df$ratio.med <- median(pprtv$`POD-HED`)/pprtv$`POD-HED`
res.ratio.df$ratio.med.lb <- median(pprtv$`POD-HED`)/pprtv$`POD-HED.lb`
res.ratio.df$ratio.med.lb.lb <- median(pprtv$`POD-HED.lb`)/pprtv$`POD-HED.lb`
# TTC
tmp.df <- ttc.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "TTC"
tmp.df$ratio.med <- ttc.dat$TTC_HEDmin/ttc.dat$`POD-HED`
tmp.df$ratio.med.lb <- ttc.dat$TTC_HEDmin/ttc.dat$`POD-HED.lb`
tmp.df$ratio.med.lb.lb <- ttc.dat$TTC_HEDmin/ttc.dat$`POD-HED.lb`
res.ratio.df <- rbind(res.ratio.df,tmp.df)
# CTV
tmp.df <- ctv.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "CTV"
tmp.df$ratio.med <- ctv.dat$HED_min/ctv.dat$`POD-HED`
tmp.df$ratio.med.lb <- ctv.dat$HED_min/ctv.dat$`POD-HED.lb`
tmp.df$ratio.med.lb.lb <- ctv.dat$HED_min.lb/ctv.dat$`POD-HED.lb`
res.ratio.df <- rbind(res.ratio.df,tmp.df)
# Aurisano
tmp.df <- aurisano.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "Aurisano"
tmp.df$ratio.med <- aurisano.dat$HEDmin/aurisano.dat$`POD-HED`
tmp.df$ratio.med.lb <- aurisano.dat$HEDmin/aurisano.dat$`POD-HED.lb`
tmp.df$ratio.med.lb.lb <- aurisano.dat$HEDmin.lb/aurisano.dat$`POD-HED.lb`
res.ratio.df <- rbind(res.ratio.df,tmp.df)
# Kvasnicka
tmp.df <- aurisano.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "Kvasnicka"
tmp.df$ratio.med <- kvasnicka.dat$pod.min/kvasnicka.dat$`POD-HED`
tmp.df$ratio.med.lb <- kvasnicka.dat$pod.min/kvasnicka.dat$`POD-HED.lb`
tmp.df$ratio.med.lb.lb <- kvasnicka.dat$pod.min.lb/kvasnicka.dat$`POD-HED.lb`
res.ratio.df <- rbind(res.ratio.df,tmp.df)
# von Borries
tmp.df <- vonBorries.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "von Borries"
tmp.df$ratio.med <- vonBorries.dat$HED.min/vonBorries.dat$`POD-HED`
tmp.df$ratio.med.lb <- vonBorries.dat$HED.min/vonBorries.dat$`POD-HED.lb`
tmp.df$ratio.med.lb.lb <- vonBorries.dat$HED.min.lb/vonBorries.dat$`POD-HED.lb`
res.ratio.df <- rbind(res.ratio.df,tmp.df)
# ToxCast
tmp.df <- toxcast.dat[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
tmp.df$name <- "ToxCast"
tmp.df$ratio.med <- toxcast.dat$HED.active.50/toxcast.dat$`POD-HED`
tmp.df$ratio.med.lb <- toxcast.dat$HED.active.50/vonBorries.dat$`POD-HED.lb`
tmp.df$ratio.med.lb.lb <- toxcast.dat$HED.active.05/vonBorries.dat$`POD-HED.lb`
res.ratio.df <- rbind(res.ratio.df,tmp.df)

res.ratio.df$name <- factor(res.ratio.df$name,levels=c("ToxCast","von Borries","Kvasnicka","CTV","Aurisano","TTC","PPRTV"))

ggplot(res.ratio.df)+geom_boxplot(aes(x=ratio.med,y=name))+scale_x_log10(limits=c(0.5e-3,2e3))+annotation_logticks(sides="b")
ggplot(res.ratio.df)+geom_boxplot(aes(x=ratio.med.lb,y=name))+scale_x_log10(limits=c(0.5e-3,2e3))+annotation_logticks(sides="b")
ggplot(res.ratio.df)+geom_boxplot(aes(x=ratio.med.lb.lb,y=name))+scale_x_log10(limits=c(0.5e-3,2e3))+annotation_logticks(sides="b")
ggplot(res.ratio.df)+geom_histogram(aes(x=ratio.med))+scale_x_log10(limits=c(0.5e-3,2e3))+facet_wrap(~name,ncol=1)+annotation_logticks(sides="b")
ggplot(res.ratio.df)+geom_histogram(aes(x=ratio.med.lb))+scale_x_log10(limits=c(0.5e-3,2e3))+facet_wrap(~name,ncol=1)+annotation_logticks(sides="b")
ggplot(res.ratio.df)+geom_histogram(aes(x=ratio.med.lb.lb))+scale_x_log10(limits=c(0.5e-3,2e3))+facet_wrap(~name,ncol=1)+annotation_logticks(sides="b")
