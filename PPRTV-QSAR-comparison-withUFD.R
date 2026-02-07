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
pprtv$`POD-HED.lb` <- pprtv$`POD-HED`/(pprtv$UF_L*pprtv$UF_S*pprtv$UF_D)

## ToxCast assays - use ctxr package to get bioactivity data
## use median and p05 AC50 for the "median" and "lb" in vitro POD in uM
## Use httk to convert from uM to mg/kg-d POD (using 95th percentile Css)

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
  print(paste(j,toxcast.dat$DTXSID[j],"getting bioactivity..."))
  tmp <- get_bioactivity_details(toxcast.dat$DTXSID[j])
  print("getting httk...")
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
ggsave(file.path(resultsfolder,"PPRTV.ToxCast-withUFD.pdf"),height=6,width=6)

## CTV - Wignall et al. 2018. The "input" file was run through the batch version
## of CTV to generate the "output"
## applicability domain z-score must be <3

ctv.input <- data.frame(CAS=pprtv$CASRN,
                        SMILES=pprtv$smiles,
                        Name.Check=pprtv$Chemical,
                        AVERAGE_MASS=pprtv$monoisotopicMass)
write.csv(ctv.input,file.path(datafolder,"PPRTV-CTV-input.csv"),row.names=FALSE)
ctv.dat <- as_tibble(read.csv(file.path(datafolder,"PPRTV-CTV_ToxValuePredictions.csv")))
ctv.dat <- rename(ctv.dat,smiles = SMILES)
ctv.dat <- rename(ctv.dat,CASRN = CAS)
ctv.dat <- subset(ctv.dat,tv %in% c("RfDNOAEL","RfDBMDL"))
ctv.dat <- subset(ctv.dat,!is.na(appl.domain)) # has appl domain
ctv.dat <- subset(ctv.dat,appl.domain < 3)
ctv.dat <- pivot_wider(ctv.dat,names_from=tv,values_from = 6:14)

# assume rat to convert predicted POD to HED - use more conservative of NOAEL and BMDL
ctv.dat$HED_RfDBMDL <- ctv.dat$prediction_RfDBMDL * 0.25
ctv.dat$HED_RfDNOAEL <- ctv.dat$prediction_RfDNOAEL * 0.25
ctv.dat$HED_RfDBMDL.lb <- ctv.dat$lower95_RfDBMDL * 0.25
ctv.dat$HED_RfDNOAEL.lb <- ctv.dat$lower95_RfDNOAEL * 0.25
ctv.dat$HED_RfDBMDL.ub <- ctv.dat$upper95_RfDBMDL * 0.25
ctv.dat$HED_RfDNOAEL.ub <- ctv.dat$upper95_RfDNOAEL * 0.25
ctv.dat$HED_min <- pmin(ctv.dat$HED_RfDNOAEL,ctv.dat$HED_RfDBMDL,na.rm=T)
ctv.dat$HED_min.lb <- pmin(ctv.dat$HED_RfDNOAEL.lb,ctv.dat$HED_RfDBMDL.lb,na.rm=T)

ctv.dat <- left_join(pprtv,ctv.dat,by="CASRN")

cat("CTV missing",sum(is.na(ctv.dat$HED_min)),"out of",nrow(pprtv),"\n")
ggplot(ctv.dat)+geom_point(aes(x=`POD-HED`,y=HED_min))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=HED_min))+
  geom_errorbar(aes(x=`POD-HED`,ymax=HED_min,ymin=HED_min.lb))+
  scale_x_log10(limits=c(1e-6,300))+scale_y_log10(limits=c(1e-6,300))+
  ggtitle("PPRTV vs. CTV")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave(file.path(resultsfolder,"PPRTV.CTV-withUFD.pdf"),height=6,width=6)

## Aurisano - obtained PODs from Aurisano et al. 2023, supplemental table
## Derive lb as p05 using the GSD from the table
## use more conservative of general non-cancer and rep/dev effects

aurisano.dat <- as_tibble(read.csv(file.path(datafolder,"Aurisano_ExcelTableS5.csv")))
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
  scale_x_log10(limits=c(1e-6,350))+scale_y_log10(limits=c(1e-6,350))+
  ggtitle("PPRTV vs. Aurisano")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave(file.path(resultsfolder,"PPRTV.Aurisano-withUFD.pdf"),height=6,width=6)

## load kvasnicka et al. 2024 precalculated values and intersect with PPRTVs
## use more conservative of general non-cancer and rep/dev effects

kvasnicka.dat <- as_tibble(read.csv(file.path(datafolder,
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
  scale_x_log10(limits=c(1e-6,350))+scale_y_log10(limits=c(1e-6,350))+
  ggtitle("PPRTV vs. Kvasnicka")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave(file.path(resultsfolder,"PPRTV.Kvasnicka-withUFD.pdf"),height=6,width=6)

## vonBorries et al. 2026
## extract PPRTV smiles and use as input to PODUAM web app at
## https://dtu-quantitative-sustainability-assessment.shinyapps.io/poduam/
## use more conservative of general non-cancer and rep/dev effects

vonBorries.input <-pprtv[,c("CASRN","smiles")]
write.csv(vonBorries.input,file.path(datafolder,"VonBorries-input.csv"))

# load externally predicted values
vonBorries.dat <- as_tibble(read.csv(file.path(datafolder,"PODUAM_results.csv")))
vonBorries.dat <- rename(vonBorries.dat,smiles = Original.SMILES)
vonBorries.dat <- rename(vonBorries.dat,POD.nc.pred = Predicted.PODnc..log10.mg.kg.d.)
vonBorries.dat <- rename(vonBorries.dat,POD.nc.pred.lb = PODnc..2.5..)
vonBorries.dat <- rename(vonBorries.dat,POD.nc.pred.ub = PODnc..97.5..)
vonBorries.dat <- rename(vonBorries.dat,POD.rd.pred = Predicted.PODrd..log10.mg.kg.d.)
vonBorries.dat <- rename(vonBorries.dat,POD.rd.pred.lb = PODrd..2.5..)
vonBorries.dat <- rename(vonBorries.dat,POD.rd.pred.ub = PODrd..97.5..)

vonBorries.dat$HED.nc <- 10^vonBorries.dat$POD.nc.pred
vonBorries.dat$HED.nc.lb <- 10^vonBorries.dat$POD.nc.pred.lb
vonBorries.dat$HED.nc.ub <- 10^vonBorries.dat$POD.nc.pred.ub
vonBorries.dat$HED.rd<- 10^vonBorries.dat$POD.rd.pred
vonBorries.dat$HED.rd.lb <- 10^vonBorries.dat$POD.rd.pred.lb
vonBorries.dat$HED.rd.ub <- 10^vonBorries.dat$POD.rd.pred.ub
vonBorries.dat$HED.min <- pmin(vonBorries.dat$HED.nc,vonBorries.dat$HED.rd,na.rm=T)
vonBorries.dat$HED.min.lb <- pmin(vonBorries.dat$HED.nc.lb,vonBorries.dat$HED.rd.lb,na.rm=T)

vonBorries.dat <- left_join(pprtv,vonBorries.dat,by="smiles")

cat("von Borries missing",sum(is.na(vonBorries.dat$HED.min.lb)),"out of",nrow(pprtv),"\n")
ggplot(vonBorries.dat)+geom_point(aes(x=`POD-HED`,y=HED.min))+
  geom_errorbarh(aes(xmin=`POD-HED.lb`,xmax=`POD-HED`,y=HED.min))+
  geom_errorbar(aes(x=`POD-HED`,ymax=HED.min,ymin=HED.min.lb))+
  scale_x_log10(limits=c(1e-6,350))+scale_y_log10(limits=c(1e-6,350))+
  ggtitle("PPRTV vs. VonBorries")+
  geom_abline(slope=1,intercept=c(0,-1,1),linetype=c("dashed","dotted","dotted"))
ggsave(file.path(resultsfolder,"PPRTV.VonBorries-withUFD.pdf"),height=6,width=6)

## Combine all PODs together in single data frame and save

res.df <- pprtv[,c("Chemical","CASRN","DTXSID","INCHIKEY")]
res.df$name <- "PPRTV"
res.df$pod.med <- pprtv$`POD-HED`
res.df$pod.med.lb <- pprtv$`POD-HED.lb`
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

res.df$name <- factor(res.df$name,levels=c("PPRTV","ToxCast","von Borries","Kvasnicka","CTV","Aurisano"))

ggplot(res.df)+geom_point(aes(x=pod.med,y=name,color=name))+
  geom_errorbarh(aes(xmin=pod.med.lb,xmax=pod.med,y=name,color=name))+
  facet_wrap(~CASRN)+scale_x_log10()
ggsave(file.path(resultsfolder,"PODs.bychem-withUFD.pdf"),height=6,width=10)
write.csv(res.df,file.path(resultsfolder,"PODs.bychem-withUFD.csv"),row.names = FALSE)

res.wide.df <- pivot_wider(res.df,names_from=name,values_from = 6:7)

write.csv(res.wide.df,file.path(resultsfolder,"PODs.bychem.source-withUFD.csv"),row.names = FALSE)
