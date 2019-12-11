# TODO: Add comment
# 
# Author: lsalas
###############################################################################


library(ggplot2);library(car);library(fitdistrplus)

load(file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/finalModelsAndData.RData")

## Calculating seal numbers using the Island model
qdf<-numSealsF
preds<-predict(mdlIsl,se.fit=T)
qdf$predRatio<-preds$fit;qdf$se.predRatio<-preds$se.fit
qdf$predRatio<-ifelse(qdf$predRatio<0,0,qdf$predRatio)

## Can only aggregate regionally by image aqcquisition date, and relate these to the nearest ground count.
## See this 
qdf$acYMDH<-format(qdf$acquisition_date,"%Y%m%d%H%M")

coldf<-data.frame()
for(q in unique(qdf$acYMDH)){
	nc<-NROW(unique(subset(qdf,acYMDH==q)$Colony))
	ymd<-substr(q,1,8)
	sat<-paste(unique(subset(qdf,acYMDH==q)$satId),collapse="::")
	tdf<-data.frame(YMDH=q,day=ymd,satId=sat,NumColonies=nc)
	coldf<-rbind(coldf,tdf)
}

## The regional counts
ground<-aggregate(jayCount~acYMDH,data=qdf,FUN=sum); names(ground)<-c("acYMDH","groundCount")
## The estimate, average ratio and average SE of the ratio
estN<-aggregate(estNumSeals~acYMDH,data=qdf,FUN=sum)
rat<-aggregate(predRatio~acYMDH,data=qdf,FUN=mean)
serat<-aggregate(se.predRatio~acYMDH,data=qdf,FUN=mean)

rdf<-merge(ground,rat,by="acYMDH")
rdf<-merge(rdf,serat,by="acYMDH")
rdf<-merge(rdf,estN,by="acYMDH")

rdf$predCount<-round(rdf$estNumSeals/rdf$predRatio)
rdf$predCount<-ifelse(rdf$predCount<0,0,rdf$predCount)

#Adding +/- SE values
rdf$pcu<-rdf$estNumSeals/(rdf$predRatio+(1.96*rdf$se.predRatio))		
rdf$pcl<-rdf$estNumSeals/(rdf$predRatio-(1.96*rdf$se.predRatio))
rdf$pcl<-ifelse(rdf$pcl<0,0,rdf$pcl)

pRegion<-ggplot(data=rdf,aes(x=groundCount,y=predCount)) + 
		geom_abline(slope=1,intercept=1,linetype="dotted",color="black", size=1.3) +
		geom_smooth(method="lm",formula=y~x-1,color="dark gray",se=F) + geom_point() +
		geom_errorbar(aes(ymin=pcl,ymax=pcu)) +
		labs(x="Ground count",y="Predicted abundance") + theme_bw() 

jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure6.jpg",width=900,height=900,res=300,quality=100)
	print(pRegion)
dev.off()

