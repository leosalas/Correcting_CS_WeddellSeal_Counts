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
qdf$predCount<-qdf$estNumSeals/qdf$predRatio
qdf$predCount<-ifelse(qdf$predCount<0,0,qdf$predCount)

#Adding +/- SE values
qdf$pcu<-qdf$estNumSeals/(qdf$predRatio+(1.96*qdf$se.predRatio))		
qdf$pcl<-qdf$estNumSeals/(qdf$predRatio-(1.96*qdf$se.predRatio))
qdf$pcl<-ifelse(qdf$pcl<0,0,qdf$pcl)

#let's fit a gamma to predRatio and use the 95% quantiles of the gamma to filter outlier values
fg<-try(fitdist(qdf$predRatio, distr="gamma",method="mme"),silent=TRUE)
plot(fg) 
print(fg$estimate)	#may need them below...
lcut<-qgamma(0.025,shape=fg$estimate["shape"],rate=fg$estimate["rate"])
ucut<-qgamma(0.975,shape=fg$estimate["shape"],rate=fg$estimate["rate"])
meandf_g<-aggregate(predCount~jayCount+Colony,subset(qdf,predRatio>=lcut & predRatio<=ucut),mean,na.rm=T);names(meandf_g)<-c("groundCount","Colony","meanPredCount")
seudf_g<-aggregate(pcu~jayCount+Colony,subset(qdf,predRatio>=lcut & predRatio<=ucut),mean,na.rm=T);names(seudf_g)<-c("groundCount","Colony","upperPredCount")
seldf_g<-aggregate(pcl~jayCount+Colony,subset(qdf,predRatio>=lcut & predRatio<=ucut),mean,na.rm=T);names(seldf_g)<-c("groundCount","Colony","lowerPredCount")

islplotdf<-merge(meandf_g,seudf_g,by=c("groundCount","Colony"),all.x=T)
islplotdf<-merge(islplotdf,seldf_g,by=c("groundCount","Colony"),all.x=T)
islplotdf$meanPredCount<-round(islplotdf$meanPredCount,1)
islplotdf$upperPredCount<-round(islplotdf$upperPredCount,1)
islplotdf$lowerPredCount<-round(islplotdf$lowerPredCount,1)

predIsland<-ggplot(data=islplotdf,aes(x=groundCount,y=meanPredCount)) +  
		geom_abline(slope=1,intercept=1,linetype="dotted",color="black", size=1.3) +
		geom_smooth(method="lm",formula=y~x-1,color="dark gray",se=F) + geom_point(aes(color=Colony)) +
		geom_errorbar(aes(ymin=lowerPredCount,ymax=upperPredCount,color=Colony)) +
		labs(x="Ground count",y="Predicted abundance",color="Location") + theme_bw() +
		annotate("text", x = 210, y = 250, label = "B",size=5)
print(predIsland)

jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure5B.jpg",width=1400,height=900,res=300,quality=100)
	print(predIsland)
dev.off()


## For the Appendix
mdlIslandCuts<-ggplot(data=qdf,aes(x=predRatio)) + geom_histogram(binwidth=0.05) + 
		geom_vline(xintercept=lcut,color="red") + geom_vline(xintercept=ucut,color="red") + 
		theme_bw() + annotate("text", x = 0.5, y = 36, label = "B",size=5) +
		scale_x_continuous(breaks=seq(0,1.4,by=0.2),limits=c(0,0.6)) + labs(x="Estimated detection rate",y="Frequency")
dev.new();print(mdlIslandCuts)
jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/FigureS4B_Reference.jpg",width=700,height=700,res=300,quality=100)
	print(mdlIslandCuts)
dev.off()


##################################################
## If we used the Colony model...
## Need to average detection rates for islands vs mainlands first, then multiply by jayCount
qdf<-numSealsF
preds<-predict(mdlCol,se.fit=T)
qdf$predRatio<-preds$fit;qdf$se.predRatio<-preds$se.fit

#let's fit a gamma to predRatio first...
fg<-try(fitdist(qdf$predRatio, distr="gamma",method="mme",start=list(shape=7.45,rate=41.9)),silent=TRUE)	#
plot(fg) 
lcut<-qgamma(0.025,shape=fg$estimate["shape"],rate=fg$estimate["rate"])
ucut<-qgamma(0.975,shape=fg$estimate["shape"],rate=fg$estimate["rate"])

mldpredRatio<-mean(subset(qdf,Island==0 & predRatio>=lcut & predRatio<=ucut)$predRatio);mldse.predRatio<-mean(subset(qdf,Island==0)$se.predRatio)
islpredRatio<-mean(subset(qdf,Island==1 & predRatio>=lcut & predRatio<=ucut)$predRatio);islse.predRatio<-mean(subset(qdf,Island==1)$se.predRatio)
avgprdf<-data.frame(Colony=c("Big Razorback","Inaccessible Island","Little Razorback","South Base","Tent Island","Turtle Rock","Hutton Cliffs","Turks Head Tryggve"),
		avgpredRatio=c(rep(mldpredRatio,6),rep(islpredRatio,2)),avgse.predRatio=c(rep(mldse.predRatio,6),rep(islse.predRatio,2)))

qdf<-merge(qdf,avgprdf,by=c("Colony"),all.x=T)

qdf$predCount<-qdf$estNumSeals/qdf$avgpredRatio
qdf$predCount<-ifelse(qdf$predCount<0,0,qdf$predCount)
qdf$pcu<-qdf$estNumSeals/(qdf$avgpredRatio+(1.96*qdf$avgse.predRatio))		
qdf$pcl<-qdf$estNumSeals/(qdf$avgpredRatio-(1.96*qdf$avgse.predRatio))
qdf$pcl<-ifelse(qdf$pcl<0,0,qdf$pcl)

meandf_g<-aggregate(predCount~jayCount+Colony,subset(qdf,predRatio>=lcut & predRatio<=ucut),mean,na.rm=T);names(meandf_g)<-c("groundCount","Colony","meanPredCount")
seudf_g<-aggregate(pcu~jayCount+Colony,subset(qdf,predRatio>=lcut & predRatio<=ucut),mean,na.rm=T);names(seudf_g)<-c("groundCount","Colony","upperPredCount")
seldf_g<-aggregate(pcl~jayCount+Colony,subset(qdf,predRatio>=lcut & predRatio<=ucut),mean,na.rm=T);names(seldf_g)<-c("groundCount","Colony","lowerPredCount")

locplotdf<-merge(meandf_g,seudf_g,by=c("groundCount","Colony"),all.x=T)
locplotdf<-merge(locplotdf,seldf_g,by=c("groundCount","Colony"),all.x=T)
locplotdf$meanPredCount<-round(locplotdf$meanPredCount,1)
locplotdf$upperPredCount<-round(locplotdf$upperPredCount,1)
locplotdf$lowerPredCount<-round(locplotdf$lowerPredCount,1)

predColony<-ggplot(data=locplotdf,aes(x=groundCount,y=meanPredCount)) +  
		geom_abline(slope=1,intercept=1,linetype="dotted",color="black", size=1.3) +
		geom_smooth(method="lm",formula=y~x-1,color="dark gray",se=F) + geom_point(aes(color=Colony)) +
		geom_errorbar(aes(ymin=lowerPredCount,ymax=upperPredCount,color=Colony)) +
		labs(x="Ground count",y="Predicted abundance",color="Location") + theme_bw() +
		annotate("text", x = 210, y = 200, label = "A",size=5)
print(predColony)
jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure5A.jpg",width=1400,height=900,res=300,quality=100)
	print(predColony)
dev.off()


## For the Appendix
mdlColonyCuts<-ggplot(data=qdf,aes(x=predRatio)) + geom_histogram(binwidth=0.05) + 
		geom_vline(xintercept=lcut,color="red") + geom_vline(xintercept=ucut,color="red") +
		theme_bw() + annotate("text", x = 0.5, y = 30, label = "A",size=5) +
		scale_x_continuous(breaks=seq(0,1.4,by=0.2),limits=c(0,0.6)) + labs(x="Estimated detection rate",y="Frequency")
dev.new();print(mdlColonyCuts)
jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/FigureS4A_Reference.jpg",width=700,height=700,res=300,quality=100)
	print(mdlColonyCuts)
dev.off()

