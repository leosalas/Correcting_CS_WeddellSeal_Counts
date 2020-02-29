# TODO: Add comment
# 
# Author: lsalas
###############################################################################


library(ggplot2);library(plyr);library(car);library(XLConnect)

## We use the estimates obtained in file calcNumberOfSeals_featureBased.R
load(file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/colonyEstimates.RData")

filterHigherSurveyQuality<-function(df,pctcut=0.2){
	colos<-unique(df$Colony)
	resdf<-ldply(.data=colos,.fun=function(x,df,pctcut){
				tdf<-subset(df,Colony==x);
				tdf<-tdf[order(tdf$survQuality),];
				nrec<-nrow(tdf)
				cutrec<-round(nrec*pctcut);
				if((nrec-cutrec)>2){
					tdf<-tdf[c(cutrec:nrec),]
				}else{tdf<-data.frame()}
				return(tdf)
			},df=df,pctcut=pctcut)
	return(resdf)
}

#####################################

#Estimates from tags
tagsOnlybyTag<-res_tagsOnlybyTag$colEstimates

exmpdf<-tagsOnlybyTag
exmpkey<-paste(exmpdf$Colony,exmpdf$acquisition_date,sep="::")
tagskey<-paste(tagStats$Colony,tagStats$ImgDate,sep="::")
nummiss<-tagskey[which(!tagskey %in% exmpkey)]
nn<-exmpkey[which(!exmpkey %in% tagskey)]; print(nn)
if(NROW(nn)==0 & NROW(nummiss)==0){
	print("All maps with tags have seal estimates!")
}

#load the counts
gCounts<-try(readWorksheetFromFile("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/ColonyData/MapsInColonies/GroundCounts_fromJRotella.xlsx",sheet="forML_colonyEstimates_v2 219012"))
names(gCounts)<-c("acquisition_date","Colony","jayCount")
gCounts$Colony<-as.character(gCounts$Colony)
gCounts$jayCount<-as.integer(as.character(gCounts$jayCount))
gCounts<-subset(gCounts,jayCount>=0 & !is.na(jayCount))

numSeals<-ldply(.data=list(tagsOnlybyTag),.fun=function(x,gCounts){	#83 satellite estimates are not linking to ground counts
			x$acquisition_date<-as.POSIXct(x$acquisition_date)
			tdf<-merge(x[,c("acquisition_date","Colony","catalogId","estNumSeals","numMaps","numTags","mapsSurveyed","corrMethod")],gCounts,by=c("acquisition_date","Colony"),all.y=T)
			return(tdf)
		},gCounts=gCounts)

load("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/dataFromGeoJson3.RData")

numSeals$diffMapsSurveyed<-numSeals$numMaps/numSeals$mapsSurveyed
numSeals$acYear<-format(numSeals$acquisition_date,"%Y")
numSeals$acMonth<-format(numSeals$acquisition_date,"%m")
numSeals$acDay<-format(numSeals$acquisition_date,"%d")
numSeals$acHour<-format(numSeals$acquisition_date,"%H")
numSeals$acMYH<-paste0(numSeals$acMonth,numSeals$acYear,numSeals$acHour)
numSeals$acYMD<-paste0(numSeals$acYear,numSeals$acMonth,numSeals$acDay)
numSeals<-merge(numSeals,tagStats,by.x=c("Colony","acquisition_date"),by.y=c("Colony","ImgDate"),all.x=TRUE)
numSeals$avgTagsPerMap<-ceiling(numSeals$NumTags/numSeals$NumMapsWtags)
names(numSeals)<-gsub("NumTags","AllTaggersTotalNumTags",names(numSeals))
names(numSeals)<-gsub("NumMapsWtags","AllTaggersNumMapsWtags",names(numSeals))

## For reference in Supplemental Information - here's the starting point to correct for the effects below...
pstart<-ggplot(data=numSeals,aes(x=jayCount,y=estNumSeals)) + geom_point() + 
		geom_smooth(method="lm",formula=y~x-1,color="red",se=F) + theme_bw() +
		labs(x="Ground count", y="Estimated number of seals")

jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/FigureS1_Reference.jpg",width=900,height=900,res=300,quality=100)
print(pstart)
dev.off()


## Add attribution about image acquisition date and sensor type
overlays$catalogId<-as.character(overlays$catalogId)
overlays$url<-as.character(overlays$url)
overlays$external_reference<-as.character(overlays$external_reference)
overlays$external_reference<-ifelse(!is.na(overlays$catalogId),overlays$catalogId,ifelse(grepl("amazonaws",overlays$external_reference,fixed=T),substr(overlays$external_reference,41,56),
				ifelse(grepl("cloudfront",overlays$external_reference),substr(overlays$external_reference,38,53),overlays$external_reference)))
overlays$satId<-ifelse((overlays$satId=="<unknown>" & substr(overlays$external_reference,1,3)=="103"),"WV02",
		ifelse((overlays$satId=="<unknown>" & substr(overlays$external_reference,1,3)=="102"),"WV01",overlays$satId))

overs<-data.frame()
for(ee in unique(overlays$external_reference)){
	tdf<-subset(overlays,external_reference==ee,select=c("overlayId","acquisition_date","satId","external_reference"))
	if(nrow(tdf)>1){
		dd<-min(tdf$acquisition_date,na.rm=TRUE)
		tdf$acquisition_date<-dd
	}
	overs<-rbind(overs,tdf)
}
nrow(unique(overs[,c("acquisition_date","satId","external_reference")]))==NROW(unique(overlays$external_reference))
names(overs)<-gsub("external_reference","catalogId",names(overs))
numSeals<-merge(numSeals,unique(overs[,c("catalogId","acquisition_date","satId")]),by=c("catalogId","acquisition_date"),all.x=TRUE)
numSeals$satId<-ifelse(!is.na(numSeals$satId),numSeals$satId,ifelse(substr(numSeals$catalogId,1,3)=="101","QB02",ifelse(substr(numSeals$catalogId,1,3)=="102","WV01","WV02")))
sum(is.na(numSeals$satId))
numSeals$f3<-substr(numSeals$catalogId,1,3)
unique(numSeals[,c("satId","f3")])

######################
## Aggregating maps by colony - omitting QB02, see below why
t1a<-aggregate(acquisition_date~Colony+acYear,data=subset(numSeals, satId!="QB02"),FUN=function(x){NROW(unique(x))})
t1b<-aggregate(jayCount~Colony+acYear,data=subset(numSeals, satId!="QB02"),FUN=function(x){NROW(unique(x))})
t1<-merge(t1a,t1b,by=c("Colony","acYear"))
names(t1)<-c("hauloutLocation","Year","numberImages","groundCounts")
write.csv(t1,file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Table3.csv")

## The problem with QB02
numSeals$seasPeriod<-ifelse(as.numeric(numSeals$acDay)<16,"Early Nov","Late Nov")
ordflm<-data.frame(acHour=sort(unique(numSeals$acHour)),ordvallm=c(4,5,6,7,8,9,10,11,1,2,3))
numSeals<-merge(numSeals,ordflm,by="acHour",all.x=T)
numSeals$acHour<-reorder(numSeals$acHour,numSeals$ordvallm)
numSeals$scaledNumTags<-as.numeric(scale(numSeals$numTags))
numSeals$scaledNumMaps<-as.numeric(scale(numSeals$numMaps))
numSeals$detRate<-numSeals$estNumSeals/numSeals$jayCount
mdl<-lm(detRate~scaledNumTags*scaledNumMaps+satId+acHour+seasPeriod+acYear,data=numSeals)
mdlsum<-summary(mdl)

## QB02 is borderline very different from WV02. Indeed, data from QB02 are really different from the other satellites...
ordf<-data.frame(acHour=sort(unique(numSeals$acHour)),localHour=c(13,14,15,16,17,18,8,9,10,11,12))
numSeals<-merge(numSeals,ordf,by="acHour",all.x=T)
numSeals$acHour<-reorder(numSeals$acHour,numSeals$localHour)
pSensor<-ggplot(data=numSeals,aes(x=detRate)) + geom_histogram(binwidth=0.1) + facet_wrap(~paste("Sensor:",satId),ncol=1) +
		scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8, 1.0),limits=c(-0.05,1.05)) +
		theme_bw() + labs(x="Detection rate",y="Frequency")

#For reference
sapply(unique(numSeals$satId),FUN=function(x,df){
			print(sd(subset(df,satId==x)$detRate))
		},df=numSeals)
jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/FigureS1_Reference.jpg",width=600,height=1200,res=300,quality=100)
	print(pSensor)
dev.off()


### The effect of hour on detection rates
## Fitting a sinusoidal to hour
numSeals$numHour<-((numSeals$localHour-15) %% 12)/12  #We tested various moduli, and 12 was best
sindf<-unique(numSeals[,c("localHour","numHour")])
sindf$sinH<-sin(2*pi*sindf$numHour)

meansdf<-aggregate(detRate~localHour,data=numSeals,FUN=mean);names(meansdf)<-c("localHour","meanRate")
meansdf<-subset(meansdf,localHour!=14)	#There is only one record for this hour, so not calculating the mean

psin<-ggplot(data=subset(numSeals,satId!="QB02"),aes(x=localHour,y=detRate)) + geom_jitter(width = 0.1) + 	#& acHour!="14"
		geom_line(data=sindf,aes(x=localHour,y=(0.2 + (sinH/8))),stat="identity") +
		geom_point(data=meansdf,aes(y=meanRate),size=2,color="red") + scale_x_continuous(breaks=c(8:18)) + 
		labs(x="Hour (GMT+13)", y="Detection rate") + theme_bw() +
		annotate("text", x = 18, y = 0.7, label = "B",size=5)
print(psin)
## But not including in paper...

## Refitting the model...
numSeals<-merge(numSeals,sindf[,c("localHour","sinH")],by="localHour",all.x=T)
numSealsF<-subset(numSeals,satId!="QB02" & acHour!="14")	#omitting QuickBird data, and the single observation for Hour 14
mdl<-lm(detRate~scaledNumTags*scaledNumMaps+satId+sinH+I(sinH^2)+seasPeriod+acYear,data=numSealsF)

## Goodness of fit...
numSealsF$predicted<-fitted(mdl); numSealsF$residuals<-residuals(mdl)
qqPlot(mdl)
p<-ggplot(data=numSealsF,aes(x=detRate,y=predicted)) + geom_point(size=1.3) + 
		geom_abline(intercept=0,slope=1,color="blue",size=1.2) + 
		labs(x="Observed detection rate",y="Predicted detection rate")

## We are still underpredicting at very low and high detection rates, adding colony effects and keeping only meaningful and influential covariates
mdlCol<-lm(detRate~scaledNumTags*scaledNumMaps+sinH+I(sinH^2)+Colony+acYear,data=numSealsF)	
summary(mdlCol)
mdlsum<-summary(mdlCol)
mdlTable<-as.data.frame(mdlsum$coef)
mdlTable$Variable<-row.names(mdlTable)
mdlTable<-mdlTable[,c(5,1:4)]
row.names(mdlTable)<-NULL
write.csv(mdlTable,file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Table4.csv")

## Goodness of fit...
layout(matrix(c(1,2,3,4),2,2)) 
plot(mdlCol) #FigureS2_Reference

numSealsF$predicted<-fitted(mdlCol); numSealsF$residuals<-residuals(mdlCol)
pF5A<-ggplot(data=numSealsF,aes(x=detRate,y=predicted)) + geom_point(size=1.3) + 
		geom_abline(intercept=0,slope=1,color="blue",size=1.2) + 
		theme_bw() +
		annotate("text", x = 0.65, y = 0.48, label = "A",size=5) +
		labs(x="Observed detection rate",y="Predicted detection rate")

jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure4A.jpg",width=900,height=900,res=300,quality=100)
	print(pF5A)
dev.off()


df<-numSealsF
df<-df[order(abs(df$residuals)),]
tail(df)
## Much better - final model for Colony!
## LRT test for a model without colony effects

## Modeling colony effects as Island vs Mainlad effects.
island.df<-data.frame(Colony=c("Hutton Cliffs","Turtle Rock","Turks Head Tryggve","Big Razorback","Little Razorback","Tent Island","Inaccessible Island","South Base"),
		Island=c(0,1,0,1,1,1,1,0))
numSealsF<-merge(numSealsF,island.df,by="Colony",all.x=T)
mdlIsl<-lm(detRate~scaledNumTags*scaledNumMaps+sinH+I(sinH^2)+acYear+Island,data=numSealsF)
summary(mdlIsl)
#Simplifying...
mdlIsl<-lm(detRate~scaledNumTags+sinH+I(sinH^2)+Island+acYear,data=numSealsF)
summary(mdlIsl)
## THIS IS the best model without colony effects. It predicts the ratio by which we inflate estNumSeals to obtain the correct number.
## Not the best fit, but will have to do...

## Creating table of coefficients and GOF plot
mdlsum<-summary(mdlIsl)
mdlTable<-as.data.frame(mdlsum$coef)
mdlTable$Variable<-row.names(mdlTable)
mdlTable<-mdlTable[,c(5,1:4)]
row.names(mdlTable)<-NULL
write.csv(mdlTable,file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Table5.csv")

## Goodness of fit...
layout(matrix(c(1,2,3,4),2,2)) 
plot(mdlIsl) # FigureS3_Reference


numSealsI<-numSealsF;numSealsI$predicted<-fitted(mdlIsl); numSealsI$residuals<-residuals(mdlIsl)
pF5B<-ggplot(data=numSealsI,aes(x=detRate,y=predicted)) + geom_point(size=1.3) + 
		geom_abline(intercept=0,slope=1,color="blue",size=1.2) + 
		theme_bw() +
		annotate("text", x = 0.65, y = 0.4, label = "B",size=5) +
		labs(x="Observed detection rate",y="Predicted detection rate")

jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure4B.jpg",width=900,height=900,res=300,quality=100)
	print(pF5B)
dev.off()

## Now we move to MethodsPaper_results3.R for the estimation of seal numbers using both these models
save(numSeals,numSealsF,numSealsI,mdlCol,mdlIsl, file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/finalModelsAndData.RData")



