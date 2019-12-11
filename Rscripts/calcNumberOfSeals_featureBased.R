# TODO: 
# How are we estimating numMaps inspected by colony?? 
# Complete overlays$catalogId from the URL, use that to aggregate -
# Add satId and correct for sat effects. Of the 4 unknowns, those with catalogId = 103... are wv02, and the other is wv01
# 
# Author: lsalas
###############################################################################


## Load the features
## See in which maps there are features, and which of these were explored by Michelle
## ProbF is the probability of a feature in the collection of maps inspected by 2 or more people - universal
## ProbS is the probability of Michelle finding a seal in a map she inspected - universal
## ProbFS is ss/(ss+ns) - for the maps inspected by ML, is her tag within 3m of a feature? - depending on CR threshold
## ProbSF is ProbFS*ProbS/ProbF - universal
## Corr factor is ProbSF/ProbFS, or ProbS/ProbF, thus independent of threshold unless... 
## We calculate ProbFS by map, and thus ProbSF by map. This means that CorrFactor has expectation ProbS/ProbF. Let's see!

libs<-c("ggplot2","raster","rgdal","plyr","sp","rgeos","fitdistrplus")
lapply(libs, require, character.only = TRUE)

######################  FUNCTIONS WE'LL NEED
## Convert to UTM
# data is the table with the geo data
# latfield is the string naming the latitude field
# lonfield is the string naming the longitude field
convertToUTM<-function(data,latfield,lonfield){
	
	coordinates(data)<-c(lonfield,latfield)
	proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
	utmdata <- spTransform(data, CRS("+proj=utm +zone=58 +south ellps=WGS84"))
	utmdata<-as.data.frame(utmdata)
	names(utmdata)<-gsub(lonfield,"easting",names(utmdata))
	names(utmdata)<-gsub(latfield,"northing",names(utmdata))
	
	return(utmdata)
}

## Get the number of map by rank of tagger
# x is the filter rank value
# df is the table that lists the maps, taggers that inspected these, and their rank value
getNumMapsByRank<-function(x,df){
	mlMaps<-unique(subset(df,taggerId==21758509,select="mapId")$mapId)
	dfml<-subset(df,mapId %in% mlMaps)
	#keep only those where ML is not the sole inspector
	mlplus<-sapply(X=dfml$mapId,FUN=function(X,dfx){
				gg<-subset(dfx,mapId==X);
				mapchk<-ifelse(NROW(unique(gg$taggerId))>1,X,NA);
				return(mapchk)
			},dfx=dfml)
	dfml<-subset(dfml,mapId %in% na.omit(mlplus))
	
	dff<-subset(df,taggerScore>=x)
	numaps<-NROW(unique(dff$mapId))
	dffml<-subset(dfml,taggerScore>=x)
	numaps<-NROW(unique(dff$mapId));numlmaps<-NROW(unique(dffml$mapId))
	nutaggers<-NROW(unique(dff$taggerId));numltaggers<-NROW(unique(dffml$taggerId))
	resdf<-data.frame(threshold=x,Nmaps=numaps,NMLmaps=numlmaps,Ntaggers=nutaggers,NMLtaggers=numltaggers)
	return(resdf)
}

# x is the value to logit-transform to
toLogit<-function(x){
	cx<-ifelse(x==0,0.00001,
			ifelse(x==1,0.99999,x))
	lgx<-log(cx)-log(1-cx)
	return(lgx)
}

# x is the value to logit-transform from
fromLogit<-function(x){
	bt<-exp(x)/(1+exp(x))
	return(bt)
}

## Find the number of seals the tager correctly identified (compared to ML)
# taggerTags is the data.frame with all the tagger's tags in maps shared with ML, including the east/north of each tag
# mlTags is the data.frame with all of ML tags, including the east/north of each tag
getTaggerSS<-function(taggerTags,mlTags){
	#loop through each tag and see if it's within 3m of a ML tag
	ss<-0
	for(rr in 1:nrow(taggerTags)){
		mlTags$taggerEast<-taggerTags[rr,"easting"]
		mlTags$taggerNorth<-taggerTags[rr,"northing"]
		mlTags$dist<-sqrt(((mlTags$easting-mlTags$taggerEast)^2)+((mlTags$northing-mlTags$taggerNorth)^2))
		featTags<-subset(mlTags,dist <= 3)
		if(nrow(featTags)>0){
			#take the closest and remove from mlTags
			ss<-ss+1
			featTags<-featTags[order(featTags$dist),]
			topTag<-as.integer(featTags[1,"tagId"])
			mlTags<-subset(mlTags,tagId!=topTag)
			if(nrow(mlTags)==0){break}
		}
	}
	return(ss)
}

## The following two functions stimate the number of seals in each map in the data, then group by colony
## One does it based on features, and the other based on tags
# gspm is the table of estimates of the correction factor for all taggers at a particular threshold
# crtm is the table of features per map, with mapId attribution
# tgvutm is the table of tags, with mapId attribution
# maps and overlays are the original data tables listing the maps and overlays
# colinfo is the colony information table, linking each map to a colony and survey date
# corrMethod is a string naming the method used to calculate the correction factor
# distfilt is a value to filter the corrFactor values used to fit the gamma
getColonyEstimates_byFeature<-function(gspm,crtm,maps,overlays,colinfo,corrMethod){
	cutrows<-nrow(gspm)-ceiling(nrow(gspm)*0.05)
	gspm<-gspm[order(gspm$corrFactor),]
	gspm<-gspm[1:cutrows,]
	gammdist<-fitdist(gspm$corrFactor, "gamma", method="mle")$estimate
	#meanCorrF<-mean(subset(gspm,corrFactor<2)$corrFactor)
	meanCorrF<-qgamma(0.5,gammdist[1],gammdist[2])
	gdupper<-qgamma(0.975,gammdist[1],gammdist[2])
	gdlower<-qgamma(0.025,gammdist[1],gammdist[2])
	propUpper<-gdupper/meanCorrF;propLower<-gdlower/meanCorrF
	
	## Estimate number of seals per map
	numSeals<-aggregate(numFeatures~mapId,crtm,sum)
	numSeals$estNumSeals<-round(numSeals$numFeatures*meanCorrF)
	numSeals$uclNumSeals<-round(numSeals$estNumSeals*propUpper)
	numSeals$lclNumSeals<-round(numSeals$estNumSeals*propLower)
	
	numSeals<-merge(numSeals,maps[,c("mapId","overlayId")],by="mapId",all.x=T)
	numSeals<-merge(numSeals,overlays[,c("overlayId","acquisition_date","satId","external_reference")],by="overlayId",all.x=T)
	numSeals$catalogId<-ifelse(grepl("amazonaws",numSeals$external_reference),substr(numSeals$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",numSeals$external_reference),substr(numSeals$external_reference,38,53),as.character(numSeals$external_reference)))
	#numSeals<-subset(numSeals,!is.na(acquisition_date))	#we lose 435 maps that do not have acDate, from 4500 maps we down to 4065
	numSeals$acDate<-format(numSeals$acquisition_date,"%Y%m%d")
	numSeals$acYear<-as.integer(format(numSeals$acquisition_date,"%Y"))
	numSeals$acHour<-as.integer(format(numSeals$acquisition_date,"%H"))
	numSeals$acquisition_date<-as.character(numSeals$acquisition_date)
	
	numCol<-merge(numSeals,colinfo[,c("mapId","Colony")],by="mapId")	#Only 786 with info!!
	colNmaps<-aggregate(estNumSeals~acquisition_date+Colony+catalogId,data=numCol,NROW)
	names(colNmaps)<-gsub("estNumSeals","numMaps",names(colNmaps))
	colest<-aggregate(estNumSeals~acquisition_date+Colony+catalogId,data=numCol,sum)
	colfeat<-aggregate(numFeatures~acquisition_date+Colony+catalogId,data=numCol,sum)
	colupper<-aggregate(uclNumSeals~acquisition_date+Colony+catalogId,data=numCol,sum)
	collower<-aggregate(lclNumSeals~acquisition_date+Colony+catalogId,data=numCol,sum)
	
	colEstimates<-merge(colNmaps,colest,by=c("acquisition_date","Colony","catalogId"))
	colEstimates<-merge(colEstimates,colfeat,by=c("acquisition_date","Colony","catalogId"))
	colEstimates<-merge(colEstimates,colupper,by=c("acquisition_date","Colony","catalogId"))
	colEstimates<-merge(colEstimates,collower,by=c("acquisition_date","Colony","catalogId"))
	
	colsurveyed<-merge(colinfo,maps,by="mapId",all.x=T)
	colsurveyed<-merge(colsurveyed,overlays,by="overlayId")
	colsurveyed$acquisition_date<-as.character(colsurveyed$acquisition_date)
	colsurveyed$catalogId<-ifelse(grepl("amazonaws",colsurveyed$external_reference),substr(colsurveyed$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",colsurveyed$external_reference),substr(colsurveyed$external_reference,38,53),as.character(colsurveyed$external_reference)))
	colsurveyed$mapsSurveyed<-1
	numInColSurveyed<-aggregate(mapsSurveyed~Colony+acquisition_date+catalogId,data=colsurveyed,sum)
	
	colEstimates<-merge(colEstimates,numInColSurveyed,by=c("acquisition_date","Colony","catalogId"),all.x=T)
	colEstimates$corrMethod<-corrMethod
	reslst<-list(numSeals=numSeals,colEstimates=colEstimates,gammdist=gammdist)
	return(reslst)
}

getColonyEstimates_byTag<-function(gspm,tgvutm,maps,overlays,colinfo,corrMethod,distfilt){
	cutrows<-nrow(gspm)-ceiling(nrow(gspm)*0.05)
	gspm<-gspm[order(gspm$corrFactor),]
	gspm<-gspm[1:cutrows,]
	gammdist<-fitdist(gspm$corrFactor, "gamma", method="mle")$estimate
	#meanCorrF<-mean(subset(gspm,corrFactor<2)$corrFactor)
	meanCorrF<-gammdist[1]/gammdist[2]
	medCorrF<-qgamma(0.5,gammdist[1],gammdist[2])
	gdupper<-qgamma(0.975,gammdist[1],gammdist[2])
	gdlower<-qgamma(0.025,gammdist[1],gammdist[2])
	propUpper<-gdupper/meanCorrF;propLower<-gdlower/meanCorrF
	
	## Estimate number of seals per map
	numSealstmp<-aggregate(numTags~mapId+taggerId,tgvutm,sum)
	numSeals<-aggregate(numTags~mapId,numSealstmp,FUN=function(x){y<-ceiling(mean(x));return(y)})
	numSeals$estNumSeals<-round(numSeals$numTags*meanCorrF)
	numSeals$uclNumSeals<-round(numSeals$estNumSeals*propUpper)
	numSeals$lclNumSeals<-round(numSeals$estNumSeals*propLower)
	
	numSeals<-merge(numSeals,maps[,c("mapId","overlayId")],by="mapId",all.x=T)
	numSeals<-merge(numSeals,overlays[,c("overlayId","acquisition_date","satId","external_reference")],by="overlayId",all.x=T)
	numSeals$catalogId<-ifelse(grepl("amazonaws",numSeals$external_reference),substr(numSeals$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",numSeals$external_reference),substr(numSeals$external_reference,38,53),as.character(numSeals$external_reference)))
	#numSeals<-subset(numSeals,!is.na(acquisition_date))	#we lose 435 maps that do not have acDate, from 4500 maps we down to 4065
	numSeals$acDate<-format(numSeals$acquisition_date,"%Y%m%d")
	numSeals$acYear<-as.integer(format(numSeals$acquisition_date,"%Y"))
	numSeals$acHour<-as.integer(format(numSeals$acquisition_date,"%H"))
	numSeals$acquisition_date<-as.character(numSeals$acquisition_date)
	
	numCol<-merge(numSeals,colinfo[,c("mapId","Colony")],by="mapId")	#Only 786 with info!!
	colNmaps<-aggregate(estNumSeals~acquisition_date+Colony+catalogId,data=numCol,NROW)
	names(colNmaps)<-gsub("estNumSeals","numMaps",names(colNmaps))
	colest<-aggregate(estNumSeals~acquisition_date+Colony+catalogId,data=numCol,sum)
	colfeat<-aggregate(numTags~acquisition_date+Colony+catalogId,data=numCol,sum)
	colupper<-aggregate(uclNumSeals~acquisition_date+Colony+catalogId,data=numCol,sum)
	collower<-aggregate(lclNumSeals~acquisition_date+Colony+catalogId,data=numCol,sum)
	
	colEstimates<-merge(colNmaps,colest,by=c("acquisition_date","Colony","catalogId"))
	colEstimates<-merge(colEstimates,colfeat,by=c("acquisition_date","Colony","catalogId"))
	colEstimates<-merge(colEstimates,colupper,by=c("acquisition_date","Colony","catalogId"))
	colEstimates<-merge(colEstimates,collower,by=c("acquisition_date","Colony","catalogId"))
	
	colsurveyed<-merge(colinfo,maps,by="mapId",all.x=T)
	colsurveyed<-merge(colsurveyed,overlays,by="overlayId")
	colsurveyed$acquisition_date<-as.character(colsurveyed$acquisition_date)
	colsurveyed$catalogId<-ifelse(grepl("amazonaws",colsurveyed$external_reference),substr(colsurveyed$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",colsurveyed$external_reference),substr(colsurveyed$external_reference,38,53),as.character(colsurveyed$external_reference)))
	colsurveyed$mapsSurveyed<-1
	numInColSurveyed<-aggregate(mapsSurveyed~Colony+acquisition_date+catalogId,data=colsurveyed,sum)
	
	colEstimates<-merge(colEstimates,numInColSurveyed,by=c("acquisition_date","Colony","catalogId"),all.x=T)
	colEstimates$corrMethod<-corrMethod
	reslst<-list(numSeals=numSeals,colEstimates=colEstimates,gammdist=gammdist)
	return(reslst)
}

## The following two functions estimate the number of seals in each map inspected by ML
## One does it based on features, and the other based on tags
# gspm is the table of estimates of the correction factor for all taggers at a particular threshold
# crtm is the table of features per map, with mapId attribution
# tgvutm is the table of tags, with mapId attribution
# maps and overlays are the original data tables listing the maps and overlays
# mlCounts is the vector of mapIds of maps inspected by ML, and the number of seals she found in each
# corrMethod is a string naming the method used to calculate the correction factor
getMLEstimates_byFeature<-function(gspm,crtm,maps,overlays,mlCounts,corrMethod){
	cutrows<-nrow(gspm)-ceiling(nrow(gspm)*0.05)
	gspm<-gspm[order(gspm$corrFactor),]
	gspm<-gspm[1:cutrows,]
	gammdist<-fitdist(gspm$corrFactor, "gamma", method="mle")$estimate
	#meanCorrF<-mean(subset(gspm,corrFactor<2)$corrFactor)
	meanCorrF<-qgamma(0.5,gammdist[1],gammdist[2])
	gdupper<-qgamma(0.975,gammdist[1],gammdist[2])
	gdlower<-qgamma(0.025,gammdist[1],gammdist[2])
	propUpper<-gdupper/meanCorrF;propLower<-gdlower/meanCorrF
	
	## Estimate number of seals per map
	numSeals<-aggregate(numFeatures~mapId,crtm,sum)
	numSeals$estNumSeals<-round(numSeals$numFeatures*meanCorrF)
	numSeals$uclNumSeals<-round(numSeals$estNumSeals*propUpper)
	numSeals$lclNumSeals<-round(numSeals$estNumSeals*propLower)
	
	numSeals<-merge(numSeals,maps[,c("mapId","overlayId")],by="mapId",all.x=T)
	numSeals<-merge(numSeals,overlays[,c("overlayId","acquisition_date","satId","external_reference")],by="overlayId",all.x=T)
	numSeals$catalogId<-ifelse(grepl("amazonaws",numSeals$external_reference),substr(numSeals$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",numSeals$external_reference),substr(numSeals$external_reference,38,53),as.character(numSeals$external_reference)))
	#numSeals<-subset(numSeals,!is.na(acquisition_date))	#we lose 435 maps that do not have acDate, from 4500 maps we down to 4065
	numSeals$acDate<-format(numSeals$acquisition_date,"%Y%m%d")
	numSeals$acYear<-as.integer(format(numSeals$acquisition_date,"%Y"))
	numSeals$acHour<-as.integer(format(numSeals$acquisition_date,"%H"))
	numSeals$acquisition_date<-as.character(numSeals$acquisition_date)
	
	numML<-merge(numSeals,mlCounts,by="mapId")	
	numML$corrMethod<-corrMethod
	
	return(numML)
}

getMLEstimates_byTag<-function(gspm,tgvutm,maps,overlays,mlCounts,corrMethod){
	cutrows<-nrow(gspm)-ceiling(nrow(gspm)*0.05)
	gspm<-gspm[order(gspm$corrFactor),]
	gspm<-gspm[1:cutrows,]
	gammdist<-fitdist(gspm$corrFactor, "gamma", method="mle")$estimate
	#meanCorrF<-mean(subset(gspm,corrFactor<2)$corrFactor)
	meanCorrF<-qgamma(0.5,gammdist[1],gammdist[2])
	gdupper<-qgamma(0.975,gammdist[1],gammdist[2])
	gdlower<-qgamma(0.025,gammdist[1],gammdist[2])
	propUpper<-gdupper/meanCorrF;propLower<-gdlower/meanCorrF
	
	## Estimate number of seals per map
	numSeals<-aggregate(numTags~mapId,tgvutm,sum)
	numSeals$estNumSeals<-round(numSeals$numTags*meanCorrF)
	numSeals$uclNumSeals<-round(numSeals$estNumSeals*propUpper)
	numSeals$lclNumSeals<-round(numSeals$estNumSeals*propLower)
	
	numSeals<-merge(numSeals,maps[,c("mapId","overlayId")],by="mapId",all.x=T)
	numSeals<-merge(numSeals,overlays[,c("overlayId","acquisition_date","satId","external_reference")],by="overlayId",all.x=T)
	numSeals$catalogId<-ifelse(grepl("amazonaws",numSeals$external_reference),substr(numSeals$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",numSeals$external_reference),substr(numSeals$external_reference,38,53),as.character(numSeals$external_reference)))
	#numSeals<-subset(numSeals,!is.na(acquisition_date))	#we lose 435 maps that do not have acDate, from 4500 maps we down to 4065
	numSeals$acDate<-format(numSeals$acquisition_date,"%Y%m%d")
	numSeals$acYear<-as.integer(format(numSeals$acquisition_date,"%Y"))
	numSeals$acHour<-as.integer(format(numSeals$acquisition_date,"%H"))
	numSeals$acquisition_date<-as.character(numSeals$acquisition_date)
	
	numML<-merge(numSeals,mlCounts[,c("mapId","mlcount")],by="mapId")
	numML$corrMethod<-corrMethod
	
	return(numML)
}

## The following function returns the number of tags in maps inspected by the tagger
# tt is the taggerId for the tagger
# tagsSel is the subset of tags inspected by taggers who share maps with Michelle
# viewsSel is the set of map views by taggers who share maps with ML
# taggerMaps is the collection of maps inspected by the tagger
getTaggerNumTags<-function(tt,tagsSel,viewsSel,taggerMaps){
	qq<-subset(tagsSel,taggerId==tt)
	qq<-merge(qq,viewsSel,by="mapViewId",all.x=TRUE)
	ntgs<-nrow(subset(qq,mapId %in% taggerMaps))
	return(ntgs)
}

## The following function returns the number of maps with tags given the crowdrank threshold
## the number of tags in those maps, and the average number of tags per map
# tags is the table with all tags and the views where these were placed
# views is the table of all map views
# maps is the table with information about each map and its link to an image (overlay)
# overlays is the table with data on the satellite images 
# colinfo is the table that links a map with a colony name
getTagStats<-function(tags,views,maps,overlays,colinfo){
	#need to know all maps with tags, and the average number of tags per map
	tgv<-merge(tags,views[,c("mapViewId","mapId")],by="mapViewId",all.x=TRUE)
	tgvc<-merge(tgv,colinfo[,c("mapId","Colony")],by="mapId",all.x=TRUE)
	tgvc<-subset(tgvc,!is.na(Colony))
	tgmp<-merge(tgvc,maps[,c("mapId","overlayId")],by="mapId",all.x=TRUE)
	tgco<-merge(tgmp,overlays[,c("overlayId","acquisition_date")],by="overlayId",all.x=TRUE)
	tgco$satdate<-format(tgco$acquisition_date,"%Y-%m-%d %H:%M:%S")
	tagsPerColonyDate<-aggregate(mapId~Colony+satdate,data=tgco,NROW)
	names(tagsPerColonyDate)<-c("Colony","ImgDate","NumTags")
	mapsPerColonyDate<-aggregate(mapId~Colony+satdate,data=tgco,FUN=function(x){y<-NROW(unique(x));return(y)})
	names(mapsPerColonyDate)<-c("Colony","ImgDate","NumMapsWtags")
	res<-merge(tagsPerColonyDate,mapsPerColonyDate,by=c("Colony","ImgDate"))
	return(res)
}


## The below is the function that calculates the correction factor defining F as the number of features in maps inspected by the specific tagger
# rankdata is the list of taggers and their crowdrank scores
# crthr is the desired crowdrank filtering threshold
# taggersSel is the list of taggers who shared maps with Michelle
# tagsSel is the subset of tags inspected by taggers who share maps with Michelle
# tagsML is the set of Michelle's tags
# crtm is the set of features with example tag, map, and view, and longlat
# viewsSel is the set of map views by taggers who share maps with ML
# ProbS is the global probability of seals in an image
# corrProbF is the value of a probability used to correct the estimated ProbF
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps
getTaggerProbabilities_tagsOnlybyTag<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,crtm,viewsSel,ProbS,corrProbF=1){
	ttdf<-subset(rankdata,taggerScore>=crthr);	#only taggers matching or besting the threshold
	taggers<-subset(taggersSel,taggersSel %in% ttdf$taggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
	
	#calculate the probabilities for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=taggers,.fun=function(tt,viewsSel,crtm,tagsSel,tagsML){
				#calculating ProbF...(number of maps with tags among the maps inspected by this tagger)
				#need the number of maps inspected by this tagger
				taggerMaps<-unique(subset(viewsSel,taggerId==tt)$mapId)
				tagger_nMaps<-NROW(taggerMaps)
				#determine how many of these have tags
				taggerTags<-subset(tagsSel,taggerId==tt);totalTags<-nrow(taggerTags)
				taggerTags<-merge(taggerTags,viewsSel[,c("mapViewId","mapId")],by="mapViewId",all.x=T)
				tagger_nTagMaps<-NROW(unique(taggerTags$mapId))
				mlTags<-subset(tagsML,mapId %in% taggerMaps)	#if this is 0, it's because ML did not find a seal in any of the maps this tagger inspected.
				totalSeals<-nrow(mlTags)
				ProbF<-totalTags/tagger_nMaps
				
				if(totalTags>0){	#features were found in maps inspected by the tagger
					#calculating ProbFS...
					#need the list of maps from this tagger also shared with ML - this is taggerMaps
					#need the list of tags in these maps
					
					if(nrow(mlTags)==0){
						ttdf<-data.frame(taggerId=tt,ss=0,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=0,ProbF=0)
					}else{
						ss<-getTaggerSS(taggerTags,mlTags)
						ProbFS<-ss/totalSeals
						
						ttdf<-data.frame(taggerId=tt,ss=ss,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=ProbFS,ProbF=ProbF)
					}
				}else{	#no features found in maps inspected by this tagger
					ttdf<-data.frame(taggerId=tt,ss=0,nTags=totalTags,nSeals=0,nMaps=tagger_nMaps,ProbFS=0,ProbF=0,ProbF=ProbF)
				}
				
				
				return(ttdf)},viewsSel=viewsSel,crtm=crtm,tagsSel=tagsSel,tagsML=tagsML)
	
	taggerProbs<-subset(taggerProbs,ProbFS>0)
	if(nrow(taggerProbs)>0){
		taggerProbs$ProbS<-ProbS
		taggerProbs$ProbSF<-taggerProbs$ProbS*taggerProbs$ProbFS/taggerProbs$ProbF	
		taggerProbs$corrFactor<-taggerProbs$ProbSF/taggerProbs$ProbFS
	}else{ #none of the computed probabilities permits the estimation of corrFactor
		taggerProbs<-NA
	}
	
	return(taggerProbs)
}


################################### DATA
## need colony info - ML had data from 7755 maps, but we have only 4065 with data
colinfo<-read.csv("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/ColonyData/MapsInColonies/mapsInColonies.csv")

mid<-21758509	## This is Michelle LaRue's tagger Id

## load the data, prepare for analyses
load("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/dataFromGeoJson3.RData")
tags<-unique(tags)
tgutm<-convertToUTM(tags,lonfield="tagcoords.x1",latfield="tagcoords.x2")

crt<-merge(crowd[,c("tagId","score","agremnt","sensor")],tgutm[,c("tagId","mapViewId","easting","northing")],by="tagId",all.x=T)
crt<-subset(crt,!is.na(mapViewId))
crtm<-merge(crt,views[,c("mapViewId","mapId")],by="mapViewId",all.x=T)
crtm<-subset(crtm,!is.na(mapId))
crtm$numFeatures<-1

tgvutm<-merge(tgutm,views[,c("mapViewId","mapId")],by="mapViewId",all.x=T)
tgvutm$numTags<-1
mltdf<-subset(tgvutm,taggerId==mid)
mlCounts<-aggregate(numTags~mapId,data=mltdf,FUN=sum)
names(mlCounts)<-c("mapId","mlcount")

## We want to know how many times a map was vewed:
aggviews<-aggregate(taggerId~mapId,data=views,NROW);names(aggviews)<-c("mapId","numViews")

## So...
#GenProbF<-NROW(unique(crtm$mapId))/sum(aggviews>0)
GenProbF<-nrow(crtm)/sum(aggviews>0)
#Calculating the overall Prob[S] and expectation for CorrFactor
mlMaps<-unique(subset(views,taggerId==mid)$mapId)
#ProbS<-NROW(unique(subset(crtm,mapId %in% mlMaps)$mapId))/NROW(mlMaps)
ProbS<-nrow(subset(crtm,mapId %in% mlMaps))/NROW(mlMaps)
expectCorrF<-ProbS/GenProbF
print(paste("Expectation for correction factor:",round(expectCorrF,2)))

## If ProbFS is ss/(ss+ns) and ProbSF is ProbFS * expectCorrF
## We would calculate a mean ProbFS, and a mean ProbSF, 
## But note that expectCorrF is a constant that we can factor out of the averages, so
## CorrFactor=mProbSF/mProbFS is equivalent to CorrFactor=ProbSF/ProbFS
## Since  ProbSF = ProbFS*expCorrF we can see that ProbFS cancels out and indeed expectation(CorrFactor)=expCorrF


############ Now calculate by map: 
## ProbFS: It is still ss/totalSeals, where for each tagger, where ss is the 11 case and totalSeals is the seals ML found in all maps inspected by both she and the tagger
## ProbF: In reality ProbF is unique to each tagger; it has several interpretations (see below) 
## ProbS: This is constant - the probability that ML found a seal
## With the above, then calculate ProbSF, and corrFactor

scores2<-read.csv("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/antarctica_exports20181102/antarctica_pilot_2_user_scores_20181102.csv",stringsAsFactors=F)
scores1<-read.csv("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/antarctica_exports20181102/antarctica_pilot_1_user_scores_20181102.csv",stringsAsFactors=F)
sc1id<-unique(scores1$tagger_id); sc2id<-unique(scores2$tagger_id)
addsc<-subset(sc1id,!sc1id %in% sc2id)
scores1<-subset(scores1,tagger_id %in% addsc)
rankdata<-rbind(scores1,scores2)
rankdata<-rankdata[,c("tagger_id","tagger_score")];names(rankdata)<-c("taggerId","taggerScore")	#tagger rank scores

taggersSel<-unique(subset(views,mapId %in% mlMaps,select="taggerId")$taggerId)
taggersSel<-subset(taggersSel,taggersSel!=mid)		#taggers who share maps with ML

## What tags to use? For ProbFS...
tagsSel<-subset(tgutm,taggerId %in% taggersSel)
tagsML<-subset(tgutm,taggerId==mid)
tagsML<-merge(tagsML,views[,c("mapId","mapViewId")],by="mapViewId",all.x=T)

## What maps can we use? For ProbF
viewsSel<-subset(views,taggerId %in% c(taggersSel,mid))

## Need the probability of a feature being found on a map and the probability of a tag being placed on a map
totalMapsInspected<-NROW(unique(views$mapId))
totalMapsWtags<-NROW(unique(tgvutm$mapId))
totalMapsWfeat<-NROW(unique(crtm$mapId))
probFeatInMap<-totalMapsWfeat/totalMapsInspected
probTagInMap<-totalMapsWtags/totalMapsInspected
probTagAsFeat<-totalMapsWfeat/totalMapsWtags

compML_feat<-list();compML_tag<-list()

##############
## ProbF = number of  tags in maps inspected/total number of maps inspected
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps 
gspm<-getTaggerProbabilities_tagsOnlybyTag(rankdata=rankdata,crthr=0.8,taggersSel=taggersSel,tagsSel=tagsSel,tagsML=tagsML,crtm=crtm,viewsSel=viewsSel,ProbS=ProbS)
res_tagsOnlybyTag<-getColonyEstimates_byTag(gspm=gspm,tgvutm=tgvutm,maps=maps,overlays=overlays,colinfo=colinfo,corrMethod="tagsOnlybyTag")
#compML_feat$trueFeatInMap<-getMLEstimates_byFeature(gspm=gspm,crtm=crtm,maps=maps,overlays=overlays,mlCounts,corrMethod="tagsOnlybyTag")
compML_tag$tagsOnlybyTag<-getMLEstimates_byTag(gspm=gspm,tgvutm=tgvutm,maps=maps,overlays=overlays,mlCounts,corrMethod="tagsOnlybyTag")

#need to know all maps with tags, and the average number of tags per map
tagStats<-getTagStats(tags=tags,views=views,maps=maps,overlays=overlays,colinfo=colinfo)

## Save all results in a single data file
save(res_tagsOnlybyTag,
		tagStats,compML_feat,compML_tag,
		file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/colonyEstimates.RData")


