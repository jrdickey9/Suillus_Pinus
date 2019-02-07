#also note that both R and Java need to be the same bit (either 32 or 64) to be compatible to run
options(java.parameters = c("‐Xss2560k", "‐Xmx8g", "‐Xms2g"))
## increase memory size of the JVM, this code may prevent memory issues of Maxent.jar

install.packages("usdm")
install.packages("dismo")
install.packages("raster")
install.packages("rgeos")
install.packages("rJava")
install.packages("devtools")
library(devtools)
install_github("bobmuscarella/ENMeval")
install.packages("parallel")
install_github("danlwarren/ENMTools")
library(raster)
library(dismo)
library(rJava)
library(usdm)
require(rgeos)
require(rJava)
library(sp)
library(rgdal)

library(ENMeval)
library(ENMTools)

setwd("/Volumes/S&V_LineaPro_Apr2012")
mycrs<-CRS("+init=epsg:4326 +proj=lcc +lat_1=49 +lat_2=77 +lon_0=-95 +datum=WGS84") #used to define
newcrs <- c("+proj=longlat +datum=WGS84") #used to transform

#RCP4.5 INMCM4 2055
sero<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/INM-CM4_rcp45_2055", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP4.5INMCM42055<-stack(sero) #creating raster brick
crs(RCP4.5INMCM42055)<-mycrs #defining
e<-extent(-875000,3010000,1800000,10418000) #creating things to crop
RCP4.5INMCM42055<-crop(RCP4.5INMCM42055,e)
collin.test<-vifstep(RCP4.5INMCM42055, th=0.7) #toss out variables with multicollinearity @ threshold 0.7 #identify collinear variables that should be excluded at
vif.sero<-vif(RCP4.5INMCM42055)
s0<-vifcor(RCP4.5INMCM42055,th=0.7) #corelates
RCP4.5INMCM42055<-exclude(RCP4.5INMCM42055,s0) #throwing out collinearity.
RCP4.5INMCM42055<-projectRaster(from=RCP4.5INMCM42055, crs = newcrs) #new projection
RCP4.5INMCM42055<-RCP4.5INMCM42055
writeRaster(RCP4.5INMCM42055, filename = "RCP4.5INMCM42055", overwrite=TRUE)

#RCP8.5 INMCM4 2055
uno<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/INM-CM4_rcp85_2055", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP8.5INMCM42055<-stack(uno) #creating rasterbrick
crs(RCP8.5INMCM42055)<-mycrs #defining
e<-extent(-875000,3010000,1800000,10418000) #creating things to crop
RCP8.5INMCM42055<-crop(RCP8.5INMCM42055,e)
collin.test<-vifstep(RCP8.5INMCM42055, th=0.7)
vif.uno<-vif(RCP8.5INMCM42055)
s1<-vifcor(RCP8.5INMCM42055,th=0.7)
RCP8.5INMCM42055<-exclude(RCP8.5INMCM42055,s1)
RCP8.5INMCM42055<-projectRaster(RCP8.5INMCM42055, crs = newcrs) #project
writeRaster(RCP8.5INMCM42055, filename = "RCP8.5INMCM42055", overwrite=TRUE)

#RCP4.5 INMCM4 2085
dos<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/INM-CM4_rcp45_2085", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP4.5INMCM42085<-stack(dos)
crs(RCP4.5INMCM42085) <- mycrs #defining
e<-extent(-875000,3010000,1800000,10418000)
RCP4.5INMCM42085<-crop(RCP4.5INMCM42085,e)
collin.test<-vifstep(RCP4.5INMCM42085, th=0.7)
vif.dos<-vif(RCP4.5INMCM42085)
s2<-vifcor(RCP4.5INMCM42085,th=0.7)
RCP4.5INMCM42085<-exclude(RCP4.5INMCM42085,s2)
RCP4.5INMCM42085<- projectRaster(RCP4.5INMCM42085, crs = newcrs) #new projection
writeRaster(RCP4.5INMCM42085, filename = "RCP4.5INMCM42085", overwrite=TRUE)


#RCP8.5 INMCM4 2085
tres<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/INM-CM4_rcp85_2085", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP8.5INMCM42085<-stack(tres)
crs(RCP8.5INMCM42085)<-mycrs
e<-extent(-875000,3010000,1800000,10418000)
RCP8.5INMCM42085<-crop(RCP8.5INMCM42085,e)
collin.test<-vifstep(RCP8.5INMCM42085, th=0.7)
vif.tres<-vif(RCP8.5INMCM42085)
s3<-vifcor(RCP8.5INMCM42085,th=0.7)
RCP8.5INMCM42085<-exclude(RCP8.5INMCM42085,s3)
RCP8.5INMCM42085<-projectRaster(RCP8.5INMCM42085, crs = newcrs)
writeRaster(RCP8.5INMCM42085, filename = "RCP8.5INMCM42085", overwrite=TRUE)


#RCP4.5 CCSM4 2055
cuatro<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/CCSM4_rcp45_2055", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP4.5CCSM42055<-stack(cuatro)
crs(RCP4.5CCSM42055)<-mycrs
e<-extent(-875000,3010000,1800000,10418000)
RCP4.5CCSM42055<-crop(RCP4.5CCSM42055,e)
collin.test<-vifstep(RCP4.5CCSM42055, th=0.7)
vif.cuatro<-vif(RCP4.5CCSM42055)
s4<-vifcor(RCP4.5CCSM42055,th=0.7)
RCP4.5CCSM42055<-exclude(RCP4.5CCSM42055,s4)
RCP4.5CCSM42055<-projectRaster(RCP4.5CCSM42055, crs = newcrs)
writeRaster(RCP4.5CCSM42055, filename = "RCP4.5CCSM42055", overwrite=TRUE)


#RCP8.5 CCSM4 2055
cinco<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/CCSM4_rcp85_2055", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP8.5CCSM42055<-stack(cinco) #brick
crs(RCP8.5CCSM42055)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP8.5CCSM42055<-crop(RCP8.5CCSM42055,e)
collin.test<-vifstep(RCP8.5CCSM42055, th=0.7)
vif.cinco<-vif(RCP8.5CCSM42055)
s5<-vifcor(RCP8.5CCSM42055,th=0.7)
RCP8.5CCSM42055<-exclude(RCP8.5CCSM42055,s5)
RCP8.5CCSM42055<-projectRaster(RCP8.5CCSM42055, crs = newcrs) #PROJECT
writeRaster(RCP8.5CCSM42055, filename = "RCP8.5CCSM42055", overwrite=TRUE)


#RCP4.5 CCSM4 2085
seis<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/CCSM4_rcp45_2085", recursive=TRUE, full.names=TRUE,pattern="\\.asc$")
RCP4.5CCSM42085<-stack(seis)
crs(RCP4.5CCSM42085)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP4.5CCSM42085<-crop(RCP4.5CCSM42085,e)
collin.test<-vifstep(RCP4.5CCSM42085, th=0.7)
vif.seis<-vif(RCP4.5CCSM42085)
s6<-vifcor(RCP4.5CCSM42085,th=0.7)
RCP4.5CCSM42085<-exclude(RCP4.5CCSM42085,s6)
RCP4.5CCSM42085<-projectRaster(RCP4.5CCSM42085, crs = newcrs) #PROJECT
writeRaster(RCP4.5CCSM42085, filename = "RCP4.5CCSM42085", overwrite=TRUE)


#RCP8.5 CCSM4 2085
siete<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/CCSM4_rcp85_2085", recursive=TRUE, full.names=TRUE, pattern="\\.asc$")
RCP8.5CCSM42085<-stack(siete)
crs(RCP8.5CCSM42085)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP8.5CCSM42085<-crop(RCP8.5CCSM42085,e)
collin.test<-vifstep(RCP8.5CCSM42085, th=0.7)
vif.siete<-vif(RCP8.5CCSM42085)
s7<-vifcor(RCP8.5CCSM42085,th=0.7)
RCP8.5CCSM42085<-exclude(RCP8.5CCSM42085,s7)
RCP8.5CCSM42085<-projectRaster(RCP8.5CCSM42085, crs = newcrs) #PROJECT
writeRaster(RCP8.5CCSM42085, filename = "RCP8.5CCSM42085", overwrite=TRUE)


#RCP4.5 GFDLCM3 2055
ocho<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/GFDL-CM3_rcp45_2055", recursive=TRUE, full.names=TRUE, pattern="\\.asc$")
RCP4.5GFDLCM32055<-stack(ocho)
crs(RCP4.5GFDLCM32055)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP4.5GFDLCM32055<-crop(RCP4.5GFDLCM32055,e)
collin.test<-vifstep(RCP4.5GFDLCM32055, th=0.7)
vif.ocho<-vif(RCP4.5GFDLCM32055)
s8<-vifcor(RRCP4.5GFDLCM32055,th=0.7)
RCP4.5GFDLCM32055<-exclude(RCP4.5GFDLCM32055,s8)
RCP4.5GFDLCM32055<-projectRaster(RCP4.5GFDLCM32055, crs = newcrs) #PROJECT
writeRaster(RCP4.5GFDLCM32055, filename = "RCP4.5GFDLCM32055", overwrite=TRUE)


#RCP8.5 GFDLCM3 2055
nueve<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/GFDL-CM3_rcp85_2055", recursive=TRUE, full.names=TRUE, pattern="\\.asc$")
RCP8.5GFDLCM32055<-stack(nueve)
crs(RCP8.5GFDLCM32055)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP8.5GFDLCM32055<-crop(RCP8.5GFDLCM32055,e) #tossing out correlated layers
collin.test<-vifstep(RCP8.5GFDLCM32055, th=0.7)
vif.nueve<-vif(RCP8.5GFDLCM32055)
s9<-vifcor(RCP8.5GFDLCM32055,th=0.7)
RCP8.5GFDLCM32055<-exclude(RCP8.5GFDLCM32055,s9)
RCP8.5GFDLCM32055<-projectRaster(RCP8.5GFDLCM32055, crs = newcrs) #PROJECT
writeRaster(RCP8.5GFDLCM32055, filename = "RCP8.5GFDLCM32055", overwrite=TRUE)


#RCP4.5 GFDLCM3 2085
diez<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/GFDL-CM3_rcp45_2085", recursive=TRUE, full.names=TRUE, pattern="\\.asc$")
RCP4.5GFDLCM32085<-stack(diez)
crs(RCP4.5GFDLCM32085)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP4.5GFDLCM32085<-crop(RCP4.5GFDLCM32085,e) #toss'em out
collin.test<-vifstep(RCP4.5GFDLCM32085, th=0.7)
vif.diez<-vif(RCP4.5GFDLCM32085)
s10<-vifcor(RCP4.5GFDLCM32085,th=0.7)
RCP4.5GFDLCM32085<-exclude(RCP4.5GFDLCM32085,s10)
RCP4.5GFDLCM32085<-projectRaster(RCP4.5GFDLCM32085, crs = newcrs) #PROJECT
writeRaster(RCP4.5GFDLCM32085, filename = "RCP4.5GFDLCM32085", overwrite=TRUE)


#RCP8.5 GFDLCM3 2085
once<-dir("/Volumes/S&V_LineaPro_Apr2012/future climate ASCII/GFDL-CM3_rcp85_2085", recursive=TRUE, full.names=TRUE, pattern="\\.asc$")
RCP8.5GFDLCM32085<-stack(once)
crs(RCP8.5GFDLCM32085)<-mycrs #DEFINE
e<-extent(-875000,3010000,1800000,10418000) #EXTENT
RCP8.5GFDLCM32085<-crop(RCP8.5GFDLCM32085,e)
#collin.test<-vifstep(RCP8.5GFDLCM32085, th=0.7)
vif.once<-vif(RCP8.5GFDLCM32085)
s11<-vifcor(RCP8.5GFDLCM32085,th=0.7)
RCP8.5GFDLCM32085<-exclude(RCP8.5GFDLCM32085,s11)
RCP8.5GFDLCM32085<-projectRaster(RCP8.5GFDLCM32085, crs = newcrs) #PROJECT
writeRaster(RCP8.5GFDLCM32085, filename = "RCP8.5GFDLCM32085")


#CURRENT SCENARIO
doce<-dir("/Volumes/S&V_LineaPro_Apr2012/CurrentNormal", recursive=TRUE, full.names=TRUE, pattern="\\.asc$")#load
CurrentNormal<-stack(doce) #brick
crs(CurrentNormal)<-mycrs #define
e<-extent(-875000,3010000,1800000,10418000) #extent
CurrentNormal<-crop(CurrentNormal,e) #crop
collin.test<-vifstep(CurrentNormal, th=0.7)
vif.once<-vif(CurrentNormal)
s11<-vifcor(CurrentNormal,th=0.7)
CurrentNormal<-exclude(CurrentNormal,s11) #collinearity
CurrentNormal<-projectRaster(CurrentNormal, crs = newcrs) #project
writeRaster(CurrentNormal, filename = "CurrentNormal") #write

save(list=ls(), file="the_future.rda")
setwd("/Volumes/S&V_LineaPro_Apr2012")
load(file="the_future.rda")

#####################################occurence data for Pinus strobus
Pinus.occ<-read.csv("newPinus.csv")
head(Pinus.occ)
names(Pinus.occ)
summary(Pinus.occ)
str(Pinus.occ)
Pinus.occ<-Pinus.occ[,3:4]
Pinus.occ<-Pinus.occ[!is.na(rowSums(extract(RCP4.5INMCM42055,Pinus.occ))),]
coordinates(object=Pinus.occ)<-~ Lon+Lat #make spatial
plot(RCP4.5INMCM42055[[1]]) #visualizing data
points(Pinus.occ) #points!
save(list=ls(), file="the_future.rda")

load(file="transformp3.rda")#for the next day, need to reload in each file for maxent to work off a perm file not a temp one. #need to do this when R restarts to reload raster files into global environment.
RCP4.5INMCM42055<-brick("RCP4.5INMCM42055.gri") #0
RCP8.5INMCM42055<-brick("RCP8.5INMCM42055.gri") #1
RCP4.5INMCM42085<-brick("RCP4.5INMCM42085.gri") #2
RCP8.5INMCM42085<-brick("RCP8.5INMCM42085.gri") #3
RCP4.5CCSM42055<-brick("RCP4.5CCSM42055.gri") #4
RCP8.5CCSM42055<-brick("RCP8.5CCSM42055.gri") #5
RCP4.5CCSM42085<-brick("RCP4.5CCSM42085.gri") #6
RCP8.5CCSM42085<-brick("RCP8.5CCSM42085.gri") #7
RCP8.5GFDLCM32055<-brick("RCP8.5GFDLCM32055.gri") #8
RCP4.5GFDLCM32055<-brick("RCP4.5GFDLCM32055.gri") #9
RCP4.5GFDLCM32085<-brick("RCP4.5GFDLCM32085.gri") #10
RCP8.5GFDLCM32085<-brick("RCP8.5GFDLCM32085.gri") #11 #these are all cropped, correlated, and converted.
CurrentNormal<-brick("CurrentNormal.gri") #reloading in
load("the_future.rda")

# download maxent.jar 3.3.3k, and place the file in the desired folder
utils::download.file(url="https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar",destfile=paste0(system.file("java", package="dismo"),"/maxent.jar"), mode="wb")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='') #most current MaxEnt 3.4.1
# also note that both R and Java need to be the same bit (either 32 or 64) to be compatible to run
options( java.parameters = c("‐Xss2560k", "‐Xmx2g") ) ## increase memory size of the JVM, this code may prevent memory issues of Maxent.jar
.jinit(classpath = NULL, parameters = getOption("java.parameters"), silent = FALSE, force.init = FALSE)
########################################future Pinus model projections
Pinus.mod.sero<-dismo::maxent(x=RCP4.5INMCM42055,p=Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_4.5INMCM42055"),args=c('responsecurves=true'))
Pinus.mod.uno<-dismo::maxent(RCP8.5INMCM42055,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_8.5INMCM42055"),args=c('responsecurves=true'))
Pinus.mod.dos<-dismo::maxent(RCP4.5INMCM42085,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_4.5INMCM42085"),args=c('responsecurves=true'))
Pinus.mod.tres<-dismo::maxent(RCP8.5INMCM42085,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_8.5INMCM42085"),args=c('responsecurves=true'))
Pinus.mod.quatro<-dismo::maxent(RCP4.5CCSM42055,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_4.5CCSM42055"),args=c('responsecurves=true'))
Pinus.mod.cinco<-dismo::maxent(RCP8.5CCSM42055,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_8.5CCSM42055"),args=c('responsecurves=true'))
Pinus.mod.seis<-dismo::maxent(RCP4.5CCSM42085,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_4.5CCSM42085"),args=c('responsecurves=true'))
Pinus.mod.siete<-dismo::maxent(RCP8.5CCSM42085,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_8.5CCSM42085"),args=c('responsecurves=true'))
Pinus.mod.ocho<-dismo::maxent(RCP4.5GFDLCM32055,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_4.5GFDLCM32055"),args=c('responsecurves=true'))
Pinus.mod.nueve<-dismo::maxent(RCP8.5GFDLCM32055,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_8.5GFDLCM32055"),args=c('responsecurves=true'))
Pinus.mod.diez<-dismo::maxent(RCP4.5GFDLCM32085,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_4.5GFDLCM32085"),args=c('responsecurves=true'))
Pinus.mod.once<-dismo::maxent(RCP8.5GFDLCM32085,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/P_8.5GFDLCM32085"),args=c('responsecurves=true')) #already saved
Pinus.mod.current<-dismo::maxent(CurrentNormal,Pinus.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/CurrentPinus"), args=c('responsecurves=true')) 

save(Pinus.mod.sero,Pinus.mod.uno,Pinus.mod.dos,Pinus.mod.tres,Pinus.mod.quatro,Pinus.mod.cinco,Pinus.mod.seis,Pinus.mod.siete,Pinus.mod.ocho,Pinus.mod.nueve,Pinus.mod.diez,Pinus.mod.once,Pinus.mod.current, file="Pinus.mods.rda")


##############################occurence data for Suillus 
suillydat<-read.csv(file="Suilly_occ.csv")
set.seed(150) #set this 
trueSuilly.occ<-suillydat[sample(nrow(suillydat), 54), ] #check to make sure n matches Pinus occ
trueSuilly.occ<-trueSuilly.occ[,1:2] #selecting out long lats, make sure data is in deg.dec
S.sprag.occ<-trueSuilly.occ[!is.na(rowSums(extract(RCP4.5INMCM42055,trueSuilly.occ))),]
coordinates(trueSuilly.occ) <-~ Long + Lat #these should be headers
plot(RCP4.5INMCM42055[[1]])
points(S.sprag.occ)

##################################future Suillus model projections 
S.sprag.sero<-dismo::maxent(RCP4.5INMCM42055,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_4.5INMCM42055"),args=c('responsecurves=true'))
S.sprag.uno<-dismo::maxent(RCP8.5INMCM42055,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_8.5INMCM42055"),args=c('responsecurves=true'))
S.sprag.dos<-dismo::maxent(RCP4.5INMCM42085,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_4.5INMCM42085"),args=c('responsecurves=true'))
S.sprag.tres<-dismo::maxent(RCP8.5INMCM42085,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_8.5INMCM42085"),args=c('responsecurves=true'))
S.sprag.cuatro<-dismo::maxent(RCP4.5CCSM42055,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_4.5CCSM42055"),args=c('responsecurves=true'))
S.sprag.cinco<-dismo::maxent(RCP8.5CCSM42055,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_8.5CCSM42055"),args=c('responsecurves=true'))
S.sprag.seis<-dismo::maxent(RCP4.5CCSM42085,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_4.5CCSM42085"),args=c('responsecurves=true'))
S.sprag.siete<-dismo::maxent(RCP8.5CCSM42085,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_8.5CCSM42085"),args=c('responsecurves=true'))
S.sprag.ocho<-dismo::maxent(RCP4.5GFDLCM32055,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_4.5GFDLCM32055"),args=c('responsecurves=true'))
S.sprag.nueve<-dismo::maxent(RCP8.5GFDLCM32055,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_8.5GFDLCM32055"),args=c('responsecurves=true'))
S.sprag.diez<-dismo::maxent(RCP4.5GFDLCM32085,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_4.5GFDLCM32085"),args=c('responsecurves=true'))
S.sprag.once<-dismo::maxent(RCP8.5GFDLCM32085,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/S_8.5GFDLCM32085"),args=c('responsecurves=true'))
S.sprag.current<-dismo::maxent(CurrentNormal,S.sprag.occ,path = paste0(getwd(),"/Volumes/S&V_LineaPro_Apr2012/SuillusCurrent"),args=c('responsecurves=true'))
save(S.sprag.sero,S.sprag.uno,S.sprag.dos,S.sprag.tres,S.sprag.cuatro,S.sprag.cinco,S.sprag.seis,S.sprag.siete,S.sprag.ocho,S.sprag.nueve,S.sprag.diez,S.sprag.once, S.sprag.current, file="S.sprag.rda")

save(list=ls(), file="the_future.rda")
source(file="the_future.rda")
load("S.sprag.rda")

################################################## POST HOC ANALYSIS
#Predict the model across you study extent and visualize it
##################
Ppred0<-predict(Pinus.mod.sero,RCP4.5INMCM42055)
#plot(Ppred0) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4) #what are these decorative arguments for points()
Ppred1<-predict(Pinus.mod.uno,RCP8.5INMCM42055)
#plot(Ppred1) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred2<-predict(Pinus.mod.dos,RCP4.5INMCM42085)
#plot(Ppred2) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred3<-predict(Pinus.mod.tres,RCP8.5INMCM42085)
#plot(Ppred3) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred4<-predict(Pinus.mod.quatro, RCP4.5CCSM42055)
#plot(Ppred4) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred5<-predict(Pinus.mod.cinco, RCP8.5CCSM42055)
#plot(Ppred5) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred6<-predict(Pinus.mod.seis, RCP4.5CCSM42085)
#plot(Ppred6) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred7<-predict(Pinus.mod.siete, RCP8.5CCSM42085)
#plot(Ppred7) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred8<-predict(Pinus.mod.ocho, RCP4.5GFDLCM32055)
#plot(Ppred8) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred9<-predict(Pinus.mod.nueve, RCP8.5GFDLCM32055)
#plot(Ppred9) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred10<-predict(Pinus.mod.diez, RCP4.5GFDLCM32085)
#plot(Ppred10) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppred11<-predict(Pinus.mod.once, RCP8.5GFDLCM32085)
#plot(Ppred11) #only useful for images and visuals 
#points(Pinus.occ, pch=16, cex=.4)
Ppredcurrent<-predict(Pinus.mod.current, CurrentNormal)
#plot(Ppredcurrent)
#points(Pinus.occ, pch=16, cex=.4)

save(Ppred0,Ppred1,Ppred2,Ppred3,Ppred4,Ppred5,Ppred6,Ppred7,Ppred8,Ppred9,Ppred10,Ppred11,Ppredcurrent, file="Ppreds.rda")

px1 <- predict(S.sprag.sero, RCP4.5INMCM42055) #prediction
final_map0<- px1 # renaming for downstream processes. 
#plot(final_map0)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.uno, RCP8.5INMCM42055) #prediction
final_map1<-px1 # renaming for downstream processes. 
#plot(final_map1)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.dos, RCP4.5INMCM42085) #prediction
final_map2<-px1 #rename
#plot(final_map2)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.tres, RCP8.5INMCM42085) #prediction
final_map3<-px1 #rename
#plot(final_map3)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.cuatro, RCP4.5CCSM42055)
final_map4<-px1
#plot(final_map4)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.cinco,RCP8.5CCSM42055)
final_map5<-px1
#plot(final_map5)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.seis, RCP4.5CCSM42085)
final_map6<-px1 
#plot(final_map6)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.siete, RCP8.5CCSM42085)
final_map7<-px1
#plot(final_map7)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.ocho, RCP4.5GFDLCM32055)
final_map8<-px1
#plot(final_map8)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.nueve, RCP8.5GFDLCM32055)
final_map9<-px1
#plot(final_map9)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.diez, RCP4.5GFDLCM32085)
final_map10<-px1
#plot(final_map10)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.once, RCP8.5GFDLCM32085)
final_map11<-px1
#plot(final_map11)
#points(S.sprag.occ, pch=16, cex=.4)
px1 <- predict(S.sprag.current, CurrentNormal)
final_mapcurrent<-px1
#plot(final_mapcurrent)
#points(S.sprag.occ, pch=16, cex=.4)

save(final_map0,final_map1,final_map2,final_map3,final_map4,final_map5,final_map6,final_map7,final_map8,final_map9,final_map10,final_map11, final_mapcurrent, file="Spreds.rda")
#########################NICHE OVERLAP BETWEEN EACH MODEL PROJECTION FOR SUILLUS AND PINUS#####################
load("Ppreds.rda") # if needed
#load("Spreds.rda")

predicted.maps0<-stack(Ppred0,final_map0)
predicted.maps1<-stack(Ppred1,final_map1)
predicted.maps2<-stack(Ppred2,final_map2)
predicted.maps3<-stack(Ppred3,final_map3)
predicted.maps4<-stack(Ppred4,final_map4)
predicted.maps5<-stack(Ppred5,final_map5)
predicted.maps6<-stack(Ppred6,final_map6)
predicted.maps7<-stack(Ppred7,final_map7)
predicted.maps8<-stack(Ppred8,final_map8)
predicted.maps9<-stack(Ppred9,final_map9)
predicted.maps10<-stack(Ppred10,final_map10)
predicted.maps11<-stack(Ppred11,final_map11)
predicted.mapsCur<-stack(Ppredcurrent,final_mapcurrent)
makentargs<-make.args(RMvalues=c(1:2),fc=c("L", "LQ"),labels=TRUE) 
niche.ovr0<-calc.niche.overlap(predicted.maps0,stat="D",maxentargs)
niche.ovr1<-calc.niche.overlap(predicted.maps1,stat="D",maxentargs)
niche.ovr2<-calc.niche.overlap(predicted.maps2,stat="D",maxentargs)
niche.ovr3<-calc.niche.overlap(predicted.maps3,stat="D",maxentargs)
niche.ovr4<-calc.niche.overlap(predicted.maps4,stat="D",maxentargs)
niche.ovr5<-calc.niche.overlap(predicted.maps5,stat="D",maxentargs)
niche.ovr6<-calc.niche.overlap(predicted.maps6,stat="D",maxentargs)
niche.ovr7<-calc.niche.overlap(predicted.maps7,stat="D",maxentargs)
niche.ovr8<-calc.niche.overlap(predicted.maps8,stat="D",maxentargs)
niche.ovr9<-calc.niche.overlap(predicted.maps9,stat="D",maxentargs)
niche.ovr10<-calc.niche.overlap(predicted.maps10,stat="D",maxentargs)
niche.ovr11<-calc.niche.overlap(predicted.maps11,stat="D",maxentargs)
niche.ovrCur<-calc.niche.overlap(predicted.mapsCur, stat="D",maxentargs)
save(predicted.maps0,predicted.maps1,predicted.maps2,predicted.maps3,predicted.maps4,predicted.maps5,predicted.maps6,predicted.maps7,predicted.maps8,predicted.maps9,predicted.maps10,predicted.maps11,predicted.mapsCur,file="predictedmaps.rda")
save(niche.ovr0,niche.ovr1,niche.ovr2,niche.ovr3,niche.ovr4,niche.ovr5,niche.ovr6,niche.ovr7,niche.ovr8,niche.ovr9,niche.ovr10,niche.ovr11,niche.ovrCur, file="nicheovr.rda")

load("nicheovr.rda")
niche.ovr.mods<-c("Current", "RCP4.5INMCM42055", "RCP8.5INMCM42055", "RCP4.5INMCM42085", "RCP8.5INMCM42085", "RCP4.5CCSM42055", "RCP8.5CCSM42055", "RCP4.5CCSM42085", "RCP8.5CCSM42085", "RCP4.5GFDLCM32055", "RCP8.5GFDLCM32055", "RCP4.5GFDLCM32085", "RCP8.5GFDLCM32085")
as.numeric(niche.ovr.dat<-c(0.8334382, 0.8554097,0.800657,0.8201132,0.7829552,0.8298264,0.8340076,0.8430086,0.8062035,0.856524,0.8597929,0.8699866,0.8294891))
str(niche.ovr.dat)
min(niche.ovr.dat)
max(niche.ovr.dat)
position<-c(12,0,1,2,3,4,5,6,7,8,9,10,11)
niche.tab<-data.frame(position,niche.ovr.mods,niche.ovr.dat)
str(niche.tab)
niche.tab[1,3]-niche.tab[5,3]
save(list=ls(), file="the_future.rda")


#OKAY MADE PROJECTION MAPS FOR EACH PREDICTED FUTURE MODEL, WHATS NEXT. ENMeval package
#######################################Pinus strobus ##ENMevaluate
PSeval0<-ENMevaluate(occ=Pinus.occ, env=RCP4.5INMCM42055, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval1<-ENMevaluate(occ=Pinus.occ, env=RCP8.5INMCM42055, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval2<-ENMevaluate(occ=Pinus.occ, env=RCP4.5INMCM42085, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval3<-ENMevaluate(occ=Pinus.occ, env=RCP8.5INMCM42085, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval4<-ENMevaluate(occ=Pinus.occ, env=RCP4.5CCSM42055, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval5<-ENMevaluate(occ=Pinus.occ, env=RCP8.5CCSM42055, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval6<-ENMevaluate(occ=Pinus.occ, env=RCP4.5CCSM42085, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval7<-ENMevaluate(occ=Pinus.occ, env=RCP8.5CCSM42085, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval8<-ENMevaluate(occ=Pinus.occ, env=RCP4.5GFDLCM32055, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval9<-ENMevaluate(occ=Pinus.occ, env=RCP8.5GFDLCM32055, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval10<-ENMevaluate(occ=Pinus.occ, env=RCP4.5GFDLCM32085, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSeval11<-ENMevaluate(occ=Pinus.occ, env=RCP8.5GFDLCM32085, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
PSevalCur<-ENMevaluate(occ=Pinus.occ, env=CurrentNormal, RMvalues=seq(1,2), fc=c("L","LQ"), method='checkerboard2')

save(PSeval0,PSeval1,PSeval2,PSeval3,PSeval4,PSeval5,PSeval6,PSeval7,PSeval8,PSeval9,PSeval10,PSeval11,PSevalCur,file="PSevals.rda")

mean(PSevalCur@results$train.AUC)
PPA0<-min(PSeval0@results[4])
PPA1<-min(PSeval1@results[4])
PPA2<-min(PSeval2@results[4])
PPA3<-min(PSeval3@results[4])
PPA4<-min(PSeval4@results[4])
PPA5<-min(PSeval5@results[4])
PPA6<-min(PSeval6@results[4])
PPA7<-min(PSeval7@results[4])
PPA8<-min(PSeval8@results[4])
PPA9<-min(PSeval9@results[4])
PPA10<-min(PSeval10@results[4])
PPA11<-min(PSeval11@results[4])
PPavgAUC<-mean(PPA0,PPA1,PPA2,PPA3,PPA4,PPA5,PPA6,PPA7,PPA8,PPA9,PPA10,PPA11)
PPavgAUC

PS0<-c(PSeval0@results[1],PSeval0@results[13]) #check to see which model is best 
PS0.p<-c(as.numeric(min(PS0$AICc)),"LQ1")
PS1<-c(PSeval1@results[1],PSeval1@results[13])
PS1.p<-c(as.numeric(min(SS1$AICc)),"LQ1")
PS2<-c(PSeval2@results[1],PSeval2@results[13])
PS2.p<-c(as.numeric(min(PS2$AICc)),"LQ1")
PS3<-c(PSeval3@results[1],PSeval3@results[13])
PS3.p<-c(as.numeric(min(PS3$AICc)),"LQ1")
PS4<-c(PSeval4@results[1],PSeval4@results[13])
PS4.p<-c(as.numeric(min(PS4$AICc)),"LQ1")
PS5<-c(PSeval5@results[1],PSeval5@results[13])
PS5.p<-c(as.numeric(min(PS5$AICc)),"LQ1")
PS6<-c(PSeval6@results[1],PSeval6@results[13])
PS6.p<-c(as.numeric(min(PS6$AICc)),"LQ1")
PS7<-c(PSeval7@results[1],PSeval7@results[13])
PS7.p<-c(as.numeric(min(PS7$AICc)),"LQ1")
PS8<-c(PSeval8@results[1],PSeval8@results[13])
PS8.p<-c(as.numeric(min(PS8$AICc)),"LQ1")
PS9<-c(PSeval9@results[1],PSeval9@results[13])
PS9.p<-c(as.numeric(min(PS9$AICc)),"LQ1")
PS10<-c(PSeval10@results[1],PSeval10@results[13])
PS10.p<-c(as.numeric(min(PS10$AICc)),"LQ1")
PS11<-c(PSeval11@results[1],PSeval11@results[13])
PS11.p<-c(as.numeric(min(PS11$AICc)),"LQ1")
PSCur<-c(PSevalCur@results[1],PSevalCur@results[13])
PSCur.p<-c(as.numeric(min(PSCur$AICc)),"LQ1")
PSAICdat<-rbind(PS0.p,PS1.p,PS2.p,PS3.p,PS4.p,PS5.p,PS6.p,PS7.p,PS8.p,PS9.p,PS10.p,PS11.p,PSCur.p)

#PSeval11@results
#modelrez<-list(PSeval0@results,PSeval1@results,PSeval2@results,PSeval3@results,PSeval4@results,PSeval5@results,PSeval6@results,PSeval7@results,PSeval8@results,PSeval9@results,PSeval10@results,PSeval11@results,PSevalCur@results)
#write.csv(modelrez,file="PSevals.csv")
######################################Suillus spraguei ##ENMevaluate
SSeval0<-ENMevaluate(occ=S.sprag.occ, env=RCP4.5INMCM42055,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval1<-ENMevaluate(occ=S.sprag.occ, env=RCP8.5INMCM42055,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval2<-ENMevaluate(occ=S.sprag.occ, env=RCP4.5INMCM42085,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval3<-ENMevaluate(occ=S.sprag.occ, env=RCP8.5INMCM42085,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval4<-ENMevaluate(occ=S.sprag.occ, env=RCP4.5CCSM42055,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval5<-ENMevaluate(occ=S.sprag.occ, env=RCP8.5CCSM42055,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval6<-ENMevaluate(occ=S.sprag.occ, env=RCP4.5CCSM42085,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval7<-ENMevaluate(occ=S.sprag.occ, env=RCP8.5CCSM42085,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval8<-ENMevaluate(occ=S.sprag.occ, env=RCP4.5GFDLCM32055,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval9<-ENMevaluate(occ=S.sprag.occ, env=RCP8.5GFDLCM32055,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval10<-ENMevaluate(occ=S.sprag.occ, env=RCP4.5GFDLCM32085,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSeval11<-ENMevaluate(occ=S.sprag.occ, env=RCP8.5GFDLCM32085,  RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
SSevalCur<-ENMevaluate(occ=S.sprag.occ, env=CurrentNormal, RMvalues= seq(1,2), fc=c("L","LQ"), method='checkerboard2')
save(SSeval0,SSeval1,SSeval2,SSeval3,SSeval4,SSeval5,SSeval6,SSeval7,file="SSevals0to7.rda")
save(SSeval8,SSeval9,SSeval10,SSeval11,SSevalCur,file="SSevals8to11.rda")

mean(SSevalCur@results$train.AUC)

SSA0<-min(SSeval0@results[4])
SSA1<-min(SSeval1@results[4])
SSA2<-min(SSeval2@results[4])
SSA3<-min(SSeval3@results[4])
SSA4<-min(SSeval4@results[4])
SSA5<-min(SSeval5@results[4])
SSA6<-min(SSeval6@results[4])
SSA7<-min(SSeval7@results[4])
SSA8<-min(SSeval8@results[4])
SSA9<-min(SSeval9@results[4])
SSA10<-min(SSeval10@results[4])
SSA11<-min(SSeval11@results[4])
SSavgAUC<-mean(SSA0,SSA1,SSA2,SSA3,SSA4,SSA5,SSA6,SSA7,SSA8,SSA9,SSS10,SSA11)
SSavgAUC

SS0<-c(SSeval0@results[1],SSeval0@results[13]) #check to see which model is best 
SS0.p<-c(as.numeric(min(SS0$AICc)),"LQ1") #save best model from evals here
SS1<-c(SSeval1@results[1],SSeval1@results[13])
SS1.p<-c(as.numeric(min(SS1$AICc)),"LQ1")
SS2<-c(SSeval2@results[1],SSeval2@results[13])
SS2.p<-c(as.numeric(min(SS2$AICc)),"LQ1")
SS3<-c(SSeval3@results[1],SSeval3@results[13])
SS3.p<-c(as.numeric(min(SS2$AICc)),"LQ1")
SS4<-c(SSeval4@results[1],SSeval4@results[13])
SS4.p<-c(as.numeric(min(SS2$AICc)),"LQ1")
SS5<-c(SSeval5@results[1],SSeval5@results[13])
SS5.p<-c(as.numeric(min(SS5$AICc)),"LQ1")
SS6<-c(SSeval6@results[1],SSeval6@results[13])
SS6.p<-c(as.numeric(min(SS6$AICc)),"LQ1")
SS7<-c(SSeval7@results[1],SSeval7@results[13])
SS7.p<-c(as.numeric(min(SS7$AICc)),"LQ1")
SS8<-c(SSeval8@results[1],SSeval8@results[13])
SS8.p<-c(as.numeric(min(SS8$AICc)),"LQ1")
SS9<-c(SSeval9@results[1],SSeval9@results[13])
SS9.p<-c(as.numeric(min(SS9$AICc)),"LQ1")
SS10<-c(SSeval10@results[1],SSeval10@results[13])
SS10.p<-c(as.numeric(min(SS10$AICc)),"LQ1")
SS11<-c(SSeval11@results[1],SSeval11@results[13])
SS11.p<-c(as.numeric(min(SS11$AICc)),"LQ1")
SSCur<-c(SSevalCur@results[1],SSevalCur@results[13])
SSCur.p<-c(as.numeric(min(SSCur$AICc)),"LQ1")
SSAICdat<-rbind(SS0.p,SS1.p,SS2.p,SS3.p,SS4.p,SS5.p,SS6.p,SS7.p,SS8.p,SS9.p,SS10.p,SS11.p,SSCur.p)

modelres<-list(SSeval0@results,SSeval1@results,SSeval2@results,SSeval3@results,SSeval4@results,SSeval5@results,SSeval6@results,SSeval7@results,SSeval8@results,SSeval9@results,SSeval10@results,SSeval11@results,SSevalCur@results)
write.csv(modelres,file="SSevals.csv")

