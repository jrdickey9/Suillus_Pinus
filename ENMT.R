#enmtools
#load(file="the_future.rda")
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
<<<<<<< HEAD
install.packages("githubinstall") #new line
library(githubinstall) #new line
gh_install_packages("ENMTools", ref="develop") #downloading from develop branch 
=======
install_github("danlwarren/ENMTools")
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
library(raster)
library(dismo)
library(rJava)
library(usdm)
require(rgeos)
require(rJava)
library(sp)
library(rgdal)
<<<<<<< HEAD
=======

>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
library(ENMeval)
library(ENMTools)
#bricks
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

#occ data for enmtools #restructuring 
Pinus.occ<-read.csv("newPinus.csv")
head(Pinus.occ)
names(Pinus.occ)
summary(Pinus.occ)
str(Pinus.occ)
Pinus.occ<-Pinus.occ[,3:4]
Pinus.occ<-Pinus.occ[!is.na(rowSums(extract(RCP4.5INMCM42055,Pinus.occ))),]
coordinates(object=Pinus.occ)<-~ Lon+Lat
Pinus.occ1<-as.data.frame(Pinus.occ)
<<<<<<< HEAD
save(Pinus.occ1,file="Pinus.occ1.rda")
=======
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665

suillydat<-read.csv(file="Suilly_occ.csv")
set.seed(150) #set this 
trueSuilly.occ<-suillydat[sample(nrow(suillydat), 54), ] #check to make sure n matches Pinus occ
trueSuilly.occ<-trueSuilly.occ[,1:2] #selecting out long lats, make sure data is in deg.dec
S.sprag.occ<-trueSuilly.occ[!is.na(rowSums(extract(RCP4.5INMCM42055,trueSuilly.occ))),]
coordinates(trueSuilly.occ) <-~ Long + Lat #these should be headers
S.sprag.occ1<-as.data.frame(S.sprag.occ)
<<<<<<< HEAD
save(S.sprag.occ1,file="S.sprag.occ1.rda")
=======
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
q1<-"Pinus_strobus" #creating character vector names for id.tests 
q2<-"Suillus_spraguei"

#RCP4.5INMCM42055
Pinus<-enmtools.species(RCP4.5INMCM42055[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP4.5INMCM42055[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test1<-identity.test(Pinus,Suillus,RCP4.5INMCM42055,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test1, file="id.test0.rda")
ps.bg.test1<-background.test(Pinus,Suillus,RCP4.5INMCM42055, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test1<-identity.test(Pinus,Suillus,RCP4.5INMCM42055,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test1, file="id.test0.rda")
ps.bg.test1<-background.test(Pinus,Suillus,RCP4.5INMCM42055, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test1, file="bg.test0.rda")

#RCP8.5INMCM42055
Pinus<-enmtools.species(RCP8.5INMCM42055[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP8.5INMCM42055[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test2<-identity.test(Pinus,Suillus,RCP8.5INMCM42055,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test2, file="id.test1.rda")
ps.bg.test2<-background.test(Pinus,Suillus,RCP8.5INMCM42055, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur.
=======
ps.id.test2<-identity.test(Pinus,Suillus,RCP8.5INMCM42055,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test2, file="id.test1.rda")
ps.bg.test2<-background.test(Pinus,Suillus,RCP8.5INMCM42055, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur.
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test2, file="bg.test1.rda")

#RCP4.5INMCM42085
Pinus<-enmtools.species(RCP4.5INMCM42085[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP4.5INMCM42085[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test3<-identity.test(Pinus,Suillus,RCP4.5INMCM42085,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test3, file="id.test2.rda")
ps.bg.test3<-background.test(Pinus,Suillus,RCP4.5INMCM42085, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur.
=======
ps.id.test3<-identity.test(Pinus,Suillus,RCP4.5INMCM42085,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test3, file="id.test2.rda")
ps.bg.test3<-background.test(Pinus,Suillus,RCP4.5INMCM42085, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur.
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test3, file="bg.test2.rda")

#RCP8.5INMCM42085
Pinus<-enmtools.species(RCP8.5INMCM42085[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP8.5INMCM42085[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test4<-identity.test(Pinus,Suillus,RCP8.5INMCM42085,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test4, file="id.test3.rda")
ps.bg.test4<-background.test(Pinus,Suillus,RCP8.5INMCM42085, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test4<-identity.test(Pinus,Suillus,RCP8.5INMCM42085,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test4, file="id.test3.rda")
ps.bg.test4<-background.test(Pinus,Suillus,RCP8.5INMCM42085, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test4, file="bg.test3.rda")

#RCP4.5CCSM42055
Pinus<-enmtools.species(RCP4.5CCSM42055[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP4.5CCSM42055[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test5<-identity.test(Pinus,Suillus,RCP4.5CCSM42055,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test5, file="id.test4.rda")
ps.bg.test5<-background.test(Pinus,Suillus,RCP4.5CCSM42055, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test5<-identity.test(Pinus,Suillus,RCP4.5CCSM42055,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test5, file="id.test4.rda")
ps.bg.test5<-background.test(Pinus,Suillus,RCP4.5CCSM42055, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test5, file="bg.test4.rda")

#RCP8.5CCSM42055
Pinus<-enmtools.species(RCP8.5CCSM42055[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP8.5CCSM42055[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test6<-identity.test(Pinus,Suillus,RCP8.5CCSM42055,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test6, file="id.test5.rda")
ps.bg.test6<-background.test(Pinus,Suillus,RCP8.5CCSM42055, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test6<-identity.test(Pinus,Suillus,RCP8.5CCSM42055,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test6, file="id.test5.rda")
ps.bg.test6<-background.test(Pinus,Suillus,RCP8.5CCSM42055, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test6, file="bg.test5.rda")

#RCP4.5CCSM42085
Pinus<-enmtools.species(RCP4.5CCSM42085[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP4.5CCSM42085[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test7<-identity.test(Pinus,Suillus,RCP4.5CCSM42085,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test7, file="id.test6.rda")
ps.bg.test7<-background.test(Pinus,Suillus,RCP4.5CCSM42085, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test7<-identity.test(Pinus,Suillus,RCP4.5CCSM42085,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test7, file="id.test6.rda")
ps.bg.test7<-background.test(Pinus,Suillus,RCP4.5CCSM42085, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test7, file="bg.test6.rda")

#RCP8.5CCSM42085
Pinus<-enmtools.species(RCP8.5CCSM42085[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP8.5CCSM42085[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test8<-identity.test(Pinus,Suillus,RCP8.5CCSM42085,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test8, file="id.test7.rda")
ps.bg.test8<-background.test(Pinus,Suillus,RCP8.5CCSM42085, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test8<-identity.test(Pinus,Suillus,RCP8.5CCSM42085,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test8, file="id.test7.rda")
ps.bg.test8<-background.test(Pinus,Suillus,RCP8.5CCSM42085, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test8, file="bg.test7.rda")

#RCP8.5GFDLCM32055
Pinus<-enmtools.species(RCP8.5GFDLCM32055[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP8.5GFDLCM32055[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test9<-identity.test(Pinus,Suillus,RCP8.5GFDLCM32055,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test9, file="id.test8.rda")
ps.bg.test9<-background.test(Pinus,Suillus,RCP8.5GFDLCM32055, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur.
=======
ps.id.test9<-identity.test(Pinus,Suillus,RCP8.5GFDLCM32055,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test9, file="id.test8.rda")
ps.bg.test9<-background.test(Pinus,Suillus,RCP8.5GFDLCM32055, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur.
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test9, file="bg.test8.rda")

#RCP4.5GFDLCM32055
Pinus<-enmtools.species(RCP4.5GFDLCM32055[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP4.5GFDLCM32055[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test10<-identity.test(Pinus,Suillus,RCP4.5GFDLCM32055,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test10, file="id.test9.rda")
ps.bg.test10<-background.test(Pinus,Suillus,RCP4.5GFDLCM32055, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test10<-identity.test(Pinus,Suillus,RCP4.5GFDLCM32055,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test10, file="id.test9.rda")
ps.bg.test10<-background.test(Pinus,Suillus,RCP4.5GFDLCM32055, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test10, file="bg.test9.rda")

#RCP4.5GFDLCM32085
Pinus<-enmtools.species(RCP4.5GFDLCM32085[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP4.5GFDLCM32085[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test11<-identity.test(Pinus,Suillus,RCP4.5GFDLCM32085,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test11, file="id.test10.rda")
ps.bg.test11<-background.test(Pinus,Suillus,RCP4.5GFDLCM32085, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test11<-identity.test(Pinus,Suillus,RCP4.5GFDLCM32085,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test11, file="id.test10.rda")
ps.bg.test11<-background.test(Pinus,Suillus,RCP4.5GFDLCM32085, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test11, file="bg.test10.rda")

#RCP8.5GFDLCM32085
Pinus<-enmtools.species(RCP8.5GFDLCM32085[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(RCP8.5GFDLCM32085[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test12<-identity.test(Pinus,Suillus,RCP8.5GFDLCM32085,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test12, file="id.test11.rda")
ps.bg.test12<-background.test(Pinus,Suillus,RCP8.5GFDLCM32085, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test12<-identity.test(Pinus,Suillus,RCP8.5GFDLCM32085,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test12, file="id.test11.rda")
ps.bg.test12<-background.test(Pinus,Suillus,RCP8.5GFDLCM32085, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test12, file="bg.test11.rda")

#CurrentNormal
Pinus<-enmtools.species(CurrentNormal[[1]], presence.points = Pinus.occ1, background.points = NA,species.name = q1, models = NA) #making enmtools species objects for id.tests
Suillus<-enmtools.species(CurrentNormal[[1]], presence.points = S.sprag.occ1,background.points=NA,species.name=q2, models=NA)
<<<<<<< HEAD
ps.id.test13<-identity.test(Pinus,Suillus,CurrentNormal,type="mx", nreps=100, low.memory=TRUE) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test13, file="id.test12.rda")
ps.bg.test13<-background.test(Pinus,Suillus,CurrentNormal, type="mx", nreps=100, low.memory=TRUE) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
=======
ps.id.test13<-identity.test(Pinus,Suillus,CurrentNormal,type="mx", nreps=100) #null hypothesis that the two species' occurrences in the environment are effectively a random draw from the same underlying distribution. Note that niche evolution is only one of many reasons why two species' realized environmental distributions might cause departures from this null hypothesis (Warren et al. 2014)
save(ps.id.test13, file="id.test12.rda")
ps.bg.test13<-background.test(Pinus,Suillus,CurrentNormal, type="mx", nreps=100) #The purpose of this test is to correct for the availability of habitat and ask whether the observed similarity between species or populations is significantly more (or less) than expected given the available set of environments in the regions in which they occur. 
>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
save(ps.bg.test13, file="bg.test12.rda")
save(list=ls(), file="ENMT.rda")

#load back in to begin analyzing data if R Session closed. #reset directory, large files may need to load a few or one at a time to view. 
<<<<<<< HEAD
=======
load("bg.test0.rda") #RCP4.5 INMCM4 2055
load("id.test0.rda") #.. 
load("bg.test1.rda") #RCP8.5 INMCM4 2055
load("id.test1.rda") # .. 
load("bg.test2.rda") #RCP4.5 INMCM4 2085
load("id.test2.rda") #..
load("bg.test3.rda") #RCP8.5 INMCM4 2085
load("id.test3.rda") #..
load("bg.test4.rda") #RCP4.5 CCSM4 2055 
load("id.test4.rda") #..
load("bg.test5.rda") #RCP8.5 CCSM4 2055 
load("id.test5.rda") #..
load("bg.test6.rda") #RCP4.5 CCSM4 2085 
load("id.test6.rda")
load("bg.test7.rda") #RCP8.5 CCSM4 2085 ..
load("id.test7.rda")
load("bg.test8.rda") #RCP8.5 GFDLCM3 2055 ..
load("id.test8.rda")
load("bg.test9.rda") #RCP4.5 GFDLCM3 2055 ..
load("id.test9.rda")
load("bg.test10.rda") #RCP4.5 GFDLCM3 2085 ..
load("id.test10.rda")
load("bg.test11.rda") #RCP8.5 GFDLCM3 2085
load("id.test11.rda")
load("bg.test12.rda") #Climate Normal
load("id.test12.rda")

>>>>>>> ba5522df1ef0cd258186c60a36367fff33ab7665
