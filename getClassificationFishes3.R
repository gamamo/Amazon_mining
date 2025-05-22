library(terra)
library(tidyterra)
library(tictoc)
library(here)
library(fs)
library(tidyverse)

# get envidata
ri <- rast("AmazonContamFromMines/AllDownstream.tif")
ri <- terra::as.polygons(ri)

ri2 <- rast("wetlands2_am_wgsUTM20S.tif")
am <- vect("amazon_basin_minewatchUTM20S.shp")

project(ri2,ri)

#get metals
as <- rast("UpdatedRastersContamination/As.tif")
as[as ==0] <- NA
as_vect <- terra::as.polygons(as)

as1 <- as
as1[as1 >1] <- NA
as1_vect <- as.polygons(as1)

as2 <- as
as2[as2 != 2] <- NA
as2_vect <- as.polygons(as2)

as3 <- as
as3[as3 != 3] <- NA
as3_vect <- as.polygons(as3)

cu <- rast("UpdatedRastersContamination/Cu.tif")
cu[cu ==0] <- NA
cu_vect <- terra::as.polygons(cu)

cu1 <- cu
cu1[cu1 >1] <- NA
cu1_vect <- as.polygons(cu1)

cu2 <- cu
cu2[cu2 != 2] <- NA
cu2_vect <- as.polygons(cu2)

cu3 <- cu
cu3[cu3 != 3] <- NA
cu3_vect <- as.polygons(cu3)

pb <- rast("UpdatedRastersContamination/Pb.tif")
pb[pb ==0] <- NA
pb_vect <- terra::as.polygons(pb)

pb1 <- pb
pb1[pb1 >1] <- NA
pb1_vect <- as.polygons(pb1)

pb2 <- pb
pb2[pb2 != 2] <- NA
pb2_vect <- as.polygons(pb2)

pb3 <- pb
pb3[pb3 != 3] <- NA
pb3_vect <- as.polygons(pb3)

zn <- rast("UpdatedRastersContamination/Zn.tif")
zn[zn ==0] <- NA
zn_vect <- terra::as.polygons(zn)

zn1 <- zn
zn1[zn1 >1] <- NA
zn1_vect <- as.polygons(zn1)

zn2 <- zn
zn2[zn2 != 2] <- NA
zn2_vect <- as.polygons(zn2)

zn3 <- zn
zn3[zn3 != 3] <- NA
zn3_vect <- as.polygons(zn3)

hg <- rast("UpdatedRastersContamination/Hg.tif")
hg[hg ==0] <- NA
hg_vect <- terra::as.polygons(hg)

hg1 <- hg
hg1[hg1 >1] <- NA
hg1_vect <- as.polygons(hg1)

hg2 <- hg
hg2[hg2 != 2] <- NA
hg2_vect <- as.polygons(hg2)

hg3 <- hg
hg3[hg3 != 3] <- NA
hg3_vect <- as.polygons(hg3)



#make the result table
comAs <- data.frame(species=NA, areaTotal = NA, areaRivers=NA, area1=NA, area2=NA, area3=NA)
comCu <- data.frame(species=NA, areaTotal = NA, areaRivers=NA, area1=NA, area2=NA, area3=NA)
comPb <- data.frame(species=NA, areaTotal = NA, areaRivers=NA, area1=NA, area2=NA, area3=NA)
comZn <- data.frame(species=NA, areaTotal = NA, areaRivers=NA, area1=NA, area2=NA, area3=NA)
comHg <- data.frame(species=NA, areaTotal = NA, areaRivers=NA, area1=NA, area2=NA, area3=NA)

#get fish species
fishes <- list.files(path = paste0("Proy_Macro/range_fishes_AmazonasTocantins_2024"),
                     pattern = c("*.shp"),recursive = T)

tic()
#for (i in 1:length(fishes)) {
  for (i in 1:10) {
  if(i == 1){next}
  if(i == 2){next}
print(i)

sp <- vect(paste0("Proy_Macro/range_fishes_AmazonasTocantins_2024/",fishes[i]))
if(length(names(sp))>1){
  sp <- sp |> select(Species) 
  names(sp) <- gsub(" ","_",sp$Species)
} else{
#names(sp) <- str_extract(fishes[i], regex("[^/]+\\.shp$"))
names(sp) <- str_extract(fishes2[z], regex('[^/]+.shp$'))
names(sp) <- gsub(".shp","",names(sp))
}
sp <- project(sp, "epsg:32720")
sp0<- aggregate(sp, by=NULL, dissolve=TRUE, fun="mean", count=TRUE)

#crop to the Amazon limits
sp1 <- crop(sp0,am)
sp1 <- terra::crop(ri2,sp1,mask=T)
sp1 <- as.polygons(sp1)


if(length(sp1)==0) {next} #if the species is not withing amazonian limits
if(terra::nrow(sp1)==0){next} #if the species has zero attributes

#clip shapes
c <- crop(sp1,ri)
#if(length(c)>0){


c1as <- terra::crop(sp1,as1_vect)
c2as <- terra::crop(sp1,as2_vect)
c3as <- terra::crop(sp1,as3_vect)

c1cu <- terra::crop(sp1,cu1_vect)
c2cu <- terra::crop(sp1,cu2_vect)
c3cu <- terra::crop(sp1,cu3_vect)

c1pb <- terra::crop(sp1,pb1_vect)
c2pb <- terra::crop(sp1,pb2_vect)
c3pb <- terra::crop(sp1,pb3_vect)

c1zn <- terra::crop(sp1,zn1_vect)
c2zn <- terra::crop(sp1,zn2_vect)
c3zn <- terra::crop(sp1,zn3_vect)

c1hg <- terra::crop(sp1,hg1_vect)
c2hg <- terra::crop(sp1,hg2_vect)
c3hg <- terra::crop(sp1,hg3_vect)

#save results
comAs[i,"species"] <- names(sp)
comAs[i,"areaTotal"] <- sum(expanse(sp1,unit="km"))
comAs[i,"areaRivers"] <- (sum(expanse(c,unit="km"))*100)/sum(expanse(sp1,unit="km"))
comAs[i,"area1"] <- if(length(c1as)>0){sum(expanse(c1as,unit="km"))} else{
  comAs[i,"area1"] <- 0
}
comAs[i,"area2"] <- if(length(c2as)>0){sum(expanse(c2as,unit="km"))} else{
  comAs[i,"area2"] <- 0
}
comAs[i,"area3"] <- if(length(c3as)>0){sum(expanse(c3as,unit="km"))} else{
  comAs[i,"area3"] <- 0
}
comAs[i,"areaTotalMetal"] <- sum(expanse(as1_vect,unit="km"))+ sum(expanse(as2_vect,unit="km"))+
  sum(expanse(as3_vect,unit="km"))

comCu[i,"species"] <- names(sp)
comCu[i,"areaTotal"] <- sum(expanse(sp1,unit="km"))
comCu[i,"areaRivers"] <- (sum(expanse(c,unit="km"))*100)/sum(expanse(sp1,unit="km"))
comCu[i,"area1"] <- if(length(c1cu)>0){sum(expanse(c1cu,unit="km"))} else{
  comCu[i,"area1"] <- 0
}
comCu[i,"area2"] <- if(length(c2cu)>0){sum(expanse(c2cu,unit="km"))} else{
  comCu[i,"area2"] <- 0
}
comCu[i,"area3"] <- if(length(c3cu)>0){sum(expanse(c3cu,unit="km"))} else{
  comCu[i,"area3"] <- 0
}
comCu[i,"areaTotalMetal"] <- sum(expanse(cu1_vect,unit="km"))+ sum(expanse(cu2_vect,unit="km"))+
  sum(expanse(cu3_vect,unit="km"))

comPb[i,"species"] <- names(sp)
comPb[i,"areaTotal"] <- sum(expanse(sp1,unit="km"))
comPb[i,"areaRivers"] <- (sum(expanse(c,unit="km"))*100)/sum(expanse(sp1,unit="km"))
comPb[i,"area1"] <- if(length(c1pb)>0){sum(expanse(c1pb,unit="km"))} else{
  comPb[i,"area1"] <- 0
}
comPb[i,"area2"] <- if(length(c2pb)>0){sum(expanse(c2pb,unit="km"))} else{
  comPb[i,"area2"] <- 0
}
comPb[i,"area3"] <- if(length(c3pb)>0){sum(expanse(c3pb,unit="km"))} else{
  comPb[i,"area3"] <- 0
}
comPb[i,"areaTotalMetal"] <- sum(expanse(pb1_vect,unit="km"))+ sum(expanse(pb2_vect,unit="km"))+
  sum(expanse(pb3_vect,unit="km"))

comZn[i,"species"] <- names(sp)
comZn[i,"areaTotal"] <- sum(expanse(sp1,unit="km"))
comZn[i,"areaRivers"] <- (sum(expanse(c,unit="km"))*100)/sum(expanse(sp1,unit="km"))
comZn[i,"area1"] <- if(length(c1zn)>0){sum(expanse(c1zn,unit="km"))} else{
  comZn[i,"area1"] <- 0
}
comZn[i,"area2"] <- if(length(c2zn)>0){sum(expanse(c2zn,unit="km"))} else{
  comZn[i,"area2"] <- 0
}
comZn[i,"area3"] <- if(length(c3zn)>0){sum(expanse(c3zn,unit="km"))} else{
  comZn[i,"area3"] <- 0
}
comZn[i,"areaTotalMetal"] <- sum(expanse(zn1_vect,unit="km"))+ sum(expanse(zn2_vect,unit="km"))+
  sum(expanse(zn3_vect,unit="km"))


comHg[i,"species"] <- names(sp)
comHg[i,"areaTotal"] <- sum(expanse(sp1,unit="km"))
comHg[i,"areaRivers"] <- (sum(expanse(c,unit="km"))*100)/sum(expanse(sp1,unit="km"))
comHg[i,"area1"] <- if(length(c1hg)>0){sum(expanse(c1hg,unit="km"))} else{
  comHg[i,"area1"] <- 0
}
comHg[i,"area2"] <- if(length(c2hg)>0){sum(expanse(c2hg,unit="km"))} else{
  comHg[i,"area2"] <- 0
}
comHg[i,"area3"] <- if(length(c3hg)>0){sum(expanse(c3hg,unit="km"))} else{
  comHg[i,"area3"] <- 0
}
comHg[i,"areaTotalMetal"] <- sum(expanse(hg1_vect,unit="km"))+ sum(expanse(hg2_vect,unit="km"))+
  sum(expanse(hg3_vect,unit="km"))


#} else {
if(F){
  
  #save results
  comAs[i,"species"] <- names(sp)
  comAs[i,"areaTotal"] <- expanse(sp1,unit="km")
  comAs[i,"areaRivers"] <- 0
  comAs[i,"area1"] <- 0
  comAs[i,"area2"] <- 0
  comAs[i,"area3"] <- 0
  
  comCu[i,"species"] <- names(sp)
  comCu[i,"areaTotal"] <- expanse(sp1,unit="km")
  comCu[i,"areaRivers"] <- 0
  comCu[i,"area1"] <- 0
  comCu[i,"area2"] <- 0
  comCu[i,"area3"] <- 0
  
  comPb[i,"species"] <- names(sp)
  comPb[i,"areaTotal"] <- expanse(sp1,unit="km")
  comPb[i,"areaRivers"] <- 0
  comPb[i,"area1"] <- 0
  comPb[i,"area2"] <- 0
  comPb[i,"area3"] <- 0
  
  comZn[i,"species"] <- names(sp)
  comZn[i,"areaTotal"] <- expanse(sp1,unit="km")
  comZn[i,"areaRivers"] <- 0
  comZn[i,"area1"] <- 0
  comZn[i,"area2"] <- 0
  comZn[i,"area3"] <- 0
  
  comHg[i,"species"] <- names(sp)
  comHg[i,"areaTotal"] <- expanse(sp1,unit="km")
  comHg[i,"areaRivers"] <- 0
  comHg[i,"area1"] <- 0
  comHg[i,"area2"] <- 0
  comHg[i,"area3"] <- 0
  
}
#}
write_csv(comAs, here("resu_fishes", "rivercrop","resuPEL_fishesAs.csv"))
write_csv(comCu, here("resu_fishes", "rivercrop","resuPEL_fishesCu.csv"))
write_csv(comPb, here("resu_fishes", "rivercrop","resuPEL_fishesPb.csv"))
write_csv(comZn, here("resu_fishes", "rivercrop","resuPEL_fishesZn.csv"))
write_csv(comHg, here("resu_fishes", "rivercrop","resuPEL_fishesHg.csv"))
}

toc()

