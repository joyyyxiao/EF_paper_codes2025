library(greenbrown)
library(raster)

# Function for calculating monthly averages
calculate_monthly_avg <- function(raster_stack) {
  # Ensure the number of layers is divisible by 12 (for monthly data)
  if (nlayers(raster_stack) %% 12 != 0) {
    stop("The raster stack does not have a number of layers divisible by 12. Please check your data.")
  }
  
  # Create a list to store monthly averages
  monthly_avg_list <- list()
  
  # Loop through each month (1 to 12)
  for (month in 1:12) {
    # Extract layers for the current month across all years
    monthly_layers <- raster_stack[[seq(month, nlayers(raster_stack), by = 12)]]
    
    # Calculate the mean for this month across all years
    monthly_avg <- calc(monthly_layers, fun = mean, na.rm = TRUE)
    
    # Add to the list
    monthly_avg_list[[month]] <- monthly_avg
  }
  
  # Combine the results into a stack
  monthly_avg_stack <- stack(monthly_avg_list)
  
  # Assign layer names to the stack
  names(monthly_avg_stack) <- paste0("Month_", 1:12)
  
  return(monthly_avg_stack)
}

m12_ndvi <- calculate_monthly_avg(ndviData)
m12_ndvi <- aggregate(m12_ndvi,fact=6,fun=mean) #aggregate to resolution 0.5 degree
m12_ndvi <- stack(m12_ndvi,m12_ndvi)
data <- array(data=NA,dim=c(360,720,24))
for (i in 1:24){
  data[,,i] <- as.matrix(subset(m12_ndvi,i))
}

###PhenoTrs function revision####
Phenfun <- function (x, approach = c("White", "Trs"), trs = NULL, min.mean = 0.1, 
                     calc.pheno = TRUE, plot = FALSE, ...) {
  if (all(is.na(x))) 
    return(c(sos = NA, eos = NA, los = NA, pop = NA, pot = NA, 
             mgs = NA, rsp = NA, rau = NA, peak = NA, trough = NA, 
             msp = NA, mau = NA))
  n <- length(x)
  avg <- mean(x, na.rm = TRUE)
  x2 <- na.omit(x)
  avg2 <- mean(x2[x2 > min.mean], na.rm = TRUE)
  peak <- max(x, na.rm = TRUE)
  trough <- min(x, na.rm = TRUE)
  ampl <- peak - trough
  pop <- median(which(x == max(x, na.rm = TRUE)))
  pot <- median(which(x == min(x, na.rm = TRUE)))
  if (!calc.pheno) {
    if (avg < min.mean) {
      return(c(sos = NA, eos = NA, los = NA, pop = NA, 
               pot = NA, mgs = NA, rsp = NA, rau = NA, peak = NA, 
               trough = NA, msp = NA, mau = NA))
    }
    else {
      return(c(sos = NA, eos = NA, los = NA, pop = pop, 
               pot = pot, mgs = avg2, rsp = NA, rau = NA, peak = peak, 
               trough = NA, msp = NA, mau = NA))
    }
  }
  approach <- approach[1]
  if (approach == "White") {
    ratio <- (x - trough)/ampl
    trs <- 0.5
    trs.low <- trs - 0.05
    trs.up <- trs + 0.05
  }
  if (approach == "Trs") {
    ratio <- (x - trough)/ampl
    a <- diff(range(ratio, na.rm = TRUE)) * 0.1
    trs.low <- trs - a
    trs.up <- trs + a
  }
  greenup <- Greenup(ratio)
  bool <- ratio >= trs.low & ratio <= trs.up
  soseos <- 1:length(x)
  sos <- median(soseos[greenup & bool], na.rm = TRUE)
  eos <- median(soseos[!greenup & bool], na.rm = TRUE)
  los <- eos - sos
  los[los < 0] <- n + (eos[los < 0] - sos[los < 0])
  mgs <- mean(x[ratio > trs], na.rm = TRUE)
  msp <- mau <- NA
  if (!is.na(sos)) {
    id <- (sos - 10):(sos + 10)
    id <- id[(id > 0) & (id < n)]
    msp <- mean(x[id], na.rm = TRUE)
  }
  if (!is.na(eos)) {
    id <- (eos - 10):(eos + 10)
    id <- id[(id > 0) & (id < n)]
    mau <- mean(x[id], na.rm = TRUE)
  }
  metrics <- c(sos = sos, eos = eos, los = los, pop = pop, 
               pot = pot, mgs = mgs, rsp = NA, rau = NA, peak = peak, 
               trough = trough, msp = msp, mau = mau)
  if (plot) {
    if (approach == "White") 
      PlotPhenCycle(x, ratio, metrics = metrics, trs = trs, 
                    ...)
    if (approach == "Trs") 
      PlotPhenCycle(ratio, metrics = metrics, trs = trs, 
                    ...)
  }
  return(metrics)
}


###Phen maps (by year)###

for(year in 1:23){
  
  phenarr <- array(data=NA,dim=c(360,720,12))  
  for (x in 1:360){
    for(y in 1:720){
      ndvi <- ndviarr[x,y,] 
      
      if(all(is.na(ndvi))){
        phen <-c(sos = NA, eos = NA, los = NA, pop = NA, pot = NA, 
                 mgs = NA, rsp = NA, rau = NA, peak = NA, trough = NA, 
                 msp = NA, mau = NA)
        print(c(x,y,"NA pixel"))
      }else{
        ndvi <- ts(data = ndvi,start=c(2000, 1), freq=12)
        smooth <- TsPP(ndvi,fpg = FillPermanentGaps,interpolate = T,tsgf = TSGFspline)[((year-1)*365+1):(year*365)]
        
        if(all(is.na(smooth))){
          phen <- c(sos = NA, eos = NA, los = NA, pop = NA, pot = NA, 
                    mgs = NA, rsp = NA, rau = NA, peak = NA, trough = NA, 
                    msp = NA, mau = NA)
          print(c(x,y,"No seasonality"))
        }else{
          
          phen <- Phenfun(smooth,approach="Trs",trs = 0.25,plot=F)
          print(c(x,y,"Succeed"))
        }
        
        phenarr[x,y,] <- as.numeric(phen)
        
      }
    }
  }
  sos <- raster(phenarr[,,1],xmn=-180,xmx=180,ymn=-90,ymx=90,
                crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  
  eos <- raster(phenarr[,,2],xmn=-180,xmx=180,ymn=-90,ymx=90,
                crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  
  sosname <- paste0("./phen_rasters/sos_",1999+year,".tif",seq="")
  eosname <- paste0("./phen_rasters/eos_",1999+year,".tif",seq="")
  writeRaster(sos,sosname,overwrite=T)
  writeRaster(eos,eosname,overwrite=T)
  
  
}


###mean Phen map##
msos <- mean(stack(list.files("./phen_rasters/",pattern = "sos_",full.names = T)),na.rm=T)
meos <- mean(stack(list.files("./phen_rasters/",pattern = "eos_",full.names = T)),na.rm=T)
writeRaster(msos,"./phen_rasters/mean_sos.tif")
writeRaster(meos,"./phen_rasters/mean_eos.tif")


###Greenup###
###Greenup function revision####
gu <- function (x, ...) 
{
  ratio.deriv <- c(NA, diff(x))
  greenup <- rep(NA, length(x))
  greenup[ratio.deriv > 0] <- 1
  greenup[ratio.deriv < 0] <- 0
  return(greenup)
}

###Ratio Stage function revision####
Ratiofun <- function (x) {
  if (all(is.na(x))) 
    return(rep(NA,12))
 
  peak <- max(x, na.rm = TRUE)
  trough <- min(x, na.rm = TRUE)
  ampl <- peak - trough
  ratio <- (x - trough)/ampl
  
  return(ratio)
}

#####

as_raster <- function(x){
  mat <- matrix(x,nrow=360,ncol=720)
  r <- raster(mat, xmn=-180,xmx=180,ymn=-90,ymx=90,crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(r)
}


greenuparr <- array(NA,dim=c(360,720,12))
for (x in 1:360){
  for(y in 1:720){
    ndvi <- data[x,y,] 
    
    if(all(is.na(ndvi))){
    greenup <- rep(NA,12)
    print(c(x,y,"NA pixel"))
    }else{
    greenup <- gu(ndvi)[13:24]
    print(c(x,y,"Succeed"))
    }
    
    greenuparr[x,y,] <- greenup
  }
}

greenup_rasters <- stack(
  as_raster(greenuparr[,,1]),
  as_raster(greenuparr[,,2]),
  as_raster(greenuparr[,,3]),
  as_raster(greenuparr[,,4]),
  as_raster(greenuparr[,,5]),
  as_raster(greenuparr[,,6]),
  as_raster(greenuparr[,,7]),
  as_raster(greenuparr[,,8]),
  as_raster(greenuparr[,,9]),
  as_raster(greenuparr[,,10]),
  as_raster(greenuparr[,,11]),
  as_raster(greenuparr[,,12])
)

monthRecog <- function(x){
  if(is.na(x))
    return(NA)
  if(x<=31)
    return(1)
  if(x>31&x<=59)
    return(2)
  if(x>59&x<=90)
    return(3)
  if(x>90&x<=120)
    return(4)
  if(x>120&x<=151)
    return(5)
  if(x>151&x<=181)
    return(6)
  if(x>181&x<=212)
    return(7)
  if(x>212&x<=243)
    return(8)
  if(x>243&x<=273)
    return(9)
  if(x>273&x<=304)
    return(10)
  if(x>304&x<=334)
    return(11)
  if(x>334&x<=365)
    return(12)
} #convert day to month

rsos <- as_raster(
  matrix(sapply(phenarr[,,1],monthRecog),
         360,720)
)

reos <- as_raster(
  matrix(sapply(phenarr[,,2],monthRecog),
         360,720)
)

writeRaster(rsos,"./sos_month.tif")
writeRaster(reos,"./eos_month.tif")
writeRaster(greenup_rasters,paste0("./green_",1:12,".tif"),bylayer=T)


###growing stage define###

msos <- as.numeric(as.matrix(rsos))
meos <- as.numeric(as.matrix(reos))
mgreenup <- matrix(NA,360*720,12)
for(i in 1:12){
  mgreenup[,i] <- as.numeric(as.matrix(subset(greenup_rasters,i)))
}
mgstage <- matrix(NA,360*720,12)


#mgstage 
#99 <- out
#1 <- in+ 
#0 <- in- 
for(i in 1:259200){
  start <- msos[i]
  end <- meos[i]
  greenup <- mgreenup[i,]
  stage <- rep(99,12)
  if(is.na(start)|is.na(end)){
    stage <- rep(NA,12)
  }else if(start<end){
    stage[1:(start-1)] <- 99
    stage[(start+1):12] <- 99
    stage[start:end] <- greenup[start:end]
  }else if(start>end){
    stage[(end+1):(start-1)] <- 99
    stage[1:end] <- greenup[1:end]
    stage[start:12] <- greenup[start:12]
  }else{
    stage <- rep(NA,12)
  }
  mgstage[i,] <- stage
  print(i)
}

for (i in 1:12){
writeRaster(as_raster(mgstage[,i]),paste0("./growing_stage_",i,".tif"))
}


#ratiostage
baseNDVI <- as_array(m12_ndvi)
mgratio <- matrix(NA,259200,12)
rownames(mgratio) <- lonlat
for(i in 1:259200){
  x <- baseNDVI[i,]
  greenup <- mgreenup[i,]
  if(any(is.na(x))|any(is.na(greenup))){
    ratio <- rep(NA,12)
  }else{
    peak <- max(x, na.rm = TRUE)
    trough <- min(x, na.rm = TRUE)
    ampl <- peak - trough
    ratio <- (x - trough) / ampl
    greenup <- greenup==1
    ratio[!greenup] <- ratio[!greenup]*-1
    ratio <- round(ratio,1)
  }
  mgratio[i,] <- ratio
  print(i)
}

