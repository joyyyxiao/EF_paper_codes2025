library(raster)
library(forcast)
library(zoo)  

# Load data
gimms_ndvi <- stack("./gimms_ndvi_resample_00-22.nc")
modis_ndvi <- stack("./MOD13A3/r_ndvi.tif")

# Function to calculate NDVI Anomaly for a single pixel
calculate_anomaly <- function(ndvi_data) {
  if (all(is.na(ndvi_data))) {
    return(rep(NA, length(ndvi_data)))  # Handle pixels with all NA values
  }
  if (any(is.na(ndvi_data))) {
    ndvi_data <- na.approx(ndvi_data, rule = 2) # Handle pixels with NA values
  }

  # Decompose
  ndvi_ts <- ts(ndvi_data, start = c(2000,2),frequency = 12)  
  ndvi_stl <- stl(ndvi_ts, s.window = "periodic",robust = TRUE)
  ndvi_detrended <- ndvi_stl$time.series[, "remainder"]
  
  # z-score
  ndvi_anomaly <- (ndvi_detrended - mean(ndvi_detrended, na.rm = TRUE)) / sd(ndvi_detrended, na.rm = TRUE)
  
  # Return the Anomaly values
  return(as.numeric(ndvi_anomaly))
}

# Apply Anomaly calculation to each pixel in the RasterStack
Ano_stack <- calc(modis_ndvi, fun = function(x) calculate_anomaly(x))
Ano_stack <- calc(gimms_ndvi, fun = function(x) calculate_anomaly(x))

# Save the results as a new raster stack
writeRaster(Ano_stack, "MODIS_ndvi_anomaly.tif", format = "GTiff", overwrite = TRUE)
writeRaster(Ano_stack, "GIMMS_ndvi_anomaly.tif", format = "GTiff", overwrite = TRUE)


