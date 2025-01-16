# Function to calculate PET using Penman-Monteith method
calculate_pet <- function(T, R_s, u_2, q, P) {
  # Constants
  gamma <- 0.665e-3 * P        # Psychrometric constant (kPa/°C)
  G <- 0                       # Soil heat flux density (assumed negligible for daily time step) MJ m-2 day-1
  albedo <- 0.23               # Albedo
  
 
  e_s <- 0.6108 * exp((17.27 * T) / (T + 237.3))  # Saturated vapor pressure (kPa)
  Delta <- (4098 * e_s) / (T + 237.3)^2           # Slope of the vapor pressure curve (kPa/°C)
  
  # Actual vapor pressure (e_a, kPa)
  e_a <- q * P / (0.622 + q)                      # Using specific humidity and air pressure
  
  # Net radiation (R_n, MJ m-2 day-1)
  R_n <- (R_s) * (1 - albedo) 
  # Ensure R_s is provided as the sum of shortwave and longwave radiation
  
  # PET Calculation (Penman-Monteith formula)
  PET <- (0.408 * Delta * (R_n - G) + gamma * (900 / (T + 273)) * u_2 * (e_s - e_a)) / 
    (Delta + gamma * (1 + 0.34 * u_2))
  
  return(PET)  # PET in mm/day
}

# Load data
library(raster)
prec <- stack("./CMFD1979_2018/prec_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]] #2000-2018
pres <- stack("./CMFD1979_2018/pres_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]]
srad <- stack("./CMFD1979_2018/srad_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]]
lrad <- stack("./CMFD1979_2018/lrad_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]]
temp <- stack("./CMFD1979_2018/temp_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]]
shum <- stack("./CMFD1979_2018/shum_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]]
u_10 <- stack("./CMFD1979_2018/wind_CMFD_V0106_B-01_01mo_010deg_197901-201812.nc")[[253:480]]


uyz_shp <-shapefile("./UYZ_bound.shp") 
prec <- mask(prec,uyz_shp)
pres <- mask(pres,uyz_shp)
srad <- mask(srad,uyz_shp)
lrad <- mask(lrad,uyz_shp)
temp <- mask(temp,uyz_shp)
shum <- mask(shum,uyz_shp)
u_10 <- mask(u_10,uyz_shp)


# Unit conversions
temp_C <- temp - 273.15                     # Convert temperature from K to °C
pres_kPa <- pres / 1000                    # Convert pressure from Pa to kPa
srad_MJ <- (srad * 3600 * 3 * 8) / 1e6         # Convert shortwave radiation from W/m² (3-hourly) to MJ/m²/day
lrad_MJ <- (lrad * 3600 * 3 * 8) / 1e6         # Convert longwave radiation from W/m² (3-hourly) to MJ/m²/day
R_s <- srad_MJ + lrad_MJ              # Total radiation (shortwave + longwave)
u_2 <- u_10 * 4.87 / log( 67.8 * 10 - 5.42 ) # Convert wind speed to 2m height according to logarithmic wind speed profile by ALLEN et al. 

# Calculate PET
PET_day <- calculate_pet(temp_C,R_s, u_2,shum,pres_kPa)


# Calculate SPEI
# Unit conversions
days_of_month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
PET_month <- PET_day*rep(days_of_month,19)
prec_month <- prec*24*rep(days_of_month,19)

Balance <- PET_month-prec_month

library(SPEI)

# Function to calculate SPEI for a single pixel
calculate_spei <- function(water_balance, scale = 3) {
  if (all(is.na(water_balance))) {
    return(rep(NA, length(water_balance)))  # Handle pixels with all NA values
  }
  
  # Fit the SPEI model for the water balance series
  spei_result <- spei(water_balance, scale = scale)
  
  # Return the SPEI values
  return(as.numeric(spei_result$fitted))
}

# Apply SPEI calculation to each pixel in the RasterStack
SPEI_stack <- calc(Balance, fun = function(x) calculate_spei(x, scale = 3))

# Save the SPEI results as a new raster stack
writeRaster(SPEI_stack, "SPEI_3month.tif", format = "GTiff", overwrite = TRUE)


# Function to calculate SPI for a single pixel
calculate_spi <- function(precipitation, scale = 3) {
  if (all(is.na(precipitation))) {
    return(rep(NA, length(precipitation)))  # Handle pixels with all NA values
  }
  
  # Fit the SPI model for the precipitation series
  spi_result <- spi(precipitation, scale = scale)
  
  # Return the SPI values
  return(as.numeric(spi_result$fitted))
}

# Apply SPI calculation to each pixel in the RasterStack
SPI_stack <- calc(prec_month, fun = function(x) calculate_spi(x, scale = 3))

# Save the SPI results as a new raster stack
writeRaster(SPI_stack, "SPI_3month.tif", format = "GTiff", overwrite = TRUE)



