#Trend of Drought Frequency

# Load data
spei_data <- read.csv("./SPEI_3month.csv")
# Drought frequency
drought_threshold <- -1
window_size <- 10 * 12  # 10 years, assuming 12 months per year
step <- 12  # 1 year step

drought_frequency <- numeric()

for (start in seq(1, length(spei_data) - window_size + 1, by = step)) {
  window_data <- spei_data[start:(start + window_size - 1)]
  
  drought_events <- sum(window_data < drought_threshold, na.rm = TRUE)
  drought_frequency <- c(drought_frequency, drought_events / window_size)
}

years <- seq(1, length(drought_frequency))  
drought_frequency_ts <- ts(drought_frequency, start = 1981, frequency = 1)  

years <- seq(1981, 2011, by = 1)  
drought_data <- data.frame(
  Year = years,
  Drought_Frequency = drought_frequency
)

write.csv(drought_data,"./drought_frequency_trend.csv")

model <- lm(Drought_Frequency ~ Year, data = drought_data)
slope <- coef(model)[2] 
r_squared <- summary(model)$r.squared  

trend::mk.test(drought_frequency)

# Trend of NDVI-SPEI correlation
modis_ndvi <- stack("./MOD13A3/r_ndvi.tif")
uyz_ndvi <- mask(modis_ndvi,uyz_shp)
ndvi_data <- cellStats(uyz_ndvi,"mean",na.rm=T)
spei_data <- spei_data[13:480]


library(forcast)
library(zoo)  
library(ggplot2)

spei_ts <- ts(spei_data,start = c(1982,1),frequency = 12)
ndvi_ts <- ts(ndvi_data, start = c(1982,1),frequency = 12)  
ndvi_stl <- stl(ndvi_ts, s.window = "periodic")
ndvi_detrended <- ndvi_stl$time.series[, "remainder"]
spei_data2 <- spei_ts[218:444]

window_size <- 120
calc_rolling_correlation <- function(x, y, window) {
  rollapplyr(1:length(x), width = window, by = 1, FUN = function(idx) {
    cor(x[idx], y[idx], use = "complete.obs", method = "spearman")
  }, align = "left", fill = NA)
}

rolling_correlation <- calc_rolling_correlation(spei_data2, ndvi_detrended, window_size)
years <- seq(from = 1982, length.out = length(rolling_correlation), by = 1/12)  
cor_data <- data.frame(Year = years, Correlation = rolling_correlation)

trend::mk.test(rolling_correlation[1:108])

write.csv(cor_data[1:108,],"./spei_ndvi_corr_trend.csv")


