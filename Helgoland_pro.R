H62 <- read.csv("dataset/Helgoland_1962.csv")
H63 <- read.csv("dataset/Helgoland_1963.csv")
H64 <- read.csv("dataset/Helgoland_1964.csv")
H65 <- read.csv("dataset/Helgoland_1965.csv")
H66 <- read.csv("dataset/Helgoland_1966.csv")
H67 <- read.csv("dataset/Helgoland_1967.csv")
H68 <- read.csv("dataset/Helgoland_1968.csv")
H69 <- read.csv("dataset/Helgoland_1969.csv")
H70 <- read.csv("dataset/Helgoland_1970.csv")
H71 <- read.csv("dataset/Helgoland_1971.csv")
H72 <- read.csv("dataset/Helgoland_1972.csv")
H73 <- read.csv("dataset/Helgoland_1973.csv")
H74 <- read.csv("dataset/Helgoland_1974.csv")
H75 <- read.csv("dataset/Helgoland_1975.csv")
H76 <- read.csv("dataset/Helgoland_1976.csv")
H77 <- read.csv("dataset/Helgoland_1977.csv")
H78 <- read.csv("dataset/Helgoland_1978.csv")
H79 <- read.csv("dataset/Helgoland_1979.csv")
H80 <- read.csv("dataset/Helgoland_1980.csv")
H81 <- read.csv("dataset/Helgoland_1981.csv")
H82 <- read.csv("dataset/Helgoland_1982.csv")
H83 <- read.csv("dataset/Helgoland_1983.csv")
H84 <- read.csv("dataset/Helgoland_1984.csv")
H85 <- read.csv("dataset/Helgoland_1985.csv")
H86 <- read.csv("dataset/Helgoland_1986.csv")
H87 <- read.csv("dataset/Helgoland_1987.csv")
H88 <- read.csv("dataset/Helgoland_1988.csv")
H89 <- read.csv("dataset/Helgoland_1989.csv")
H90 <- read.csv("dataset/Helgoland_1990.csv")
H91 <- read.csv("dataset/Helgoland_1991.csv")
H92 <- read.csv("dataset/Helgoland_1992.csv")
H93 <- read.csv("dataset/Helgoland_1993.csv")
H94 <- read.csv("dataset/Helgoland_1994.csv")
H98 <- read.csv("dataset/Helgoland_1998.csv")
H2002 <- read.csv("dataset/Helgoland_2002.csv")

library(dplyr)
library(ggplot2)
library(patchwork)

# Combine all datasets into one
df <- bind_rows(H62, H63, H64, H65, H66, H67, H68, H69, H70, H71, H72, H73, H74, H75, H76,
                H77, H78, H79, H80, H81, H82, H83, H84, H85, H86, H87, H88, H89, H90, H91,
                H92, H93, H94, H98)

library(lubridate)

# add year as a column
df$year <- year(ymd_hm(df$Date))

# See changes in temp throughout the years 
## Calculate yearly average temperature
yearly_avg <- df %>%
  group_by(year) %>%
  summarise(avg_temp = mean(Temp, na.rm = TRUE))
### plot the average temp for each year
ggplot(yearly_avg, aes(x = year, y = avg_temp)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Temperature per Year",
       x = "Year",
       y = "Average Temperature")


# see how nitrite differ over the years 
library(dplyr)
library(e1071)  # for skewness and kurtosis
## variance
yearly_stats <- df %>%
  group_by(year) %>%
  summarise(
    mean_nitrite = mean(Nitrite, na.rm = TRUE),
    var_nitrite = var(Nitrite, na.rm = TRUE),
    skew_nitrite = skewness(Nitrite, na.rm = TRUE),
    kurt_nitrite = kurtosis(Nitrite, na.rm = TRUE),
    mean_nitrate = mean(Nitrate, na.rm = TRUE),
    var_nitrate = var(Nitrate, na.rm = TRUE),
    skew_nitrate = skewness(Nitrate, na.rm = TRUE),
    kurt_nitrate = kurtosis(Nitrate, na.rm = TRUE),
    mean_phospate = mean(Phosphate, na.rm = TRUE),
    var_phosphate = var(Phosphate, na.rm = TRUE),
    skew_phosphate = skewness(Phosphate, na.rm = TRUE),
    kurt_phosphate = kurtosis(Phosphate, na.rm = TRUE)
  )
##visualize for variance 
ni_v <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = var_nitrite), color = 'blue') +
  labs(title = "Variance of Nitrite Over Years", y = "Variance")
##visualize for skewness 
ni_s <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = skew_nitrite), color = 'blue') +
  labs(title = "Skewness of Nitrite Over Years", y = "Skewness")
##visualize for kurtosis 
ni_k <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = kurt_nitrite), color = 'blue') +
  labs(title = "Kurtosis of Nitrite Over Years", y = "Kurtosis")
###stack them together 
(ni_v / ni_s / ni_k) + plot_layout(ncol = 1)


# see how nitrate differ over the years 
##visualize for variance 
na_v <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = var_nitrate), color = 'goldenrod') +
  labs(title = "Variance of Nitrate Over Years", y = "Variance")
##visualize for skewness 
na_s <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = skew_nitrate), color = 'goldenrod') +
  labs(title = "Skewness of Nitrate Over Years", y = "Skewness")
##visualize for kurtosis 
na_k <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = kurt_nitrate), color = 'goldenrod') +
  labs(title = "Kurtosis of Nitrate Over Years", y = "Kurtosis")
###stack them together 
(na_v / na_s / na_k) + plot_layout(ncol = 1)


# see how phosphate differ over the years 
##visualize for variance 
p_v <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = var_phosphate), color = 'darkolivegreen') +
  labs(title = "Variance of Phosphate Over Years", y = "Variance")
##visualize for skewness 
p_s <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = skew_phosphate), color = 'darkolivegreen') +
  labs(title = "Skewness of Phosphate Over Years", y = "Skewness")
##visualize for kurtosis 
p_k <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = kurt_phosphate), color = 'darkolivegreen') +
  labs(title = "Kurtosis of Phosphate Over Years", y = "Kurtosis")
###stack them together 
(p_v / p_s / p_k) + plot_layout(ncol = 1)

