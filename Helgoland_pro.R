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

# Combine all datasets into one
df <- bind_rows(H62, H63, H64, H65, H66, H67, H68, H69, H70, H71, H72, H73, H74, H75, H76,
                H77, H78, H79, H80, H81, H82, H83, H84, H85, H86, H87, H88, H89, H90, H91,
                H92, H93, H94, H98)

library(lubridate)

# add year as a column
df$year <- year(ymd_hm(df$Date))

