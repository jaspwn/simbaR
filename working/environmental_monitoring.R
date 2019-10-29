library(data.table)
library(ggplot2)
library(scales)
library(cowplot)


solar_file <- file.choose()

raw_solar <- fread(solar_file, drop = c(1,4), col.names = c("time", "solar(w/m2)"))

raw_solar$time <- as.POSIXct(raw_solar$time, format = '%y-%m-%d %H:%M:%S', tz = "GMT")
## round time to nearest whole minute
raw_solar$time <- format(round(raw_solar$time, units = "mins"), format = "%y-%m-%d %H:%M:%S")


anemometer_file <- file.choose()



raw_anemometer <- fread(anemometer_file, fill = TRUE, col.names = c("time", "windspeed", "rh", "temp"))
raw_anemometer$time <- paste(strftime(raw_solar$time[1], format = "%y-%m-%d"), raw_anemometer$time, sep = " ")

raw_anemometer$time <- as.POSIXct(raw_anemometer$time, format = '%y-%m-%d %H:%M:%S', tz = "GMT")

cuts = cut(raw_anemometer$time, breaks= "min", labels=FALSE)

raw_anemometer[, cuts := cuts]

anemometer_dt <- raw_anemometer[, .(time = min(time, na.rm = TRUE),
                                    windspeed = mean(windspeed, na.rm = TRUE),
                                    rh = mean(rh, na.rm = TRUE),
                                    temp = mean(temp, na.rm = TRUE)),
                                by = cuts]

anemometer_dt[, cuts := NULL]
## round time to nearest whole minute
anemometer_dt$time <- format(round(anemometer_dt$time, units = "mins"), format = "%y-%m-%d %H:%M:%S")

## merge solar and anemometer data by Time
dt <- merge(raw_solar, anemometer_dt, by = "time")
dt$time <- as.POSIXct(dt$time, format = "%y-%m-%d %H:%M:%S", tz = "GMT")

## plot data

## melt data

molten_dt <- melt(dt, id = c('time'), variable.name = 'env.cond', value.name = 'value')

env_plot <- ggplot(data = molten_dt, aes(x = time, y = value, colour = env.cond, group = env.cond)) +
  geom_line(size = 1) +
  scale_x_datetime(breaks = scales::date_breaks(width = '5 min'),
                   expand = c(0,0),
                   labels = date_format("%H:%M", tz = "GMT")) +
  facet_wrap( ~ env.cond, nrow = 4, scales = "free_y", labeller = label_wrap_gen(multi_line=FALSE),
              strip.position = "right") +
  #theme_cowplot() +
  theme_minimal_hgrid(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "none") +
  ggtitle(paste0(strftime(min(dt$time), format = "%y-%m-%d %H:%M", tz = "GMT"),
                "-",
                strftime(max(dt$time), format = "%H:%M", tz = "GMT"),
                " Minepa (field - primary)"))


env_plot


ggsave(filename = paste0(dirname(solar_file), "/",
                         strftime(min(dt$time), format = "%y-%m-%d %H:%M", tz = "GMT"),
                         "-",
                         strftime(max(dt$time), format = "%H:%M", tz = "GMT"),
                         " Minepa (field - primary).png"),
       plot = env_plot,
       device = "png")



