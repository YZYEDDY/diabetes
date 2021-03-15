library(data.table)
library(parallel) # makeCluster detectCores stopCluster
library(doParallel) # registerDoParallel
library(lubridate)
library(sf)
library(raster)
library(rgdal)
library(survival)
library(splines)
library(ggplot2)
library(egg)
library(ggpubr)
library(DescrTab2)
library(dvmisc)

save.image(".RData")
#Tibbling-----------------------------------------------------------------------------------------
##Hospital Admission
admission <- readRDS("gzid.t2dm.rds") #.rds preserves data types and classes
##Air pollution and meteorology 
load("D:/Research Projects/Datasets/mc.rda")
##Define type of variables
admission[, id := as.character(id)][, sex := factor(sex)][, race := factor(race)][, doha := as.Date(doha)]
##Age.cut
adm.coord$age.cut <- quant_groups(adm.coord$age, 2)
##Season
ggplot(gzid.ap.tmp, aes(x = date, y = tmp)) +
  geom_point() +
  geom_smooth()
adm.coord[, season1 := "c"]#######################
adm.coord[month(doha) %in% 5:10, season1 := "w"]###season1
adm.coord[, season1 := factor(season1)]###########
##Geo-code
library(baidumap)
options(baidumap.key = 'GB7Gxf3ddT0NKE54Ch9HQfdODsZtQyeP')
a <- getCoordinate("广东省广州市番禺区市桥街道长堤东路十七巷5号", format = T)
address.coordinates[id == 661492,]

system.time({
  address <- unique(adm$address)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(address), .packages = c("baidumap", "data.table"), .combine = "rbind") %dopar%{
    options(baidumap.key = 'GB7Gxf3ddT0NKE54Ch9HQfdODsZtQyeP')
    coordinate <- getCoordinate(address[i], formatted = T)
    print(data.table(address = address[i], lon = coordinate[1], lat = coordinate[2]))
  }
  stopCluster(cl)
})
b <- merge(adm, a, all = T)
adm <- adm[order(adm$id, adm$doha),]
address.coordinates <- b

crs =  "+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
adm <- st_transform(st_as_sf(address.coordinates, coords = c('lon', 'lat'), crs = 4326), crs = crs)

#Exposure assessment---------------------------------------------------------------------------------
## Function: find control date
function(doha, lag = 0) {
  death.date.flag <- seq(floor_date(doha, "month"), ceiling_date(doha, "month") - 1, by = "days")
  return(death.date.flag[weekdays(death.date.flag) == weekdays(doha) & death.date.flag != doha] - lag)
}

##Assessment AP
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  gzid.ap <- foreach (i = 1:nrow(adm.sf), .packages = c("sf", "data.table", "lubridate", "raster", "rgdal"), .combine = "rbind") %dopar% {
    dates <- c(adm.sf$doha[i], find.control.dates(adm.sf$doha[i]))
    dates.lag <- data.table(sid = adm.sf$sid[i],
                            id = adm.sf$id[i],
                            date = c(dates, dates - 1, dates - 2, dates - 3), 
                            lag = rep(0:3, each = length(dates)), 
                            nday = rep(1:length(dates), 4))
    idw.points <- apc.study.area[apc.study.area$date %in% dates.lag$date, ]
    idw.points <- idw.points[lengths(st_intersects(idw.points, st_buffer(adm.sf[i, ], dist = 50000))) > 0, ]
    if (nrow(idw.points) > 0) {
      idw.dist <- data.table(data.table(idw.points)[, c(4:11)], dist = as.numeric(st_distance(adm.sf[i,], idw.points)))
      idw.dist[, weights := 1 / (dist ** 2)]
      idw.dist <- idw.dist[, lapply(.SD, function(x) sum(x * weights, na.rm = T) / sum(weights[!is.na(x)])), .SDcols = 2:8, keyby = date]
      idw.dist <- merge(dates.lag, idw.dist, by = "date", all.x = T)
      #idw.dist$tmp<- extract(subset(tmp, paste0("tmp.", format(idw.dist[, date], "%Y%m%d"))), 
              #adm.sf[i,], method = "bilinear")[1, ]
      return(idw.dist)
    }
  }
  stopCluster(cl)
})
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  gzid.ap.sup <- foreach (i = 1:nrow(adm.sf), .packages = c("sf", "data.table", "lubridate", "raster", "rgdal"), .combine = "rbind") %dopar% {
    dates <- c(adm.sf$doha[i], find.control.dates(adm.sf$doha[i]))
    dates.lag <- data.table(sid = adm.sf$sid[i],
                            id = adm.sf$id[i],
                            date = c(dates - 4, dates - 5, dates - 6, dates - 7), 
                            lag = rep(4:7, each = length(dates)), 
                            nday = rep(1:length(dates), 4))
    idw.points <- apc.study.area[apc.study.area$date %in% dates.lag$date, ]
    idw.points <- idw.points[lengths(st_intersects(idw.points, st_buffer(adm.sf[i, ], dist = 50000))) > 0, ]
    if (nrow(idw.points) > 0) {
      idw.dist <- data.table(data.table(idw.points)[, c(4:11)], dist = as.numeric(st_distance(adm.sf[i,], idw.points)))
      idw.dist[, weights := 1 / (dist ** 2)]
      idw.dist <- idw.dist[, lapply(.SD, function(x) sum(x * weights, na.rm = T) / sum(weights[!is.na(x)])), .SDcols = 2:8, keyby = date]
      idw.dist <- merge(dates.lag, idw.dist, by = "date", all.x = T)
      #idw.dist$tmp<- extract(subset(tmp, paste0("tmp.", format(idw.dist[, date], "%Y%m%d"))), 
      #adm.sf[i,], method = "bilinear")[1, ]
      return(idw.dist)
    }
  }
  stopCluster(cl)
})
##Assessment Meteorology
mc.tmp <- mc[["tmp"]]
tmp <- stack("tmp.tif")
writeRaster(mc.tmp, "tmp", format = 'GTiff', overwrite = T)
names(tmp) <- names(mc.tmp)

system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  gzid.tmp.sup <- foreach (i = 1:nrow(adm.sf), .packages = c("sf", "data.table", "lubridate", "raster", "rgdal"), .combine = "rbind") %dopar% {
    dates <- c(adm.sf$doha[i], find.control.dates(adm.sf$doha[i]))
    dates.lag <- data.table(sid = adm.sf$sid[i],
                            id = adm.sf$id[i],
                            date = c(dates - 4, dates - 5, dates - 6, dates - 7), 
                            lag = rep(4:7, each = length(dates)), 
                            nday = rep(1:length(dates), 4))
    id.tmp <- extract(subset(tmp, paste0("tmp.", format(dates.lag[, date], "%Y%m%d"))), 
              adm.sf[i,], method = "bilinear")[1, ]
    id.tmp <- data.table(dates.lag, tmp = id.tmp)
  }
  stopCluster(cl)
})

gzid.ap.tmp.sup <- merge(gzid.ap.sup, gzid.tmp.sup, all = T)
gzid.ap.tmp.sup$admission <- 0
gzid.ap.tmp.sup[nday == 1, admission := 1]
gzid.ap.tmp.sup[,tmp := tmp - 273.15]
gzid.ap.tmp.sup[, co := 1000*co]
gzid.ap.tmp.sup <- merge(gzid.ap.tmp.sup, adm.coord[, c("sid", "icd.code0", "sex", "age", "season1", "age.cut")], by = "sid")
gzid.ap.tmp.bind <- rbind(gzid.ap.tmp, gzid.ap.tmp.sup)
gzid.ap.tmp.bind <- gzid.ap.tmp.bind[order(sid, date),]
##Creat lag list--------------------------------------------------------------------------------------------------------------------
gzid.ap.tmp.bind[, lag := factor(lag, levels = c(0:7))]
summary(gzid.ap.tmp.bind$lag)
gzid.ap.tmp.bind[, age.cut := factor(age.cut, levels = c("[60,71]", "(71,103]"), labels = c("l", "h"))]

gzid.by.lag <- split(gzid.ap.tmp.bind, gzid.ap.tmp.bind$lag)

mean.ap <- function(lag, ap){rowMeans(sapply(gzid.by.lag[lag], `[[`, ap))}

gzid.by.lag[["01"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:2), simplify = F, USE.NAMES = T)))

gzid.by.lag[["02"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:3), simplify = F, USE.NAMES = T)))

gzid.by.lag[["03"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:4), simplify = F, USE.NAMES = T)))

gzid.by.lag[["04"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:5), simplify = F, USE.NAMES = T)))

gzid.by.lag[["05"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:6), simplify = F, USE.NAMES = T)))

gzid.by.lag[["06"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:7), simplify = F, USE.NAMES = T)))

gzid.by.lag[["07"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  do.call(data.table, sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp"), 
                                                             mean.ap, lag = c(1:8), simplify = F, USE.NAMES = T)))

gzid.by.lag[["03"]] <- data.table(gzid.by.lag[[1]][, c("sid", "id", "nday", "admission", "icd.code0", "sex", "age", "age.cut", "season1")], 
                                  pm2.5 = (gzid.by.lag[['0']]$pm2.5 + gzid.by.lag[['1']]$pm2.5 + gzid.by.lag[['2']]$pm2.5 + gzid.by.lag[['3']]$pm2.5)/4,
                                  pm10 = (gzid.by.lag[['0']]$pm10 + gzid.by.lag[['1']]$pm10 + gzid.by.lag[['2']]$pm10 + gzid.by.lag[['3']]$pm10)/4,
                                  so2 = (gzid.by.lag[['0']]$so2 + gzid.by.lag[['1']]$so2 + gzid.by.lag[['2']]$so2 + gzid.by.lag[['3']]$so2)/4,
                                  no2 = (gzid.by.lag[['0']]$no2 + gzid.by.lag[['1']]$no2 + gzid.by.lag[['2']]$no2 + gzid.by.lag[['3']]$no2)/4,
                                  co = (gzid.by.lag[['0']]$co + gzid.by.lag[['1']]$co + gzid.by.lag[['2']]$co + gzid.by.lag[['3']]$co)/4,
                                  o3.1h = (gzid.by.lag[['0']]$o3.1h + gzid.by.lag[['1']]$o3.1h + gzid.by.lag[['2']]$o3.1h + gzid.by.lag[['3']]$o3.1h)/4,
                                  o3.8h = (gzid.by.lag[['0']]$o3.8h + gzid.by.lag[['1']]$o3.8h + gzid.by.lag[['2']]$o3.8h + gzid.by.lag[['3']]$o3.8h)/4,
                                  tmp = (gzid.by.lag[['0']]$tmp + gzid.by.lag[['1']]$tmp + gzid.by.lag[['2']]$tmp + gzid.by.lag[['3']]$tmp)/4)

#Modeling-----------------------------------------------------------------------------------------------------
##Basic model

model <- clogit(admission ~ I(pm2.5/10) + ns(tmp, df = 3) + strata(sid), data = gzid.by.lag[["7"]]) %>% 
  summary()
exp.percentchange <- round(100*(model[["conf.int"]] - 1), 2)
p <- round(model[["coefficients"]][1,5], 4)
outcome <- paste0(exp.percentchange[1,1], " (", exp.percentchange[1, 3], ", ", exp.percentchange[1, 4], ")")
outcome <- data.table(outcome = outcome, p = p)
outcome

##Linear association
##Stratified: icd.code0 == "E11.900"; sex == "Male"; age.cut == "[60, 71]"; season1 == "c"
system.time({          ###Output to excel table
a.bind <- data.table() 
for (l in names(gzid.by.lag)) {
  a <- lapply(
    c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h"),
    function(i){
      model <- clogit(admission ~ I(get(i)/10) + ns(tmp, df = 3) + strata(sid), subset = season1 != "c", data = gzid.by.lag[[l]]) %>% 
        summary()
      exp.percentchange <- round(100*(model[["conf.int"]] - 1), 2)
      p <- round(model[["coefficients"]][1,5], 4)
      outcome <- paste0(exp.percentchange[1,1], " (", exp.percentchange[1, 3], ", ", exp.percentchange[1, 4], ")")
      outcome <- data.table(exposure = i, outcome = outcome, p = p)
      return(outcome)}
  )
  a <- rbindlist(a)
  a.bind <- cbind(a.bind, a[,c(2,3)])
}
})
write.csv(a.bind, "a.csv")

##non-linear association
dev.off()
par(mfrow = c(2,4))

a <- termplot(clogit(admission ~ ns(pm2.5, df = 3) + ns(tmp, df = 3) + strata(sid), data = gzid.by.lag[["0"]]), 
         terms = "ns(pm2.5, df = 3)", se = T, plot = F) #Single plot
a <- a[[1]]
ggplot(a, aes(x)) +
  geom_line(aes(y = exp(y))) +
  geom_ribbon(aes(ymin = exp(y - 1.96*se), ymax = exp(y + 1.96*se)), fill = "grey70", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")

###By lag, just change the dataset
termplot.lag <- function (ap) {
  term <- termplot(clogit(as.formula(paste0("admission ~ ns(", ap, ", df = 3) + ns(tmp, df = 3) + strata(sid)")), data = gzid.by.lag[["01"]]), 
         terms = paste0("ns(", ap, ", df = 3)"), se = T, plot = F)
  term <- data.table(term[[1]])
  setnames(term, "x", ap)
  return(term)}
term.list <- sapply(c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h"), termplot.lag, simplify = F, USE.NAMES = T) #Save parameters of all regressions
plot.list <- sapply(term.list, function(i) {
  ggplot(i, aes(x = get(names(i)[1]))) +
    geom_line(aes(y = exp(y))) +
    geom_ribbon(aes(ymin = exp(y - 1.96*se), ymax = exp(y + 1.96*se)), fill = "grey70", alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = paste0(toupper(names(i)[1]), ", ", expression(μg/m^3)), y = "Percent change")
}, simplify = F, USE.NAMES = T)
do.call(ggarrange, plot.list)

#Spearman correlation-----------------------------------------------------------
library(Hmisc)
a <- rcorr(as.matrix(gzid.ap.tmp[, c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h", "tmp")]), type = "spearman")
a <- round(a$r, 2)
write.csv(a, "a.csv")

#Plot for linear relationship across lag days
##Output for plots
system.time({           
  a.bind <- data.table()
  for (l in names(gzid.by.lag)) {
    a <- lapply(
      c("pm2.5", "pm10", "so2", "no2", "co", "o3.1h", "o3.8h"),
      function(i){
        model <- clogit(admission ~ I(get(i)/10) + ns(tmp, df = 3) + strata(sid), data = gzid.by.lag[[l]]) %>% 
          summary()
        exp.percentchange <- data.table(round(100*(model[["conf.int"]] - 1), 4))[1, c(1,3,4)]
        #p <- round(model[["coefficients"]][1,5], 4)
        #outcome <- paste0(exp.percentchange[1,1], " (", exp.percentchange[1, 3], ", ", exp.percentchange[1, 4], ")")
        #outcome <- data.table(exposure = i, outcome = outcome, p = p)
        outcome <- data.table(ap = i, lag = l, exp.percentchange)
        return(outcome)}
    )
    a <- rbindlist(a)
    a.bind <- rbind(a.bind, a)
    colnames(a.bind) <- c("ap", "lag", "estimate", "lower", "upper")
    a.bind[, lag := factor(lag, levels = names(gzid.by.lag))][, ap := factor(ap)]
  }
})
##Plot
a.bind.list <- split(a.bind, a.bind$ap)
plot.list <- lapply(a.bind.list, function(i){ggplot(i, aes(x = lag, y = estimate)) +
  geom_point() +
  labs(x = NULL, y = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper))})
figure <- ggarrange(plotlist = plot.list, ncol = 1, nrow = 7)
annotate_figure(figure,
                left = text_grob(expression(paste("Percent change in diabetes hospitalizations per 10", mu, "g/m"^"3", "increase of air pollutants")), 
                                 rot = 90, size = 8, face = "bold"),
                bottom = text_grob("Lags", size = 10, face = "bold", vjust = -0.6))
