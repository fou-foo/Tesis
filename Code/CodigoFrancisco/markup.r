
remove(list=ls())

# libraries
library(vars)
library(pls)
library(psych)

# path
dt.file <- "C:\\Users\\fou-f\\Documents\\GitHub\\MCE2\\4\\Tesina\\Code\\CodigoFrancisco\\"

# config model

# two models
h <- 13
lag.max <- 12

# var-pls model
runs <- 100
crit <- "FPE(n)"
season <- NULL
ec.det <- c("none", "const", "trend", "both")

# var model
seas.min <- 2
seas.max <- 12
lag.min <- 2
type.model <- c("none", "const", "trend", "boot")

# functions name
source(paste(dt.file, "funcion_pls_vars.r", sep = ""))

# load data
data <- read.csv(paste(dt.file, "inflacion.csv", sep = ""), row.names = 1)

# real exchange
data[,"p"] <- c(data[,"p"]/data[1,"p"])*100
data[,"peu"] <- c(data[,"peu"]/data[1,"peu"])*100
data[,"e"] <- c(data[,"e"]*data[,"peu"])/data[,"p"]

# markup model P=f(ER,W)
data <- log(as.matrix(data[,c("p","e","w")]))


# delete data for forecast
n <- nrow(data)-h
data.into <- data[1:n,]
data.fore <- data[c(n+1):c(n+h), ]
K <- ncol(data)

# p order
p <- VARselect(data.into, lag.max = lag.max, type = "const",
      season = season)$selection[crit]

# johansen test
summary(ca.jo(data.into, K = p, ecdet = "const", season = season))

list.fore.var.pls <- list.mape.var.pls <- list()

for(ncomp in 1 : c(K*p)){    print(ncomp)
  # VAR-PLS estimation
  var.pls1 <- var.pls(data.into, p = p, ncomp = ncomp, season = season)

  # forecast to h head with ci via bootstrap
  fore.var.pls <- ci.var.pls.boot(var.pls1, n.ahead = h, runs = runs)

  # mape stastitic
  mape.var.pls <- round(sapply(1:K, function(x) abs((fore.var.pls$Forecast[,x] -
                  data.fore[,x])/data.fore[,x]))*100, 2)
  colnames(mape.var.pls) <- colnames(data.into)

  # save results
  list.fore.var.pls[[ncomp]] <- fore.var.pls
  list.mape.var.pls[[ncomp]] <- mape.var.pls
}

# two models var: estimation with delete h observations in the sample and
# to make the forecast

# combination models
c.num.cols <- ncol(data.into)
idx.matrix  <- as.data.frame(matrix(0, (lag.max - 1)*(seas.max - 1)*
                length(ec.det), 3))
colnames(idx.matrix) <- c("seas", "lag", "ec.det")
idx.matrix[,"seas"] <- rep(2:seas.max, each = (length(ec.det))*(lag.max - 1))
idx.matrix[,"lag"] <- rep(c(2:lag.max), nrow(idx.matrix)/(lag.max - 1))
idx.matrix[,"ec.det"] <- rep(rep(ec.det, each = lag.max - 1), seas.max - 1)

# benchmark model ar(1)
ar1 <- sapply(1:ncol(data.into), function(x)
        predict(arima(data.into[1:c(nrow(data.into)-h), x],
        order = c(1,0,0)), n.ahead = h)$pred)
colnames(ar1) <- colnames(data.into)

# save data for test
test.data <- data.into[-c(1:c(nrow(data.into)-h)), ]

# matrix test
mat.test <- matrix(Inf, nrow(idx.matrix), 7)
colnames(mat.test) <- c("mape", "mdape", "rmspe", "rmdspe", "mrae",
                        "mdrae", "gmrae")

# runs all models
for(i in 1 : nrow(idx.matrix)){
  # specification models
  i.row <- idx.matrix[i, ]
  i.lag <- i.row[,"lag"]
  ec.det <- i.row[,"ec.det"]
  seas <- i.row[,"seas"]

  # season 2 is null
  if(seas == 2)
    seas <- NULL

  # test forecast 1
  var1.test <- VAR(data.into[1:c(nrow(data.into)-h), ], p = i.lag,
                type = ec.det, season = seas)
  fore.test <- predict(var1.test, n.ahead = h)

  # matrix forecast
  fore <- sapply(1:ncol(data.into), function(x) fore.test$fcst[[x]][,"fcst"])

  # error with model 1
  p.test <- sapply(1:ncol(data.into), function(x) c(test.data[,x]-fore[,x])/
              test.data[,x]*100)
  colnames(p.test) <- colnames(data.into)

  # statsitics for error (model 1)
  mape <- colMeans(abs(p.test))
  mdape <- apply(p.test, 2, function(x) median(abs(x)))
  rmspe <- sqrt(colMeans(p.test^2))
  rmdspe <- sqrt(apply(p.test, 2, function(x) median(x^2)))

  # relative errors for model 1
  r.test <- sapply(1:ncol(data.into), function(x) c(test.data[,x] -
              fore[,x])/c(test.data[,x] - ar1[,x]))

  # statistics for relative error (model 1)
  mrae <- colMeans(abs(r.test))
  mdrae <- apply(r.test, 2, function(x) median(abs(x)))
  gmrae <- apply(r.test, 2, function(x) geometric.mean(abs(x)))

  # statistics models
  mat.test[i,"mape"] <- mape[1]
  mat.test[i,"mdape"] <- mdape[1]
  mat.test[i,"rmspe"] <- rmspe[1]
  mat.test[i,"rmdspe"] <- rmdspe[1]
  mat.test[i,"mrae"] <- mrae[1]
  mat.test[i,"mdrae"] <- mdrae[1]
  mat.test[i,"gmrae"] <- gmrae[1]
}

# selected model
var.i <- apply(mat.test, 2, function(x) which(x == min(x)))

# save results
list.fore <- list.ci <- list()

# vars models selected
for(i in 1 : length(var.i)){

  # specification optimal
  seas <- idx.matrix[var.i[i],"seas"]
  lags <- idx.matrix[var.i[i],"lag"]
  type <- idx.matrix[var.i[i],"ec.det"]

  # var model
  var1 <- VAR(data.into, p = lags, type = type, season = seas)
  fore.var <- predict(var1, n.ahead = h)

  # fore
  list.fore[[i]] <- sapply(1:ncol(data.into),
                      function(x) fore.var$fcst[[x]][,"fcst"])
  # ci
  list.ci[[i]] <- sapply(1:ncol(data.into),
                      function(x) fore.var$fcst[[x]][,"CI"])
  colnames(list.fore[[i]]) <- colnames(list.ci[[i]]) <- colnames(data.into)
}
names(list.fore) <- names(list.ci) <- colnames(mat.test)

# integral forecast
fore.var <- sapply(1:length(var.i), function(x) list.fore[[x]][,"p"])
colnames(fore.var) <- colnames(mat.test)

# unique integral (mean & distribution quantile)
fore.int.var <- apply(fore.var, 1, function(x) quantile(x, probs = 0.5))

# mape with integral forecast
mape.int.var <- as.matrix(round(c(abs(c(data.fore[,"p"] -
              fore.int.var)/data.fore[,"p"])), 4)*100)
colnames(mape.int.var) <- "p"

# mape var models
mape.vars <- sapply(1:length(var.i), function(x)
              round(c(abs(c(data.fore[,"p"] - fore.var[,x])/
              data.fore[,"p"])), 4)*100)
colnames(mape.vars) <- colnames(mat.test)


mape.var.pls <- sapply(1:c(K*p), function(x) list.mape.var.pls[[x]][,"p"])

# models comparative
# all VAR-PLS vs integral VAR
comp1 <- sapply(1:c(K*p), function(x) mape.int.var > mape.var.pls[,x])
round(sum(comp1)/length(comp1), 4)*100

# all VAR selected vs all VAR-PLS
comp2 <- list()
for(i in 1 : length(var.i)){
  comp2[[i]] <- sapply(1:c(K*p),
    function(x) mape.vars[,i] > mape.var.pls[,x])
}
names(comp2) <- colnames(mat.test)

round(sapply(1:length(var.i), function(x) sum(comp2[[x]])/
  length(comp2[[x]])), 4)*100

# integral VAR-PLS
fore.var.pls <- sapply(1:c(K*p),
                function(x) list.fore.var.pls[[x]]$Forecast[,"p"])

# forecast
fore.int.var.pls <- apply(fore.var.pls, 1, function(x) quantile(x, probs = 0.5))

# mape with integral forecast
mape.int.var.pls <- round(abs(c(data.fore[,"p"] - fore.int.var.pls)/
                    data.fore[,"p"]), 4)*100

# integral VAR vs integral VAR-PLS
comp3 <- mape.int.var > mape.int.var.pls

# alls VAR cs integral VAR-PLS
comp4 <- mape.vars > mape.int.var.pls
round(sum(comp4)/length(comp4), 4)*100


mean(mape.int.var.pls)


mat <- cbind(data[,"p"], c(data.into[,"p"], fore.int.var.pls))


ts.plot(mat[-c(1:100),], col = c(1,2))

