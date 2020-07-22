############################################
# Implementacion de la metodologia VAR-PLS #
# J. Antonio Garcia jose.ramirez@cimat.mx  #
############################################

remove(list=ls()) # removemos todos los objetos del enviroment
{
###########################################
# librerias                               #
{
t1 <- Sys.time()
library(vars) # tiene varias dependencias (implicitamente carga en el ambiente otras) util para modelos VAR
library(pls)  # para estimacion pls
library(psych) # solo ocupamos una funcion de aqui CHECAR CUAL ES
library(ggplot2) # libreria de graficos
library(lubridate)  # libreria para manejo de fechas
library(reshape2) #manipulacion de dataframes
library(parallel)
}
###########################################
# Parametros                              #
{
path <- "C:\\Users\\Sandybell Ferrer\\Desktop\\Tesis\\VincFacCienciasUNAM2020\\" # ubicacion del archivo 'Funciones_VARPLSParallel.R' y los datos
h <- 16 # numero de steps a pronostricar
lag.max <- 6 # lag maximo para la determinacion inicial del AR(p)
runs <- 100  # numero de iteraciones bootstrap para los intervalos de confianza
crit <- "FPE(n)" # criterio con cual elegir el orden inicial del VAR(p)
confianza <- .95
}
source(paste0(path, "Funciones_VARPLSParallel.R"))# cargar funciones auxiliares
##########################################
# Lectura de datos                       #
# Se espera un dataframe donde la primer columna sean las fechas de las series y la segunda la variable de interes a pronosticar
data <- read.csv(paste0(path, "Compendio07abril2020Dolar.csv"), row.names = 1)
##########################################
# imputacion de datos                    #
data <- na.omit(data)
temp <- row.names(data)
data <- as.data.frame(sapply(data, log))
row.names(data) <- temp
data2 <- data
write.csv(data2, file=paste0(path, 'input.csv') )
#data[, 3:17] <- NULL
##########################################
# visualizacion
data.s <- data
data.s$time <- row.names(data)
data.s$time <- dmy(data.s$time)
data.s <- melt(data.s, id='time')
ggplot(data.s, aes(x= time, y =exp(value), color=variable)) + geom_line() +
  facet_wrap(variable~., scales = "free") +  theme_minimal() + xlab('') +
  ylab('') +  theme(legend.position = "bottom", legend.title = element_text(color = "white")) +
  ggtitle('Variables macroeconómicas') + guides( color=FALSE)
ggplot( subset(data.s, time >= dmy('01/01/2010')),
        aes(x= time, y =exp(value), color=variable)) + geom_line() +
  facet_wrap(variable~., scales = "free") +  theme_minimal() + xlab('') +
  ylab('') +  theme(legend.position = "bottom", legend.title = element_text(color = "white")) +
  ggtitle('Variables macroeconómicas') + guides( color=FALSE)
###########################################
fechas <- row.names(data) # guardamos las fechas
row.names(data) <- fechas
# division de la muestra como sugiere Frances
n <- dim(data)[1] # tamaÃ±o total de la muestra
k <- dim(data)[2] # numero de componentes
Y <- tail(data, n - h + 1 )
X <- head(data, n - h + 1 )
# test de cointegracion de





CarloMagno(data, lag.max=lag.max)
p <- VARselect(y= X, lag.max = lag.max)$selection[crit]# determinacion del orden del VAR(p)
Y <- Y[ dim(Y)[1]:1, ] # Frances propone este acomodo de las observaciones
X <- X[ dim(X)[1]:1, ]
Y <- as.data.frame(Y)
X <- as.data.frame(X)
colnames(Y) <- colnames(X) <- colnames(data)
X.lags <- SpanMatrix(X, p=(p))# generamos la matriz extendida con todos los lags
Y.lags <- SpanMatrix(X, p=(h-1))# generamos la matriz extendida con todos los lags
model <- plsr(Y.lags~ X.lags, method = "simpls", x=TRUE, y=TRUE)
componentes.practicas <- model$ncomp # por los nulos se reducen
t1 <- Sys.time()
Pronosticos <- mclapply(FUN=function(x) Predict.PLS(modelo=model, original=Y, ncomp=x, h=h),                  1:componentes.practicas)
t2 <- Sys.time()
t2 -t1


temp1 <- temp <- matrix(nrow = componentes.practicas, ncol = dim(data)[2])
temp1 <- temp <- as.data.frame(temp)
colnames(temp1) <- colnames(temp) <- colnames(data)
row.names(temp1) <- row.names(temp) <- paste0('Comp.', 1:componentes.practicas)
for(i in 1:componentes.practicas)
{
    temp[i, ] <- Error.relativo(Pronosticos[[i]], ncomp = i, data, h)
    temp1[i, ] <- MAPE(x=Pronosticos[[i]], ncomp=i, data, h)
}
temp1$id <- 1:componentes.practicas
Ncomp <- which.min(temp[, 4])
t2 <- Sys.time()
}
# grafica de mape
z <- melt(temp1, id='id')
Ncomp.mape <- which.min(temp1[, 4])
ggplot(z, aes(x=id, y=value, color=variable)) + geom_line() +
    facet_wrap(variable~., scales = "free") +  theme_minimal() + xlab('') +
    ylab('') +  theme(legend.position = "bottom", legend.title = element_text(color = "white")) +
    ggtitle('MAPE') +guides( color=FALSE)
###################
###################
###################
# intervalos de confianza
set.seed(0)
runs <- 100
t2 <- Sys.time()
Intervalos <- mclapply(rep(p, runs), function(x) Bootstrap(x=X, X.lags = X.lags, p=x, Y=Y, Y.lags = Y.lags, Ncomp.mape=Ncomp.mape)  )
t3 <- Sys.time()
t3 - t2

intervalos <- do.call('rbind', Intervalos)
intervalos <- as.data.frame(intervalos)
significancia <- (1 - confianza)/2
c.i <- lapply(intervalos, function(x)quantile(x, probs=c(significancia, 1-significancia )) )
c.i <- as.data.frame(do.call('rbind', c.i))
c.i$l <- c.i[, 2] -c.i[, 1]
write.csv(c.i, paste0(path, 'resultados_var_pls.csv') )
#######################################################
### estimacion del VAR a comparar
var <- VAR(ts(data2, start = c(2015, 1), frequency = 12), p=p, ic='FPE',
           lag.max = lag.max )
tabla.var <- data.frame(Pronostico=predict(var, n.ahead = h)$fcst$Tipo.de.cambio[,'fcst'],
                        y = tail(data2[,4],h ))
tabla.var$Error.relativo <- abs(tabla.var$y - tabla.var$Pronostico)/abs(tabla.var$Pronostico)*100
tabla.var
mean(tabla.var$Error.relativo)
######################################
# plot de pronostico conjunto
Final <- data.frame(Valor.Real=tail(data2[,4], h) ,
                    VAR.PLS = Predict.PLS(modelo=model, original=Y, h = h, ncomp=Ncomp.mape)[,4],
                    VAR=predict(var, n.ahead = h)$fcst$Tipo.de.cambio[,'fcst'])
tabla <- Final <- exp(Final)
Final <- cbind(Final, c.i[, 1:2])
Final[ , c(4,5)] <- Final[ , c(4,5)] - 2
Final[, 4] <- Final$VAR.PLS - Final[, 4]
Final[, 5] <- Final$VAR.PLS + Final[, 5]
Final$t <- dmy(row.names(tail(data, h)))
write.csv(Final, file=paste0(path, 'Final_var_y_pls.csv'))
data <- exp(data2)
data$t <- dmy(row.names(data))
t <- melt(Final, id='t')
tabla$Error.relativo.pls <- abs(tabla$Valor.Real - tabla$VAR.PLS )/abs(tabla$Valor.Real)*100
tabla$Error.relativo.var <- abs(tabla$Valor.Real - tabla$VAR )/abs(tabla$Valor.Real)*100
tabla
tabla$Presicio.VAR.PLS <- 100 - tabla$Error.relativo.pls
tabla$Presicio.VAR <- 100 - tabla$Error.relativo.var
library(xtable)
xtable(tabla)
sapply(tabla, mean)
xtable(Final[, c(2, 4,5)])
write.csv(tabla, file='error.csv')
ggplot(data = data, aes(y=Tipo.de.cambio, x=t)) + geom_line() +
    geom_line(data=subset(t, variable != 'VAR'),
              aes(x=t, y=value, color=variable)) + theme_minimal() +
    ylab('') + xlab('') +xlim(ymd('2016-06-01'), max(data$t)) +
    ylim(c(16,22)) + guides( color=FALSE) +
    ggtitle(paste0('Pronóstico VAR(', p, ')-PLS(h=', h, ', k=', Ncomp.mape, ')'))  +
    scale_color_manual(values = c('Valor.Real'= 'navy', 'VAR.PLS' = 'red',
                                  '2.5%' = 'orange', '97.5%' = 'orange'))
t4 <- Sys.time()
tabla
