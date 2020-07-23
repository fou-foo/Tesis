path <- "C:\\Users\\Sandybell Ferrer\\Desktop\\Tesis\\VincFacCienciasUNAM2020\\" # ubicacion del archivo 'Funciones_VARPLSParallel.R' y los datos
setwd(path)
dir()
data <- read.csv(file='Final_var_y_pls - Final_var_y_pls.csv', stringsAsFactors = FALSE)
data <- na.omit(data )
library(lubridate)
first <- unlist(strsplit(data$t[1], split = '/'))
first.date <- ymd( paste0(first[3], '-', first[2], '-01'))
data$time <- ymd(first.date) + months(0:11) # harcoding h params code 
library(ggplot2)
names(data)
# amarillo amazon FF9A00
ggplot(data , aes(x= time, y = Valor.Real))  + geom_line()+
  geom_line(aes(x= time, y= VAR.PLS, color=I('#005DA1'))) + 
  geom_line(aes(x= time, y= X2.50., color=I('#005DA1'))) + 
  geom_line(aes(x= time, y= X97.50., color=I('#005DA1'))) + 

  geom_line(aes(x= time, y= aws3, color=I('#FF9A00'))) + 
  geom_line(aes(x= time, y= awst, color=I('#FF9A00'))) + 
  geom_line(aes(x= time, y= aws98, color=I('#FF9A00'))) + 
  
  theme_minimal() + xlab('') + ylab('') + ylim(c(15.5, 23)) 
