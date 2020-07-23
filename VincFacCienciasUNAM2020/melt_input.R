input.csv <- read.csv("C:/Users/usuario/Downloads/input.csv.txt")
library(reshape2)
names(input.csv)
first <- input.csv$Date[1]
first.date <- unlist(strsplit(first, split = '/'))
first.date.string <- paste0(first.date[[3]], '-', first.date[[2]], '-01')
input.csv <- head(input.csv, dim(input.csv)[1] - 12) # el 12 es por el horizonte a comparar con el codigo de VAR-PLS 
input.csv$t <- ymd(first.date.string) + months(0: (dim(input.csv)[1] -1) )
input.csv$Date <- NULL
output.data <- melt(input.csv, id = 't')
colnames(output.data) <- c('timestamp ', 'item_id', 'demand')
output.data <- output.data[, c('item_id', 'timestamp ', 'demand'   )]
class(output.data$`timestamp ` )
library(lubridate)
write.csv(output.data, 'C:/Users/usuario/Downloads/ouput_melt_justo.csv', row.names = FALSE)
