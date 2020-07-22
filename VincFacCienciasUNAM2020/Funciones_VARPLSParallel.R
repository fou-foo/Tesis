# Closure para poder paralelizar la funcion 'CarloMagno'
Rellena <- function(dataframe, dataframe.pares, index, crit, lag.max, significancia)
{
  # Entradas:
    # dataframe (data.frame): Conjunto de datos original
    # dataframe.pares (data.frame): Tabla con todas las posibles combinaciones por pares de variables
    # index (int): variable anonima
    # crit (character): Sttring con el criterio para determinar el orden del VAR(p)
    # lag.max (int): Parametro para determinacion del orden del VAR(p)
    # significancia (string): Significancia de la prueba '1pct', '5pct', '10pct'

  # Salidas: Funcion que determina si dos variables son cointegradas
  data <- dataframe
  Cointegracion.pares <- dataframe.pares
  crit <- crit
  lag.max <- lag.max
  significancia <- significancia
  function(index)
  {
    series <- Cointegracion.pares[index, names(Cointegracion.pares)[1:2] ]
    series <- unlist(series)
    datos <- data[, series ]
    p.mini <- VARselect(y= datos, lag.max = lag.max, type = "const")$selection[crit]
    l <- ca.jo(datos, K= p.mini)
    rango1 <- l@cval[ 1  , significancia] # valor critico para rango <=1 El numero magico es porque siempre hay solo dos series
    estadistico1 <- l@teststat[1]
    resultado <- ifelse( estadistico1 < rango1, 'Si', 'No')
    return(resultado)
  }
}

# Funcion para determinar si existe cointegracion en un VAR con mas de once variables, basado en el paper de
# Discovering common trends in a large set of disaggregates: statistical procedures and their properties
# En el mismo paper se da otra propuesta para el caso de cientos
# requiere de las librerias 'vars' y parallel
CarloMagno <- function(data, significancia='1pct', crit = 'FPE(n)', lag.max)
{
  # Entradas:
  # data (data.frame): Conjunto de datos original
  # significancia (string): Significancia de la prueba '1pct', '5pct', '10pct'
  # crit (string): criterio para seleccion del orden del VAR(p) FIJO DE MOMENTO
  # lag.max (int): Limite para los lags de todos los VAR(p) se se estiman
  # Salidas:
    # cointegracion.grande (character): Mensaje indicando si hay o no cointegracion
  K <- dim(data)[2]
  # La primer etapa consite en encontrar cointegracion por pares
  Cointegracion.pares <- as.data.frame( t(combn(K, 2))) #las combinaciones por pares de variables
  Cointegracion.pares$Cointegran <- ''
  Rellena.cointegracion <- Rellena(dataframe=data, dataframe.pares=Cointegracion.pares,  crit=crit, lag.max=lag.max, significancia = significancia )
  Cointegracion.pares$Cointegran <- mapply(function(x){Rellena.cointegracion(x)}, 1:dim(Cointegracion.pares)[1])
  if( !('No' %in% Cointegracion.pares$Cointegran ) )
  {
    cointegracion.grande <- 'Todas las variables cointegran'
    print(cointegracion.grande)
    return(cointegracion.grande)
  } else{
    cointegracion.grande <- 'No todas las variables cointegran'
    print(cointegracion.grande)
    return(cointegracion.grande)
  }
}

# Funcion para generar las matrices con todos los lags
Span <- function(x, lag=NULL)
{
    # Entradas:
        #  x (matrix): matriz con variables originales
        #  lag (int): Lag a calcular
    # Salidas:
        # m (matrix): Matriz con las mismas dimensiones que 'x' pero con el 'lag'
    m <- sapply(x, function(x) dplyr::lead(x, n = lag))
    return(m)
}
# Vectorizacion de la funcion anterior
SpanMatrix <- function(x, p=NULL)
{
    # Entradas:
    #  x (matrix): matriz con variables originales
    #  p (int): Orden del VAR(p)
    # Salidas:
    # m (matrix): Matriz con todos los lags de dimensiones (n, p*k)
    nombres.variables <- colnames(x)
    M <- lapply(X=0:p, function(i) Span(x, lag=i   ) )
    M <- do.call('cbind', M)
    colnames(M) <- paste0(rep(nombres.variables, p+1),'.', rep(0:p, each=length(nombres.variables)))
    return(M)
}
# Funcion para pronosticar h pasos hacia adelante
Predict.PLS <- function(modelo, original, ncomp, h)
{
    # Entradas
        #  modelo (pls): Modelo PLS entrenado
        #  original (data.frame): Conjunto de Observaciones a pronosticar
        #  ncomp (int): Numero de componentes con las que se desea pronosticar
        #  h (int): Numero de pronosticos a realizar
    # Salidas
        #  salidas (vector) : Los h pronosticos para cada serie
    Y.lags <- modelo$y
    X.lags <-modelo$x
    Y <- original
    salidas <- matrix(nrow = h, ncol=dim(Y.lags)[2])
    entradas <- matrix(nrow = dim(X.lags)[1], ncol=dim(X.lags)[2])
    colnames(salidas) <- colnames(Y.lags)
    # primer pronostico
    pronosticos <- predict( modelo, newdata=head(X.lags,1) )[,,ncomp]# como lo recomienda Frances en el paper
    salidas[1, ] <- pronosticos
    X.aux <- SpanMatrix(Y, p=p)
    for (i in 2:h)
    {
        # se pronostica iteritivamente los demas pronosticos, en cada iteracion se reajusta el modelo
        Y.h <- rbind(head(salidas,i-1), head(Y.lags, nrow(Y.lags) - i+1))
        X.h <- rbind( head(X.aux, i-1), head(X.lags, nrow(X.lags) - i+1)  )
        model <- plsr( Y.h ~ X.h, ncomp= ncomp, method='simpls') #
        pronosticos <- predict(model, newdata = head(X.h, 1 ))[,,ncomp]
        salidas[i, ] <- pronosticos
    }
    return(salidas)
}

# Funcion para calcular el error relativo porcentual del pronostico
Error.relativo <- function(x, ncomp, data, h )
{
    # Entradas :
        #  x (data.frame): Tabla con los h pronosticos de todas las variables de Y.lags
        #  ncomp (int): Numero de componentes con el quue se pronostico
        #  data (data.frame): Valores originales de las series
        #  h (int): Numero de pronosticos que se realizaron
    # Salidas:
        #  Error.relativo.medio (data.frame): Tabla con un solo registro con el Error relativo de la variable a pronosticar
                                  #           para cada uno de los h orizontes
    Res1 <- x
    index.InflacionNacional <- grep('InflacionNacional', colnames(Res1))
    index.Salarios<- grep('Salarios', colnames(Res1))
    index.Precios.EUA <- grep('Precios.EUA', colnames(Res1))
    index.Tipo.de.cambio <- grep('Tipo.de.cambio', colnames(Res1))
    index.Industria.EUA <- grep('Industria.EUA', colnames(Res1))
    index.Actividades.Primarioas <- grep('Actividades.Primarioas', colnames(Res1))
    index.Actividades.Secundarias<- grep('Actividades.Secundarias', colnames(Res1))
    index.Actividades.Terciarias<- grep('Actividades.Terciarias', colnames(Res1))
    index.Mineria <- grep('Mineria', colnames(Res1))
    index.Manufactura <- grep('Manufactura', colnames(Res1))
    index.Construccion <- grep('Construccion', colnames(Res1))
    index.Energia<- grep('Energia', colnames(Res1))
    index.Desempleo<- grep('Desempleo', colnames(Res1))
    index.Billetes.y.monedas<- grep('Billetes.y.monedas', colnames(Res1))
    index.Ingreso.Mexico<- grep('Ingreso.Mexico', colnames(Res1))
    index.TIIE28<- grep('TIIE28', colnames(Res1))
    index.TIIE91<- grep('TIIE91', colnames(Res1))
    reorden <- c(index.InflacionNacional, index.Salarios, index.Precios.EUA,
             index.Tipo.de.cambio, index.Industria.EUA, index.Actividades.Primarioas,
             index.Actividades.Secundarias, index.Actividades.Terciarias,
             index.Mineria, index.Manufactura, index.Construccion,
             index.Energia, index.Desempleo, index.Billetes.y.monedas,
             index.Ingreso.Mexico, index.TIIE28, index.TIIE91)
    Res1 <- Res1[, reorden]
    marcas <- seq(1, dim(Res1)[2], by=h)
    Res1 <- Res1[h:1, marcas]
    verdadero <-  tail(data, h)[h:1,]
    Error.relativo.medio <- abs(Res1 - verdadero)/abs(verdadero)*100
    Error.relativo.medio <- as.data.frame(t(apply(Error.relativo.medio, 2, mean)))
    row.names(Error.relativo.medio) <- ncomp
    return(Error.relativo.medio)
}

# Funcion para calcular el mape relativo porcentual del pronostico
MAPE <- function(x, ncomp, data, h )
{
    # Entradas :
    #  x (data.frame): Tabla con los h pronosticos de todas las variables de Y.lags
    #  ncomp (int): Numero de componentes con el quue se pronostico
    #  data (data.frame): Valores originales de las series
    #  h (int): Numero de pronosticos que se realizaron
    # Salidas:
    #  mape (data.frame): Tabla con un solo registro con el mape de la variable a pronosticar
    #           para cada uno de los h orizontes
    Res1 <- x
    index.InflacionNacional <- grep('InflacionNacional', colnames(Res1))
    index.Salarios<- grep('Salarios', colnames(Res1))
    index.Precios.EUA <- grep('Precios.EUA', colnames(Res1))
    index.Tipo.de.cambio <- grep('Tipo.de.cambio', colnames(Res1))
    index.Industria.EUA <- grep('Industria.EUA', colnames(Res1))
    index.Actividades.Primarioas <- grep('Actividades.Primarioas', colnames(Res1))
    index.Actividades.Secundarias<- grep('Actividades.Secundarias', colnames(Res1))
    index.Actividades.Terciarias<- grep('Actividades.Terciarias', colnames(Res1))
    index.Mineria <- grep('Mineria', colnames(Res1))
    index.Manufactura <- grep('Manufactura', colnames(Res1))
    index.Construccion <- grep('Construccion', colnames(Res1))
    index.Energia<- grep('Energia', colnames(Res1))
    index.Desempleo<- grep('Desempleo', colnames(Res1))
    index.Billetes.y.monedas<- grep('Billetes.y.monedas', colnames(Res1))
    index.Ingreso.Mexico<- grep('Ingreso.Mexico', colnames(Res1))
    index.TIIE28<- grep('TIIE28', colnames(Res1))
    index.TIIE91<- grep('TIIE91', colnames(Res1))
    reorden <- c(index.InflacionNacional, index.Salarios, index.Precios.EUA,
                 index.Tipo.de.cambio, index.Industria.EUA, index.Actividades.Primarioas,
                 index.Actividades.Secundarias, index.Actividades.Terciarias,
                 index.Mineria, index.Manufactura, index.Construccion,
                 index.Energia, index.Desempleo, index.Billetes.y.monedas,
                 index.Ingreso.Mexico, index.TIIE28, index.TIIE91)
    Res1 <- Res1[, reorden]
    marcas <- seq(1, dim(Res1)[2], by=h)
    Res1 <- Res1[h:1, marcas]
    verdadero <-  tail(data, h)[h:1,]
    mape <- abs(Res1 - verdadero)/abs(verdadero)
    mape <- as.data.frame(t(apply(mape, 2, function(x) sum(x)*(100/h)  )))
    row.names(mape) <- ncomp
    return(mape)
}

# Funcion para efectuar una sola replica de Bootstrap
Bootstrap <- function(x, X.lags, p, Y, Y.lags, Ncomp.mape)
{
    # Entradas:
        #   x (matrix): Matriz con los datos originales
        #   X.lags (matrix): Matriz con todos los lags
        #   p (int): Orden del VAR(p)
        #   Y  (matrix): Matriz con las observaciones a pronosticar en 'Y'
        #   Y.lags (matrix): Matriz con los h lags de 'Y'
        #   Ncomp.mape (int): Numero de componentes a considear en los pronosticos
    # Salidas:
        #   Res1 (vector): Vector de longitud h, la replica bootstrap del pronostico
    x <- X # copia de la matriz 'x' con las variables originales sin lags
    X.lags <- X.lags # copia de la matriz con todos los lags
    p <- p  # copia del orden del AR(p)
    Y <- Y
    Y.lags <- Y.lags
    # Matrices y elementos
    K <- ncol(x)
    total <- nrow(x)
    modelo <- plsr(Y.lags~ X.lags, method = 'simpls', ncomp=Ncomp.mape)
    # Coeficientes de modelo PLS
    # Paso 1 del procedimiento bootstrap
    # estimar B y fijar los primeros 'p' elementos
    B <- coef(modelo)[,,1]
    residuos <- scale(modelo$resid[ , , Ncomp.mape], scale = FALSE)*( (total-p)/ (total-2*p))**(.5) # guarda los rresiduos y los centra y scale por el factor que dice el paper
    iniciales <- tail(as.data.frame(na.omit(X.lags)), p)
    Y.muestra.boot <- matrix( nrow = total, ncol = dim(Y.lags)[2]) # aqui guarda todas las muestras de bootstrap
    X.muestra.boot <- matrix( nrow = total, ncol = dim(X.lags)[2]) # aqui guarda todas las muestras de bootstrap
    colnames(Y.muestra.boot) <- colnames(Y.lags)
    colnames(X.muestra.boot) <- colnames(X.lags)
    index.res <- sample(1:dim(residuos)[1], p, replace = TRUE)
    Y.muestra.boot[1:p, ] <- as.matrix(iniciales) %*% B + residuos[index.res, ]
    X.muestra.boot[1:p, ] <- as.matrix(iniciales)
    # paso 2
    # generar una muestra de bootstrap de longitud 't'
    for(i in 2:(total-p+1))
    {
        index.res <- sample(1:dim(residuos)[1], p, replace = TRUE)
        Y.muestra.boot[i:(i+p-1), ] <- X.muestra.boot[(i-1):(i+p-2), ] %*% B + residuos[index.res, ]
        X.muestra.boot[i:(i+p-1), ] <- X.lags[(i-1):(i+p-2), ]
    }
    # paso 3
    # Restimacion de coeficientes B y prediccion
    model <- plsr(Y.muestra.boot ~ X.muestra.boot , method = 'simpls', ncomp=Ncomp.mape, x=TRUE, y=TRUE)
    Res1 <- Predict.PLS(modelo=model, original=Y, ncomp = Ncomp.mape, h=h)
    Res1 <- Res1[, 4]
    return(Res1)
}
