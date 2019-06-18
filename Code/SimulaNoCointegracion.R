l <- 5
n <- 166
set.seed(0)
y <- matrix(nrow = n, ncol=l)
y[1, ] <- runif(l, -2, 2)
y[2, ] <- runif(l, -2, 2)
a1 <- runif(l*l, -.5, .5)
dim(a1) <- c(l,l)
a2 <- runif(l*l, -.5, .5 )
dim(a2) <- c(l,  l)
# ar
for(i in 3:n)
{
    y[i, ] <- a1%*%y[i-1,] + a2%*%y[i-2, ] + rnorm(l)
}
library(vars)
plot(ts((y)))
p <- VARselect(y)
p <- p$selection[4]
colnames(y) <- paste0('V', 1:(l))
a <-ca.jo(y)
summary(a)
simulacion <- as.data.frame(y)

a <- CarloMagno(simulacion, lag.max = 10)
