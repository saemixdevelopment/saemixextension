library(tibble)

df1 <- tibble(x = runif(10), y = x * 2)

df2<-as.data.frame(df1)
