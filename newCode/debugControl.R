

fun1<-function(model,data,control=list()) {
  print(control)
}

fun1(1,2,list(save.graphs=TRUE))

fun1(1,2,list(save.graphs=TRUE, nbiter.saemix=c(600,200)))


