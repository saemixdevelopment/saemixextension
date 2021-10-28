myfunc<-function(x, firstmiss, secondmiss) {
  if(missing(firstmiss)) firstmiss<-0
  if(missing(secondmiss)) secondmiss<-0
  sum(x)+firstmiss+secondmiss
}

myfunc(1:10)
myfunc(1:10,1)
myfunc(1:10,1,2)
myfunc(1:10,secondmiss=2)
myfunc(1:10,secondmiss=2,firstmiss=1)
