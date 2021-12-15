# Using autonpde within a function

testfun <- function(namobs, namsim)
  npdeObject<-autonpde(namobs=namobs, namsim=namsim, detect=TRUE)

# Works
datDir<-"/home/eco/work/npde/npde30/npde/data"
y1 <- testfun(file.path(datDir, "theopp.tab"), file.path(datDir, "simtheopp.tab"))

# Doesn't work
obsdat<-read.table(file.path(datDir, "theopp.tab"), header=T)
y1 <- try(testfun(obsdat, file.path(datDir, "simtheopp.tab")))

# Works
y2 <- autonpde(obsdat, namsim= file.path(datDir, "simtheopp.tab"), detect=TRUE)

# Doesn't work
# writing to a temporary file
testfun2 <- function(namobs, namsim) {
  if(is(namobs,"data.frame")) {
    namfile<-tempfile()
    write.csv(namobs, namfile)
    cat(namfile,"\n")
    npdeObject<-autonpde(namobs=namfile, namsim=namsim, detect=TRUE)
  } else
  npdeObject<-autonpde(namobs=namobs, namsim=namsim, detect=TRUE)
}

y1 <- try(testfun2(namobs=obsdat, namsim=file.path(datDir, "simtheopp.tab")))

############################################################################
# requirements 
# function accepting both dataframes and 
# can be used both interactively or within another function (tried a double encapsulation, still works)

encapsfun <- function(namobs, namsim) {
  basefun(namobs=namobs, namsim=namsim)
}

encapsfun2 <- function(namobs, namsim) {
  encapsfun(namobs=namobs, namsim=namsim)
}

basefun <- function(namobs, namsim) {
  # req: namobs, namsim can be either dataframes or files on disk
  if(is(namobs, "data.frame")) dat1<-namobs else 
    dat1<-read.table(namobs, header = T)
  if(is(namsim, "data.frame")) dat2<-namsim else 
    dat2<-read.table(namsim, header = T)
  print(summary(dat1))
  print(summary(dat2))
}

# Testing
data1<-data.frame(id=1:10, x=runif(10))
data2<-data.frame(time=1:10, x=rnorm(10))
basefun(data1, data2)
encapsfun(data1,data2)

namfile<-tempfile()
write.table(data1, namfile, row.names=F)

basefun(namfile, data2)
encapsfun(namfile,data2)

# Works
for(i in 1:3) {
  data1<-data.frame(id=1:10, x=runif(10))
  data2<-data.frame(time=1:10, x=rnorm(10))
  encapsfun2(data1,data2)
  
  namfile<-tempfile()
  write.table(data1, namfile, row.names=F)
  encapsfun2(namfile,data2)
}

# So... try to recode autonpde to read the data and feed it to npdeData ?



############################################################################