## Current definitions: two different generics
if(FALSE) {
  
# npde
setGeneric(name="read",
           def=function(object,name.data,header=TRUE,sep="",na.strings=c("NA","."),detect=TRUE,verbose=FALSE){standardGeneric("read")})

# saemix
setGeneric(name="read",
             def=function(object, dat=NULL){standardGeneric("read")}
  )

# read.delim, read.csv,...
# file, header = TRUE, sep = "\t", quote = "\"", dec = ",",  fill = TRUE, comment.char = "", ...
# but read.table has different definition

####################### Common signature, mini-classes
datDir<-"/home/eco/work/saemix/saemixextension/data"


# Mini-classes
setClass(
  Class="MiniSaemixData",
  representation=representation(
    name.data="character",	# name of dataset
    data="data.frame"),
  validity=function(object){
    return(object)
  }
)

setClass(
  Class="MiniNpdeData",
  representation=representation(
    name.data="character",	# name of dataset
    data="data.frame"),
  validity=function(object){
    return(object)
  }
)


# Generic
setGeneric(name="read",
           def=function(object,name.data,header=TRUE,sep="",na.strings=c("NA","."),detect=TRUE,verbose=FALSE){standardGeneric("read")}
)
}


# Methods
setMethod("myread",
          signature="MiniSaemixData",
          function(object, dat = NULL) {
          }
)

setMethod("myread",
          signature="MiniNpdeData",
          function(object, dat = NULL) {
          }
)

setMethod("myread",
          signature=c("MiniNpdeData","data.frame"),
          function(object, dat = NULL) {
          }
)

####################### Common signature for full classes
npdeDir<-"/home/eco/work/npde/npde30/npde/R"
saemixDir<-"/home/eco/work/saemix/saemixextension/R"
source(file.path(npdeDir, "NpdeData.R"))
source(file.path(saemixDir, "SaemixData.R"))

# Methods
setMethod("myread",
          signature="SaemixData",
          function(object, dat = NULL) {
          }
)

setMethod("myread",
          signature="NpdeData",
          function(object, dat = NULL) {
          }
)

setMethod("myread",
          signature=c("NpdeData","data.frame"),
          function(object, dat = NULL) {
          }
)

# Testing
theo.saemix<-read.table(file.path(datDir, "theo.saemix.tab"), header = T)
