testExist <- function(x) {
  if(exists("x")) cat("Test passed\n")
  if(!exists("x")) cat("x does not exist\n") else  print(x)
}

x1<-2
testExist(x1)

testExist(x)

testExist2 <- function(x) {
  if(is(x,"character")) {
    if(!exists("x")) cat("x does not exist\n") else return(get(x))
  } else cat("Pass a character\n")
}

x1<-2
testExist2(x1) # should ask to Pass a character
testExist2("x1") # should return the value of x1 (2)
try(testExist2(thisdoesnotexist)) # should give an error message

################# Testing whether an object actually exists when it's passed to a function - tricky (pby side-effects when nested within loops)

testExist3 <- function(x) {
  if(missing(x)) return("missing") else print(x)
}

x1<-2
testExist3(x1) # should print 2
testExist3() # should print missing
try(testExist3(thisdoesnotexist)) # should give an error message

# Suggested solution - looks for variable in the global environment which may fail when encapsulated within eg loops or another function
## also can't pass directly eg a number
library(rlang)

f1 <- function(x){ 
  arg <- quo_name(enquo(x))
  if(exists(arg, where = .GlobalEnv)){
    return(2*x)
  } else {
    cat('variable ', '`', arg, '`', ' does not exist', sep = "")
  }
}
x1<-2
f1(thisdoesnotexist) # diagnoses non-existence correctly
f1(3) # lol, as a side-effect, this doesn't work either
f1(x1)
