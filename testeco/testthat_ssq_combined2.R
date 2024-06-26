context("Testing auxiliary functions - combined 2 error model (g=a^2 + b^2 f^2) \n")

test_that("Error model - one response type", {
  f1<-rep(10,20)
  yt1<-rep(1,20)
  ab1<-c(1,0)
  ab2<-c(0,1)
  ab3<-c(1,1)
  expect_equal(sum(error(f1,ab1,yt1)),20)
  expect_equal(sum(error(f1,ab2,yt1)),200)
#  expect_equal(sum(error(f1,ab3,yt1)),sum((1+f1)))
  expect_equal(sum(error(f1,ab3,yt1)),sum(sqrt(1+f1**2)))
})

test_that("Error model - two response types", {
  f1<-c(rep(c(5,10),10))
  yt1<-rep(c(1,2),10)
  ab1<-c(1,0,0,1)
  ab2<-c(1,0,1,1)
  ab3<-c(1,0,2,0)
  ab4<-c(0,1,1,1)
  expect_equal(sum(error(f1,ab1,yt1)),110)
  xt1<-sum(yt1==1)+sum(sqrt(1+f1[yt1==2]**2))
  expect_equal(sum(error(f1,ab2,yt1)),xt1)
  expect_equal(sum(error(f1,ab3,yt1)),30)
  xt2<-sum(f1[yt1==1])+sum(sqrt(1+f1[yt1==2]**2))
  expect_equal(sum(error(f1,ab4,yt1)),xt2)
})

test_that("Sum of squares", {
  f1<-rep(10,20)
  y1<-rep(11,20)
  yt1<-rep(1,20)
  ab1<-c(1,0)
  ab2<-c(0,1)
  ab3<-c(1,1)
  expect_equal(ssq(ab1,y1,f1,yt1),20)
  expect_equal(ssq(ab2,y1,f1,yt1),(0.2+40*log(10)))
  xt1<-((y1-f1)**2)/(ab3[1]**2+(ab3[2]*f1)**2) + 2*log(sqrt((ab3[1]**2+(ab3[2]*f1)**2)))
  expect_equal(ssq(ab3,y1,f1,yt1),sum(xt1))
})
