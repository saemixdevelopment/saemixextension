context("Testing auxiliary functions\n")

test_that("Error model - one response type", {
  f1<-rep(10,20)
  yt1<-rep(1,20)
  ab1<-c(1,0)
  ab2<-c(0,1)
  ab3<-c(1,1)
  expect_equal(sum(error(f1,ab1,yt1)),20)
  expect_equal(sum(error(f1,ab2,yt1)),200)
  expect_equal(sum(error(f1,ab3,yt1)),220)
})

test_that("Error model - two response types", {
  f1<-c(rep(c(5,10),10))
  yt1<-rep(c(1,2),10)
  ab1<-c(1,0,0,1)
  ab2<-c(1,0,1,1)
  ab3<-c(1,0,2,0)
  ab4<-c(0,1,1,1)
  expect_equal(sum(error(f1,ab1,yt1)),110)
  expect_equal(sum(error(f1,ab2,yt1)),120)
  expect_equal(sum(error(f1,ab3 ,yt1)),30)
  expect_equal(sum(error(f1,ab4,yt1)),160)
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
  expect_equal(ssq(ab3,y1,f1,yt1),(20*(1/11)**2+40*log(11)))
})

