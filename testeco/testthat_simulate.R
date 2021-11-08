context("Testing simulate for a structural model \n")

test_that("Simulate from a model fit for the dataset in the object", {
  simdat<-simulate(theo.fit1, nsim=2)
  expect_equal(dim(theo.fit1@sim.data@data)[1],0)
  expect_equal(dim(simdat@sim.data@datasim)[1],theo.fit1@data@ntot.obs*2)
  par(mfrow=c(1,3))
  plot(simdat@data, new=FALSE)
  plot(simdat@sim.data, irep=1, new=FALSE)
  plot(simdat@sim.data, irep=2, new=FALSE)
})

if(FALSE) {
  par(mfrow=c(1,3))
  plot(simdat@data)
  plot(simdat@sim.data, irep=1)
  plot(simdat@sim.data, irep=2)
  
  plot.saemix.mirrorplot(theo.fit1)
}
