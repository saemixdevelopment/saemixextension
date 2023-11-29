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

test_that("Simulate from a model fit for the dataset in the object, without residual variability", {
  simdat<-simulate(theo.fit1, nsim=2, res.var=FALSE)
  expect_equal(dim(theo.fit1@sim.data@data)[1],0)
  expect_equal(dim(simdat@sim.data@datasim)[1],theo.fit1@data@ntot.obs*2)
  par(mfrow=c(1,3))
  plot(simdat@data, new=FALSE)
  plot(simdat@sim.data, irep=1, prediction=TRUE, new=FALSE) # plot predictions instead of simulations
  plot(simdat@sim.data, irep=2, prediction=TRUE, new=FALSE)
})


test_that("Simulate parameters from a model fit for the dataset in the object", {
  simdat<-simulate(theo.fit1, nsim=50, predictions=FALSE)
  expect_equal(dim(theo.fit1@sim.data@data)[1],0)
  expect_equal(dim(simdat@sim.data@datasim)[1],theo.fit1@data@ntot.obs*simdat@sim.data@nsim)
  expect_equal(dim(simdat@sim.data@datasim)[2],2) # no predictions or simulations, just idsim and irep
  # Graph overlaying the histogram of individual parameters and the density from 20 simulations of the same parameters
  par(mfrow=c(1,3))
  for(icol in 1:3) {
    hist(theo.fit1@results@map.psi[,icol],freq=F)
    lines(density(simdat@sim.data@sim.psi[,(icol+1)]), lwd = 2, col = 'red')
  }
})

# TODO
context("Testing simulate for a binary model \n")


if(FALSE) {
  par(mfrow=c(1,3))
  plot(simdat@data)
  plot(simdat@sim.data, irep=1)
  plot(simdat@sim.data, irep=2)
  
  plot.saemix.mirrorplot(theo.fit1)
}
