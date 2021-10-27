############################### Simulation - MCMC kernels (E-step) #############################

estep<-function(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed)
	
	# Function to perform MCMC simulation
	nb.etas<-length(varList$ind.eta)
	domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
	omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
	omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
	chol.omega<-try(chol(omega.eta))
	somega<-solve(omega.eta)
	
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
	
	VK<-rep(c(1:nb.etas),2)
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	
	U.y<-compute.LLy(phiM,Uargs,Dargs,DYF,varList$pres)

	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM
	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		Uc.y<-compute.LLy(phiMc,Uargs,Dargs,DYF,varList$pres)
		deltau<-Uc.y-U.y
		ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
		etaM[ind,]<-etaMc[ind,]
		U.y[ind]<-Uc.y[ind]
	}
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	
	# Second stage
	
	if(opt$nbiter.mcmc[2]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		for (u in 1:opt$nbiter.mcmc[2]) {
			for(vk2 in 1:nb.etas) {
				etaMc<-etaM
				#				cat('vk2=',vk2,' nrs2=',nrs2,"\n")
				etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				Uc.y<-compute.LLy(phiMc,Uargs,Dargs,DYF,varList$pres)
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,]<-etaMc[ind,]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}
		}
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
	}
	
	if(opt$nbiter.mcmc[3]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-kiter%%(nb.etas-1)+2
		if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
		for (u in 1:opt$nbiter.mcmc[3]) {
			if(nrs2<nb.etas) {
				vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
				nb.iter2<-nb.etas
			} else {
				vk<-0:(nb.etas-1)
				#        if(nb.etas==1) vk<-c(0)
				nb.iter2<-1
			}
			for(k2 in 1:nb.iter2) {
				vk2<-VK[k2+vk]
				etaMc<-etaM
				etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2])
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				Uc.y<-compute.LLy(phiMc,Uargs,Dargs,DYF,varList$pres)
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				ind<-which(deltu<(-log(runif(Dargs$NM))))
				etaM[ind,]<-etaMc[ind,]
				#        if(kiter<20 | (kiter>150 & kiter<170)) {
				#        	cat("kiter=",kiter,length(ind),"  varList$ind.eta=",varList$ind.eta,"  nrs2=",nrs2,"\n")
				#        	print(head(etaMc))
				#        }
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}
		}
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
	}


	if(opt$nbiter.mcmc[4]>0 & kiter<opt$nbiter.map) {
		etaMc<-etaM
		propc <- U.eta
		prop <- U.eta
		phi.map<-mean.phi
	  	i1.omega2<-varList$ind.eta
	    iomega.phi1<-solve(omega.eta[i1.omega2,i1.omega2])

	 	# Setup for MAP calculation (MAP is identical no matter the chain)
	  	id<-Dargs$IdM[1:Dargs$nobs]
	  	xind<-Dargs$XM[1:Dargs$nobs,]
	  	yobs<-Dargs$yM[1:Dargs$nobs]
	  	id.list<-unique(id)

	  	if(Dargs$modeltype=="structural"){
	  		for(i in 1:length(id.list)) {
  	  	  		isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM[i,i1.omega2]
			    phii<-phiM[i,]
			    phi1<-phii[i1.omega2]
			    phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=Dargs$transform.par, model=Dargs$structural.model, pres=varList$pres, err=Dargs$error.model)
			    phi.map[i,i1.omega2]<-phi1.opti$par
			}
			
			# Repeat the map nchains time
			phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ]

		 	map.psi<-transphi(phi.map,Dargs$transform.par)
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)
			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map - mean.phiM

			fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
			gradf <- matrix(0L, nrow = length(fpred1), ncol = nb.etas) 

			## Compute gradient of structural model (gradf)
			for (j in 1:nb.etas) {
				psi_map2 <- psi_map
				psi_map2[,j] <- psi_map[,j]+psi_map[,j]/1000
				fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
				fpred2<-Dargs$structural.model(psi_map2, Dargs$IdM, Dargs$XM)
				for (i in 1:(Dargs$NM)){
					r = which(Dargs$IdM==i)
					gradf[r,j] <- (fpred2[r] - fpred1[r])/(psi_map[i,j]/1000)
				}
			}


			## Compute gradient of mapping (psi to phi) function (gradh)
			gradh <- list(omega.eta,omega.eta)
			for (i in 1:Dargs$NM){
				gradh[[i]] <- gradh[[1]]
			}
			for (j in 1:nb.etas) {
				phi_map2 <- phi_map
				phi_map2[,j] <- phi_map[,j]+phi_map[,j]/1000
				psi_map2 <- transphi(phi_map2,Dargs$transform.par) 
				for (i in 1:(Dargs$NM)){
					gradh[[i]][,j] <- (psi_map2[i,] - psi_map[i,])/(phi_map[i,]/1000)
				}
			}
			
			## Calculation of the covariance matrix of the proposal
			Gamma <- chol.Gamma <- inv.chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
			for (i in 1:(Dargs$NM)){
				r = which(Dargs$IdM==i)
		        temp <- gradf[r,]%*%gradh[[i]]
				Gamma[[i]] <- solve(t(temp)%*%temp/(varList$pres[1])^2+solve(omega.eta))
				chol.Gamma[[i]] <- chol(Gamma[[i]])
				inv.chol.Gamma[[i]] <- solve(chol.Gamma[[i]])
				inv.Gamma[[i]] <- solve(Gamma[[i]])
			}

		} else {
		  for(i in 1:length(id.list)) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM[i,i1.omega2]
			    phii<-phiM[i,]
			    phi1<-phii[i1.omega2]
			    phi1.opti<-optim(par=phi1, fn=conditional.distribution_d, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=Dargs$transform.par, model=Dargs$structural.model)
			    phi.map[i,i1.omega2]<-phi1.opti$par
			}
			#rep the map nchains time
			phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ] 
		  	map.psi<-transphi(phi.map,Dargs$transform.par)
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)

			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map[,varList$ind.eta] - mean.phiM[,varList$ind.eta]
			
			#gradient at the map estimation
			gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas) 

			for (j in 1:nb.etas) {
				phi_map2 <- phi_map
				phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
				psi_map2 <- transphi(phi_map2,Dargs$transform.par) 
				fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM]<- fpred1
				l1<-colSums(DYF)
				fpred2<-Dargs$structural.model(psi_map2, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM]<- fpred2
				l2<-colSums(DYF)

				for (i in 1:(Dargs$NM)){
					gradp[i,j] <- (l2[i] - l1[i])/(phi_map[i,j]/100)
				}
			}

			#calculation of the covariance matrix of the proposal
			fpred<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- fpred
			denom <- colSums(DYF)
			
			Gamma <- chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
			z <- matrix(0L, nrow = length(fpred), ncol = 1) 
			for (i in 1:(Dargs$NM)){
				Gamma[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
				chol.Gamma[[i]] <- chol(Gamma[[i]])
				inv.Gamma[[i]] <- solve(Gamma[[i]])
			}

		}

		etaM <- eta_map
	  	phiM<-etaM+mean.phiM
	  	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	  	U.y<-compute.LLy(phiM,Uargs,Dargs,DYF,varList$pres)

	  	for (u in 1:opt$nbiter.mcmc[4]) {
			
			#generate candidate eta with new proposal
			for (i in 1:(Dargs$NM)){
				Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
				etaMc[i,varList$ind.eta]<- eta_map[i,varList$ind.eta] + Mi
			}

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc[,varList$ind.eta]
			Uc.y<-compute.LLy(phiM,Uargs,Dargs,DYF,varList$pres)
			Uc.eta<-0.5*rowSums(etaMc[,varList$ind.eta]*(etaMc[,varList$ind.eta]%*%somega))

			for (i in 1:(Dargs$NM)){
				propc[i] <- 0.5*rowSums((etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
				prop[i] <- 0.5*rowSums((etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
			}
			
			deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			etaM[ind,varList$ind.eta]<-etaMc[ind,varList$ind.eta]
			U.y[ind]<-Uc.y[ind]
			U.eta[ind]<-Uc.eta[ind]

  		}
	}

	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM))
}
