# Mails Mélanie 07/05
cov.eta <- fit1@results@omega
se.cov.eta <- fit1@results@se.cov

omega<-sqrt(diag(cov.eta))
se.omega<- diag(se.cov.eta)/(2*omega)

corr.eta <- t(t(cov.eta/omega)/omega)
se.corr.eta <- se.cov.eta/omega%*%t(omega) - t(t(cov.eta)*se.omega)/omega^2%*%t(omega) - cov.eta*se.omega/omega%*%t(omega)^2

Soit, sans création de variable intermédiaire :

corr.eta <- t(t( fit1@results@omega / sqrt(diag(fit1@results@omega)) )/ sqrt(diag(fit1@results@omega)) )


se.corr.eta <- sqrt( (se.cov.eta/omega%*%t(omega))^2 + (t(t(cov.eta)*se.omega)/omega^2%*%t(omega))^2 + (cov.eta*se.omega/omega%*%t(omega)^2)^2)

se.corr.eta <- sqrt( (fit1@results@se.cov / sqrt(diag(fit1@results@omega)) %*%t( sqrt(diag(fit1@results@omega))))^2 + (t(t( fit1@results@omega )* diag(fit1@results@se.cov )/(2*sqrt(diag(fit1@results@omega))) )/ sqrt(diag(fit1@results@omega))^2%*%t( sqrt(diag(fit1@results@omega))))^2 + (fit1@results@omega * diag(fit1@results@se.cov )/(2*sqrt(diag(fit1@results@omega))))/ sqrt(diag(fit1@results@omega)) %*%t( sqrt(diag(fit1@results@omega)) )^2)^2)

