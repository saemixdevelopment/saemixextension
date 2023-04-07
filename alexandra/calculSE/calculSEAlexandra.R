library("survival")

# Essai 1 : loi exponentielle

N = 100
h0 = 0.05

T = rexp(N,h0)
delta = T<30
T = pmin(T,30)

mod = survreg(Surv(T,delta)~1,dist="exponential")

h0 	# Valeur exacte
sum(delta)/sum(T) # Formule explicite (maximum de la log-vraisemblance)
exp(-coefficients(mod)) # valeur estimée par survreg

sum(delta)/sum(T)^2 # Variance asymptotique (formule explicite)
vcov(mod)[1,1]*exp(-2*coefficients(mod)) # variance estimée par survreg (méthode delta)

# Essai 2 : loi de Gompertz (deux paramètres, forme proche de ce que tu utilises)

library("flexsurv")

N = 100000
h0 = 0.02
alpha = 0.05

T = rgompertz(N,alpha,h0)
delta = T<30
T = pmin(T,30)

mod = flexsurvreg(Surv(T,delta)~1,dist="gompertz")

c(h0,alpha)	# coefficients exacts
c(exp(coef(mod)[2]),coef(mod)[1]) 	# valeur estimé

x = exp(alpha*T)
M11 = h0/alpha^3 * sum(-2+2*x+alpha^2*T^2*x-2*alpha*T*x)
M12 = M21 = 1/alpha^2 * sum(1-x+alpha*T*x)
M22 = 1/h0^2 *sum(delta)
V = solve(matrix(c(M11,M12,M21,M22),nrow=2))	# Matrice de variance asymptotique (formule)

S = matrix(c(1,0,0,h0),nrow=2)
V2 = S%*%vcov(mod)%*%t(S)	# Matrice de variance estimée (méthode delta car la paramétrisation est différente)




