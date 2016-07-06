N <- 10000
err <- 1

parms <- c(c1x = 0.1, c2x = 0.1, c3m = 0.1, c1m = 0.1,
           c1y = 0.2, c2y = 0.2, c3y = 0.3,
           xy = 0.5, xm = 0.2, xl = 0.6,
           my = 0.5, ly = 0.2,
           lm = 0.2)

c_probs <- c(c1 = 0.5, c2 = 0.5, c3 = 0.5)
for(i in 1:3) assign(paste0("C", i), rbinom(N, 1, c_probs[i]))

X <- rbinom(N, 1,cbind(C1, C2)%*% parms[c("c1x","c2x")])
L <- rbinom(N, 1,0.2 + cbind(X)%*% parms[c("xl")])
M <- rbinom(N, 1,cbind(C1, C3, X, L)%*% parms[c("c1m","c3m","xm","lm")])
Y <- rnorm(N, cbind(C1, C2, C3, X, L, M)%*% parms[c("c1y","c2y","c3y","xy","ly","my")], sd = err)

lm(Y ~ X )
lm(Y ~ X + C1 + C2) # total effect
lm(Y ~ X + C1 + C2 + C3 + M + L) #direct effect


