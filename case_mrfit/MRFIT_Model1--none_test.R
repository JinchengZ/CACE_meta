setwd('/linux/jinchenz/research')
source("start.R")
#########################################
####  1. no random effect
modelString="
model{
for (i in 1:n) {

R[i, ] ~ dmulti(prob[i, 1:8], Ntol[i])

prob[i, 1] <- (1-r[i])*(pi_n[i]*(1-s1[i]) + pi_c[i]*(1-v1[i]))
prob[i, 2] <- (1-r[i])*(pi_n[i]*s1[i] + pi_c[i]*v1[i])
prob[i, 3] <- (1-r[i])*(pi_a[i]*(1-b1[i]))
prob[i, 4] <- (1-r[i])*(pi_a[i]*b1[i])
prob[i, 5] <- r[i]*(pi_n[i]*(1-s1[i]))
prob[i, 6] <- r[i]*(pi_n[i]*s1[i])
prob[i, 7] <- r[i]*(pi_c[i]*(1-u1[i])+pi_a[i]*(1-b1[i]))
prob[i, 8] <- r[i]*(pi_c[i]*u1[i]+pi_a[i]*b1[i])

probit(u1[i]) <- alpha_u
probit(v1[i]) <- alpha_v
# u1[i] <- phi(alpha_u)
# v1[i] <- phi(alpha_v)

z1[i] <- alpha_z1
z2[i] <- alpha_z2 
pi_n[i] <- exp(z1[i])/(1+exp(z1[i])+exp(z2[i]))
pi_a[i] <- exp(z2[i])/(1+exp(z1[i])+exp(z2[i]))
pi_c[i] <- 1-pi_a[i]-pi_n[i]

logit(s1[i]) <- alpha_s1 
logit(b1[i]) <- alpha_b1
} 

CACE <- phi(alpha_u)-phi(alpha_v)

pin <- exp(alpha_z1)/(1+exp(alpha_z1)+exp(alpha_z2))
pia <- exp(alpha_z2)/(1+exp(alpha_z1)+exp(alpha_z2))
pic <- 1-pia-pin
u1out <- phi(alpha_u)
v1out <- phi(alpha_v)
s1out <- ilogit(alpha_s1)
b1out <- ilogit(alpha_b1)

# priors
alpha_z1 ~  dnorm(0, 0.16)
alpha_z2 ~ dnorm(0, 0.16)

alpha_s1 ~  dnorm(0, 0.25)
alpha_b1 ~  dnorm(0, 0.25)

alpha_u ~  dnorm(0, 0.25)
alpha_v ~  dnorm(0, 0.25)
}"
  
writeLines(modelString, con="MRFIT_m1.txt")

params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia") 


seed <- 123
n.chains <- 3
n.adapt <- 1000
n.burnin <- 1000
n.iter <- 5000
n.thin <- 1

init.jags <- vector("list", n.chains)
set.seed(seed)
init.seeds <- sample(1:100000, n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}
################# rjags code #########################################################
jags.m1 <- jags.model(file="MRFIT_m1.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m1, n.iter=n.burnin) # burn in
samps.m1 <- coda.samples.dic(jags.m1, variable.names=params, n.iter=n.iter, thin=n.thin)

setwd('/linux/jinchenz/research')
write.table(summary(samps.m1$samples)[[1]], "out1_1.txt", sep="\t")
write.table(summary(samps.m1$samples)[[2]], "out2_1.txt", sep="\t")

out3_1 <- rep(NA, 3)
out3_1[1] <- samps.m1$dic[[1]] # mean deviance
out3_1[2] <- samps.m1$dic[[2]] # pD
out3_1[3] <- samps.m1$dic[[1]] + samps.m1$dic[[2]]  # DIC
write.table(out3_1, "out3_1.txt", sep="\t")
