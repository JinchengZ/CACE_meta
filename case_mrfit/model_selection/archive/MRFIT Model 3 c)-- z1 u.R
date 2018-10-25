setwd('/linux/jinchenz/research')
source("start.R")
#########################################
####  3. c) random delta_z1, delta_u
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

probit(u1[i]) <- alpha_u + delta_u[i]
delta_u[i] ~ dnorm(0, tau_u)
probit(v1[i]) <- alpha_v

z1[i] <- alpha_z1 + delta_z1[i]
delta_z1[i] ~ dnorm(0, tau_z1)
z2[i] <- alpha_z2 
pi_n[i] <- exp(z1[i])/(1+exp(z1[i])+exp(z2[i]))
pi_a[i] <- exp(z2[i])/(1+exp(z1[i])+exp(z2[i]))
pi_c[i] <- 1-pi_a[i]-pi_n[i]

logit(s1[i]) <- alpha_s1
logit(b1[i]) <- alpha_b1
} 

CACE <- phi(alpha_u/sqrt(1+sigma_u^2))-phi(alpha_v)

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
tau_z1 ~ dgamma(2, 2)
sigma_z1 <- 1/sqrt(tau_z1)

alpha_s1 ~  dnorm(0, 0.25)
alpha_b1 ~  dnorm(0, 0.25)

alpha_u ~  dnorm(0, 0.25)
tau_u ~ dgamma(2, 2)
sigma_u <- 1/sqrt(tau_u)
alpha_v ~  dnorm(0, 0.25)
}"
  
writeLines(modelString, con="MRFIT_m3c.txt")

#params <- c("CACE", "alpha_z1", "alpha_z2", "alpha_s1", "alpha_b1", "alpha_u", "alpha_v") 
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_u") 

Inits1 <- list("alpha_z1"= 1, "alpha_z2"= 1, "alpha_s1"= 1, "alpha_b1"= 1, "alpha_u"=1, "alpha_v"=2, 
               "tau_z1"=0.5, "tau_u"=0.5, .RNG.name="base::Wichmann-Hill", .RNG.seed=1) 
Inits2 <- list("alpha_z1"= 0, "alpha_z2"= 0, "alpha_s1"= 0, "alpha_b1"= 0, "alpha_u"=0, "alpha_v"=1, 
               "tau_z1"=2, "tau_u"=2, .RNG.name="base::Wichmann-Hill", .RNG.seed=1)
Inits=list(Inits1,Inits2)

################# rjags code #########################################################
jags.m3c <- jags.model(file="MRFIT_m3c.txt", data=data, inits=Inits, n.chains=2, n.adapt=2000)
update(jags.m3c, n.iter=5000) # burn in
samps.m3c <- coda.samples.dic(jags.m3c, variable.names=params, n.iter=20000, thin=5)

setwd('/linux/jinchenz/research')
write.table(summary(samps.m3c$samples)[[1]], "out1_3c.txt", sep="\t")
write.table(summary(samps.m3c$samples)[[2]], "out2_3c.txt", sep="\t")

out3_3c <- rep(NA, 3)
out3_3c[1] <- samps.m3c$dic[[1]] # mean deviance
out3_3c[2] <- samps.m3c$dic[[2]] # pD
out3_3c[3] <- samps.m3c$dic[[1]] + samps.m3c$dic[[2]]  # DIC
write.table(out3_3c, "out3_3c.txt", sep="\t")

