rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)

source('../fun_coda_dic.R')
obs <- read.table("../epidural.txt",header=T)
attach(obs)

Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
is.vector(Ntol)
N0 <- n000+n001+n010+n011
N1 <- n100+n101+n110+n111
is.vector(N0)
R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
n <- length(Ntol)
pi <- pi
r <- (n100+n101+n110+n111)/Ntol
data <- list(N0=N0, N1=N1, R=R, n=n, pi=pi)

seed <- 2017
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 1
set.seed(seed)
init.seeds <- sample(1:1000000, n.chains)

### set initials ##########################
init.jags <- vector("list", n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}


################# Model 4 a) random delta_z1 delta_z2 delta_s ####################
modelString="
model{
for (i in 1:n) {

prob[i, 1] <- (pi_n[i]*(1-s1[i]) + pi_c[i]*(1-v1[i]))
prob[i, 2] <- (pi_n[i]*s1[i] + pi_c[i]*v1[i])
prob[i, 3] <- (pi_a[i]*(1-b1[i]))
prob[i, 4] <- (pi_a[i]*b1[i])
prob[i, 5] <- (pi_n[i]*(1-s1[i]))
prob[i, 6] <- (pi_n[i]*s1[i])
prob[i, 7] <- (pi_c[i]*(1-u1[i])+pi_a[i]*(1-b1[i]))
prob[i, 8] <- (pi_c[i]*u1[i]+pi_a[i]*b1[i])

R[i, 1:4] ~ dmulti(prob[i, 1:4], N0[i])
R[i, 5:8] ~ dmulti(prob[i, 5:8], N1[i])

probit(u1[i]) <- alpha_u
probit(v1[i]) <- alpha_v

z1[i] <- alpha_z1 + delta_z1[i]
delta_z1[i] ~ dnorm(0, tau_z1)
z2[i] <- alpha_z2 + delta_z2[i]
delta_z2[i] ~ dnorm(0, tau_z2)
pi_n[i] <- exp(z1[i])/(1+exp(z1[i])+exp(z2[i]))
pi_a[i] <- exp(z2[i])/(1+exp(z1[i])+exp(z2[i]))
pi_c[i] <- 1-pi_a[i]-pi_n[i]

logit(s1[i]) <- alpha_s1 + delta_s[i]
delta_s[i] ~ dnorm(0, tau_s)
logit(b1[i]) <- alpha_b1
} 

CACE <- phi(alpha_u)-phi(alpha_v)

pin <- exp(alpha_z1)/(1+exp(alpha_z1)+exp(alpha_z2))
pia <- exp(alpha_z2)/(1+exp(alpha_z1)+exp(alpha_z2))
pic <- 1-pia-pin
u1out <- phi(alpha_u)
v1out <- phi(alpha_v)
s1out <- ilogit(alpha_s1/sqrt(1+(16*sqrt(3)/(15*pi))^2*(sigma_s^2)))
b1out <- ilogit(alpha_b1)

# priors
alpha_z1 ~  dnorm(0, 0.16)
tau_z1 ~ dgamma(2, 2)
sigma_z1 <- 1/sqrt(tau_z1)
alpha_z2 ~ dnorm(0, 0.16)
tau_z2 ~ dgamma(2, 2)
sigma_z2 <- 1/sqrt(tau_z2)

alpha_s1 ~  dnorm(0, 0.25)
tau_s ~ dgamma(2, 2)
sigma_s <- 1/sqrt(tau_s)
alpha_b1 ~  dnorm(0, 0.25)

alpha_u ~  dnorm(0, 0.25)
alpha_v ~  dnorm(0, 0.25)
}"
  
writeLines(modelString, con="../models/epidural_m4a_integration.txt")


params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pi_n", "pi_a", "pi_c", "pin", "pia", "pic", 
            "alpha_z1", "alpha_z2", "sigma_z1", "sigma_z2", "sigma_s") 
jags.m4a <- jags.model(file="../models/epidural_m4a_integration.txt", data=data, inits=init.jags, 
                       n.chains=n.chains, n.adapt=n.adapt)
update(jags.m4a, n.iter=n.burnin) # burn in
samps.m4a <- coda.samples.dic(jags.m4a, variable.names=params, n.iter=n.iter, thin=n.thin)
summary(samps.m4a$samples)[[1]] #  [[1]] is mean sd,  [[2]] is quantiles  

samps.m4a$samples[1:5, ]  # there are 3 chains, each parameter is one column
samples <- rbind(samps.m4a$samples[[1]], samps.m4a$samples[[2]], samps.m4a$samples[[3]])

sampledata <- as.data.frame(samples)
sampledata[1:5,]

sampledata$mean_pia <- rowMeans(sampledata[,5:14])
sampledata$mean_pic <- rowMeans(sampledata[,15:24])
sampledata$mean_pin <- rowMeans(sampledata[,25:34])


sampledata <- sampledata[, -c(5:37)]
sampledata[1:5,]
table_m4a <- rbind(
    quantile(sampledata[, 1], c(.025, .50,  .975)), # CACE
    quantile(sampledata[, 11], c(.025, .50,  .975)), # pia
    quantile(sampledata[, 13], c(.025, .50,  .975)), # pin
    quantile(sampledata[, 12], c(.025, .50,  .975)), # pic
    quantile(sampledata[, 5], c(.025, .50,  .975)), # s1
    quantile(sampledata[, 4], c(.025, .50,  .975)), # b1
    quantile(sampledata[, 9], c(.025, .50,  .975)), # u1
    quantile(sampledata[, 10], c(.025, .50,  .975)), # v1
    quantile(sampledata[, 2], c(.025, .50,  .975)), # alpha_z1
    quantile(sampledata[, 3], c(.025, .50,  .975)), # alpha_z2
    quantile(sampledata[, 7], c(.025, .50,  .975)), # sigma_z1
    quantile(sampledata[, 8], c(.025, .50,  .975)), # sigma_z2
    quantile(sampledata[, 6], c(.025, .50,  .975))  # sigma_s
)
rownames(table_m4a) <- c("CACE", "pia", "pin", "pic", "s1", "b1", "u1", "v1", 
                         "alpha_z1", "alpha_z2", "sigma_z1", "sigma_z2", "sigma_s")

write.table(table_m4a, "table_m4a.txt", sep="\t")
