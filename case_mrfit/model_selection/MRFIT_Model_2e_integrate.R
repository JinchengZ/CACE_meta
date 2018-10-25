source("../start.R")

seed <- 2017
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 40000
n.thin <- 1
set.seed(seed)
init.seeds <- sample(1:100000, n.chains)

### set initials ##########################
init.jags <- vector("list", n.chains)
for(i in 1:n.chains){
  init.jags[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
}



################# Model 2 e) random delta_z1 ###########################

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

z1[i] <- alpha_z1 + delta_z1[i]
delta_z1[i] ~ dnorm(0, tau_z1)
z2[i] <- alpha_z2 
pi_n[i] <- exp(z1[i])/(1+exp(z1[i])+exp(z2[i]))
pi_a[i] <- exp(z2[i])/(1+exp(z1[i])+exp(z2[i]))
pi_c[i] <- 1-pi_a[i]-pi_n[i]

logit(s1[i]) <- alpha_s1 
logit(b1[i]) <- alpha_b1
} 

CACE <- phi(alpha_u)-phi(alpha_v)
CACEbystudy <- u1-v1

pin <- exp(alpha_z1)/(1+exp(alpha_z1)+exp(alpha_z2))
pia <- exp(alpha_z2)/(1+exp(alpha_z1)+exp(alpha_z2))
pic <- 1-pia-pin
u1out <- phi(alpha_u)
v1out <- phi(alpha_v)
s1out <- ilogit(alpha_s1)
b1out <- ilogit(alpha_b1)

# priors
alpha_z1~  dnorm(0, 0.16)
tau_z1 ~ dgamma(2, 2)
sigma_z1 <- 1/sqrt(tau_z1)
alpha_z2 ~ dnorm(0, 0.16)

alpha_s1 ~  dnorm(0, 0.25)
alpha_b1 ~  dnorm(0, 0.25)

alpha_u ~  dnorm(0, 0.25)
alpha_v ~  dnorm(0, 0.25)
}"
  
writeLines(modelString, con="../models/MRFIT_m2e_integration.txt")

params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pi_n", "pi_a", "pi_c", "pin", "pia", "pic", "alpha_z1", "alpha_z2", "sigma_z1") 
jags.m2e <- jags.model(file="../models/MRFIT_m2e_integration.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2e, n.iter=n.burnin) # burn in
samps.m2e <- coda.samples.dic(jags.m2e, variable.names=params, n.iter=n.iter, thin=n.thin)
summary(samps.m2e$samples)[[1]] #  [[1]] is mean sd,  [[2]] is quantiles  

samps.m2e$samples[1:5, ]  # there are 3 chains, each parameter is one column
samples <- rbind(samps.m2e$samples[[1]], samps.m2e$samples[[2]], samps.m2e$samples[[3]])

sampledata <- as.data.frame(samples)

sampledata$mean_pia <- rowMeans(sampledata[,5:26])
sampledata$mean_pic <- rowMeans(sampledata[,27:48])
sampledata$mean_pin <- rowMeans(sampledata[,49:70])

sampledata <- sampledata[, -c(5:73)]
table_m2e <- rbind(
    quantile(sampledata[, 1], c(.025, .50,  .975)), # CACE
    quantile(sampledata[, 9], c(.025, .50,  .975)), # pia
    quantile(sampledata[, 11], c(.025, .50,  .975)), # pin
    quantile(sampledata[, 10], c(.025, .50,  .975)), # pic
    quantile(sampledata[, 5], c(.025, .50,  .975)), # s1
    quantile(sampledata[, 4], c(.025, .50,  .975)), # b1
    quantile(sampledata[, 7], c(.025, .50,  .975)), # u1
    quantile(sampledata[, 8], c(.025, .50,  .975)), # v1
    quantile(sampledata[, 6], c(.025, .50,  .975))  # sigma_z1
)
rownames(table_m2e) <- c("CACE", "pia", "pin", "pic", "s1", "b1", "u1", "v1", "sigma_z1")

write.table(table_m2e, "../outputs/table_m2e.txt", sep="\t")
