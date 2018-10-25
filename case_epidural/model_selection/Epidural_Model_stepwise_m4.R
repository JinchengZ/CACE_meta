source("../start.R")

seed <- 1234
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

################# Model 4 a) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2", "sigma_s")  
jags.m4a <- jags.model(file="../models/m4a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m4a, n.iter=n.burnin) # burn in
samps.m4a <- coda.samples.dic(jags.m4a, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m4a$samples)[[1]], "../outputs/out1_4a.txt", sep="\t")
write.table(summary(samps.m4a$samples)[[2]], "../outputs/out2_4a.txt", sep="\t")

out3_4a <- rep(NA, 3)
out3_4a[1] <- samps.m4a$dic[[1]] # mean deviance
out3_4a[2] <- samps.m4a$dic[[2]] # pD
out3_4a[3] <- samps.m4a$dic[[1]] + samps.m4a$dic[[2]]  # DIC
write.table(out3_4a, "../outputs/out3_4a.txt", sep="\t")

################# Model 4 a_a) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_s") 
jags.m4a_a <- jags.model(file="m4a_a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m4a_a, n.iter=n.burnin) # burn in
samps.m4a_a <- coda.samples.dic(jags.m4a_a, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m4a_a$samples)[[1]], "../outputs/out1_4a_a.txt", sep="\t")
write.table(summary(samps.m4a_a$samples)[[2]], "../outputs/out2_4a_a.txt", sep="\t")

out3_4a_a <- rep(NA, 3)
out3_4a_a[1] <- samps.m4a_a$dic[[1]] # mean deviance
out3_4a_a[2] <- samps.m4a_a$dic[[2]] # pD
out3_4a_a[3] <- samps.m4a_a$dic[[1]] + samps.m4a_a$dic[[2]]  # DIC
write.table(out3_4a_a, "../outputs/out3_4a_a.txt", sep="\t")

################# Model 4 a_b) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z2", "sigma_s") 
jags.m4a_b <- jags.model(file="../models/m4a_b.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m4a_b, n.iter=n.burnin) # burn in
samps.m4a_b <- coda.samples.dic(jags.m4a_b, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m4a_b$samples)[[1]], "../outputs/out1_4a_b.txt", sep="\t")
write.table(summary(samps.m4a_b$samples)[[2]], "../outputs/out2_4a_b.txt", sep="\t")

out3_4a_b <- rep(NA, 3)
out3_4a_b[1] <- samps.m4a_b$dic[[1]] # mean deviance
out3_4a_b[2] <- samps.m4a_b$dic[[2]] # pD
out3_4a_b[3] <- samps.m4a_b$dic[[1]] + samps.m4a_b$dic[[2]]  # DIC
write.table(out3_4a_b, "../outputs/out3_4a_b.txt", sep="\t")

################# Model full #########################################################
data_full <- list(N0=N0, N1=N1, R=R, n=n,
               II = structure(.Data = c(1, 0, 0, 1), .Dim = c(2, 2))  )
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "Sigma_Z", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
jags.mfull <- jags.model(file="../models/m_full.txt", data=data_full, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.mfull, n.iter=n.burnin) # burn in
samps.mfull <- coda.samples.dic(jags.mfull, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.mfull$samples)[[1]], "../outputs/out1_full.txt", sep="\t")
write.table(summary(samps.mfull$samples)[[2]], "../outputs/out2_full.txt", sep="\t")

out3_full <- rep(NA, 3)
out3_full[1] <- samps.mfull$dic[[1]] # mean deviance
out3_full[2] <- samps.mfull$dic[[2]] # pD
out3_full[3] <- samps.mfull$dic[[1]] + samps.mfull$dic[[2]]  # DIC
write.table(out3_full, "../outputs/out3_full.txt", sep="\t")



############## save out3 tables ###################
out3_4a <- read.table('../outputs/out3_4a.txt')
out3_4a_a <- read.table('../outputs/out3_4a_a.txt')
out3_4a_b <- read.table('../outputs/out3_4a_b.txt')
out3_full <- read.table('../outputs/out3_full.txt')

out4_stepwise <- cbind(out3_4a, out3_4a_a, out3_4a_b, out3_full)

write.table(out4_stepwise, "../outputs/out4_stepwise.txt", sep="\t")
