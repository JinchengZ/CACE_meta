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

# ################# Model b_1 full #########################################################
### NOTE: this model is the same as the full model in forward, so don't need to run.
data_full <- list(N0=N0, N1=N1, R=R, n=n,
                   II = structure(.Data = c(1, 0, 0, 1), .Dim = c(2, 2))  )
# params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
#             "pic", "pin", "pia", "Sigma_Z", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
# jags.mfull <- jags.model(file="m_full.txt", data=data_full, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
# update(jags.mfull, n.iter=n.burnin) # burn in
# samps.mfull <- coda.samples.dic(jags.mfull, variable.names=params, n.iter=n.iter, thin=n.thin)
# 
# write.table(summary(samps.mfull$samples)[[1]], "../outputs/out1_full.txt", sep="\t")
# write.table(summary(samps.mfull$samples)[[2]], "../outputs/out2_full.txt", sep="\t")
# 
# out3_full <- rep(NA, 3)
# out3_full[1] <- samps.mfull$dic[[1]] # mean deviance
# out3_full[2] <- samps.mfull$dic[[2]] # pD
# out3_full[3] <- samps.mfull$dic[[1]] + samps.mfull$dic[[2]]  # DIC
# write.table(out3_full, "../outputs/out3_full.txt", sep="\t")

################# Model b_2a drop delta_s #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "Sigma_Z", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_2a <- jags.model(file="../models/bm_2a.txt", data=data_full, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2a, n.iter=n.burnin) # burn in
samps.bm_2a <- coda.samples.dic(jags.bm_2a, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2a$samples)[[1]], "../outputs/out1_bm_2a.txt", sep="\t")
write.table(summary(samps.bm_2a$samples)[[2]], "../outputs/out2_bm_2a.txt", sep="\t")

out3_bm_2a <- rep(NA, 3)
out3_bm_2a[1] <- samps.bm_2a$dic[[1]] # mean deviance
out3_bm_2a[2] <- samps.bm_2a$dic[[2]] # pD
out3_bm_2a[3] <- samps.bm_2a$dic[[1]] + samps.bm_2a$dic[[2]]  # DIC
write.table(out3_bm_2a, "../outputs/out3_bm_2a.txt", sep="\t")

################# Model b_2b drop delta_b #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "Sigma_Z", "sigma_s", "sigma_u", "sigma_v") 
jags.bm_2b <- jags.model(file="../models/bm_2b.txt", data=data_full, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2b, n.iter=n.burnin) # burn in
samps.bm_2b <- coda.samples.dic(jags.bm_2b, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2b$samples)[[1]], "../outputs/out1_bm_2b.txt", sep="\t")
write.table(summary(samps.bm_2b$samples)[[2]], "../outputs/out2_bm_2b.txt", sep="\t")

out3_bm_2b <- rep(NA, 3)
out3_bm_2b[1] <- samps.bm_2b$dic[[1]] # mean deviance
out3_bm_2b[2] <- samps.bm_2b$dic[[2]] # pD
out3_bm_2b[3] <- samps.bm_2b$dic[[1]] + samps.bm_2b$dic[[2]]  # DIC
write.table(out3_bm_2b, "../outputs/out3_bm_2b.txt", sep="\t")

################# Model b_2c drop delta_u #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "Sigma_Z", "sigma_s", "sigma_b", "sigma_v") 
jags.bm_2c <- jags.model(file="../models/bm_2c.txt", data=data_full, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2c, n.iter=n.burnin) # burn in
samps.bm_2c <- coda.samples.dic(jags.bm_2c, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2c$samples)[[1]], "../outputs/out1_bm_2c.txt", sep="\t")
write.table(summary(samps.bm_2c$samples)[[2]], "../outputs/out2_bm_2c.txt", sep="\t")

out3_bm_2c <- rep(NA, 3)
out3_bm_2c[1] <- samps.bm_2c$dic[[1]] # mean deviance
out3_bm_2c[2] <- samps.bm_2c$dic[[2]] # pD
out3_bm_2c[3] <- samps.bm_2c$dic[[1]] + samps.bm_2c$dic[[2]]  # DIC
write.table(out3_bm_2c, "../outputs/out3_bm_2c.txt", sep="\t")

################# Model b_2d drop delta_v #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "Sigma_Z", "sigma_s", "sigma_b", "sigma_u") 
jags.bm_2d <- jags.model(file="../models/bm_2d.txt", data=data_full, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2d, n.iter=n.burnin) # burn in
samps.bm_2d <- coda.samples.dic(jags.bm_2d, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2d$samples)[[1]], "../outputs/out1_bm_2d.txt", sep="\t")
write.table(summary(samps.bm_2d$samples)[[2]], "../outputs/out2_bm_2d.txt", sep="\t")

out3_bm_2d <- rep(NA, 3)
out3_bm_2d[1] <- samps.bm_2d$dic[[1]] # mean deviance
out3_bm_2d[2] <- samps.bm_2d$dic[[2]] # pD
out3_bm_2d[3] <- samps.bm_2d$dic[[1]] + samps.bm_2d$dic[[2]]  # DIC
write.table(out3_bm_2d, "../outputs/out3_bm_2d.txt", sep="\t")

################# Model b_2e drop pho_z #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_2e <- jags.model(file="m7a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2e, n.iter=n.burnin) # burn in
samps.bm_2e <- coda.samples.dic(jags.bm_2e, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2e$samples)[[1]], "../outputs/out1_bm_2e.txt", sep="\t")
write.table(summary(samps.bm_2e$samples)[[2]], "../outputs/out2_bm_2e.txt", sep="\t")

out3_bm_2e <- rep(NA, 3)
out3_bm_2e[1] <- samps.bm_2e$dic[[1]] # mean deviance
out3_bm_2e[2] <- samps.bm_2e$dic[[2]] # pD
out3_bm_2e[3] <- samps.bm_2e$dic[[1]] + samps.bm_2e$dic[[2]]  # DIC
write.table(out3_bm_2e, "../outputs/out3_bm_2e.txt", sep="\t")

################# Model b_2f drop delta_z1 #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z2", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_2f <- jags.model(file="../models/bm_2f.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2f, n.iter=n.burnin) # burn in
samps.bm_2f <- coda.samples.dic(jags.bm_2f, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2f$samples)[[1]], "../outputs/out1_bm_2f.txt", sep="\t")
write.table(summary(samps.bm_2f$samples)[[2]], "../outputs/out2_bm_2f.txt", sep="\t")

out3_bm_2f <- rep(NA, 3)
out3_bm_2f[1] <- samps.bm_2f$dic[[1]] # mean deviance
out3_bm_2f[2] <- samps.bm_2f$dic[[2]] # pD
out3_bm_2f[3] <- samps.bm_2f$dic[[1]] + samps.bm_2f$dic[[2]]  # DIC
write.table(out3_bm_2f, "../outputs/out3_bm_2f.txt", sep="\t")

################# Model b_2g drop delta_z2 #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_2g <- jags.model(file="../models/bm_2g.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_2g, n.iter=n.burnin) # burn in
samps.bm_2g <- coda.samples.dic(jags.bm_2g, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_2g$samples)[[1]], "../outputs/out1_bm_2g.txt", sep="\t")
write.table(summary(samps.bm_2g$samples)[[2]], "../outputs/out2_bm_2g.txt", sep="\t")

out3_bm_2g <- rep(NA, 3)
out3_bm_2g[1] <- samps.bm_2g$dic[[1]] # mean deviance
out3_bm_2g[2] <- samps.bm_2g$dic[[2]] # pD
out3_bm_2g[3] <- samps.bm_2g$dic[[1]] + samps.bm_2g$dic[[2]]  # DIC
write.table(out3_bm_2g, "../outputs/out3_bm_2g.txt", sep="\t")


############## save out3 tables ###################
setwd('/linux/jinchenz/epidural')
out3_bm_1 <- read.table('../outputs/out3_full.txt')
out3_bm_2a <- read.table('../outputs/out3_bm_2a.txt')
out3_bm_2b <- read.table('../outputs/out3_bm_2b.txt')
out3_bm_2c <- read.table('../outputs/out3_bm_2c.txt')
out3_bm_2d <- read.table('../outputs/out3_bm_2d.txt')
out3_bm_2e <- read.table('../outputs/out3_bm_2e.txt')
out3_bm_2f <- read.table('../outputs/out3_bm_2f.txt')
out3_bm_2g <- read.table('../outputs/out3_bm_2g.txt')

out3_bm12 <- cbind(out3_bm_1, out3_bm_2a, out3_bm_2b, out3_bm_2c, 
                  out3_bm_2d, out3_bm_2e, out3_bm_2f, out3_bm_2g)

write.table(out3_bm12, "../outputs/out3_bm12.txt", sep="\t")
