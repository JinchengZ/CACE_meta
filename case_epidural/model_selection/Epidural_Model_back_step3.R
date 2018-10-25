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

################# Model b_3a drop delta_s #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_3a <- jags.model(file="bm_3a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_3a, n.iter=n.burnin) # burn in
samps.bm_3a <- coda.samples.dic(jags.bm_3a, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_3a$samples)[[1]], "../outputs/out1_bm_3a.txt", sep="\t")
write.table(summary(samps.bm_3a$samples)[[2]], "../outputs/out2_bm_3a.txt", sep="\t")

out3_bm_3a <- rep(NA, 3)
out3_bm_3a[1] <- samps.bm_3a$dic[[1]] # mean deviance
out3_bm_3a[2] <- samps.bm_3a$dic[[2]] # pD
out3_bm_3a[3] <- samps.bm_3a$dic[[1]] + samps.bm_3a$dic[[2]]  # DIC
write.table(out3_bm_3a, "../outputs/out3_bm_3a.txt", sep="\t")

################# Model b_3b drop delta_b #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2", "sigma_s", "sigma_u", "sigma_v") 
jags.bm_3b <- jags.model(file="../models/m6b.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_3b, n.iter=n.burnin) # burn in
samps.bm_3b <- coda.samples.dic(jags.bm_3b, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_3b$samples)[[1]], "../outputs/out1_bm_3b.txt", sep="\t")
write.table(summary(samps.bm_3b$samples)[[2]], "../outputs/out2_bm_3b.txt", sep="\t")

out3_bm_3b <- rep(NA, 3)
out3_bm_3b[1] <- samps.bm_3b$dic[[1]] # mean deviance
out3_bm_3b[2] <- samps.bm_3b$dic[[2]] # pD
out3_bm_3b[3] <- samps.bm_3b$dic[[1]] + samps.bm_3b$dic[[2]]  # DIC
write.table(out3_bm_3b, "../outputs/out3_bm_3b.txt", sep="\t")


################# Model b_3c drop delta_u #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2", "sigma_s", "sigma_b", "sigma_v") 
jags.bm_3c <- jags.model(file="bm_3c.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_3c, n.iter=n.burnin) # burn in
samps.bm_3c <- coda.samples.dic(jags.bm_3c, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_3c$samples)[[1]], "../outputs/out1_bm_3c.txt", sep="\t")
write.table(summary(samps.bm_3c$samples)[[2]], "../outputs/out2_bm_3c.txt", sep="\t")

out3_bm_3c <- rep(NA, 3)
out3_bm_3c[1] <- samps.bm_3c$dic[[1]] # mean deviance
out3_bm_3c[2] <- samps.bm_3c$dic[[2]] # pD
out3_bm_3c[3] <- samps.bm_3c$dic[[1]] + samps.bm_3c$dic[[2]]  # DIC
write.table(out3_bm_3c, "../outputs/out3_bm_3c.txt", sep="\t")

################# Model b_3d drop delta_v #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2", "sigma_s", "sigma_b", "sigma_u") 
jags.bm_3d <- jags.model(file="../models/m6a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_3d, n.iter=n.burnin) # burn in
samps.bm_3d <- coda.samples.dic(jags.bm_3d, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_3d$samples)[[1]], "../outputs/out1_bm_3d.txt", sep="\t")
write.table(summary(samps.bm_3d$samples)[[2]], "../outputs/out2_bm_3d.txt", sep="\t")

out3_bm_3d <- rep(NA, 3)
out3_bm_3d[1] <- samps.bm_3d$dic[[1]] # mean deviance
out3_bm_3d[2] <- samps.bm_3d$dic[[2]] # pD
out3_bm_3d[3] <- samps.bm_3d$dic[[1]] + samps.bm_3d$dic[[2]]  # DIC
write.table(out3_bm_3d, "../outputs/out3_bm_3d.txt", sep="\t")

################# Model b_3e drop delta_z1 #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z2", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_3e <- jags.model(file="bm_3e.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_3e, n.iter=n.burnin) # burn in
samps.bm_3e <- coda.samples.dic(jags.bm_3e, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_3e$samples)[[1]], "../outputs/out1_bm_3e.txt", sep="\t")
write.table(summary(samps.bm_3e$samples)[[2]], "../outputs/out2_bm_3e.txt", sep="\t")

out3_bm_3e <- rep(NA, 3)
out3_bm_3e[1] <- samps.bm_3e$dic[[1]] # mean deviance
out3_bm_3e[2] <- samps.bm_3e$dic[[2]] # pD
out3_bm_3e[3] <- samps.bm_3e$dic[[1]] + samps.bm_3e$dic[[2]]  # DIC
write.table(out3_bm_3e, "../outputs/out3_bm_3e.txt", sep="\t")

################# Model b_3f drop delta_z2 #######################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_s", "sigma_b", "sigma_u", "sigma_v") 
jags.bm_3f <- jags.model(file="bm_3f.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.bm_3f, n.iter=n.burnin) # burn in
samps.bm_3f <- coda.samples.dic(jags.bm_3f, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.bm_3f$samples)[[1]], "../outputs/out1_bm_3f.txt", sep="\t")
write.table(summary(samps.bm_3f$samples)[[2]], "../outputs/out2_bm_3f.txt", sep="\t")

out3_bm_3f <- rep(NA, 3)
out3_bm_3f[1] <- samps.bm_3f$dic[[1]] # mean deviance
out3_bm_3f[2] <- samps.bm_3f$dic[[2]] # pD
out3_bm_3f[3] <- samps.bm_3f$dic[[1]] + samps.bm_3f$dic[[2]]  # DIC
write.table(out3_bm_3f, "../outputs/out3_bm_3f.txt", sep="\t")


############## save out3 tables ###################
out3_bm_3a <- read.table('../outputs/out3_bm_3a.txt')
out3_bm_3b <- read.table('../outputs/out3_bm_3b.txt')
out3_bm_3c <- read.table('../outputs/out3_bm_3c.txt')
out3_bm_3d <- read.table('../outputs/out3_bm_3d.txt')
out3_bm_3e <- read.table('../outputs/out3_bm_3e.txt')
out3_bm_3f <- read.table('../outputs/out3_bm_3f.txt')

out3_bm3 <- cbind(out3_bm_3a, out3_bm_3b, out3_bm_3c, 
                   out3_bm_3d, out3_bm_3e, out3_bm_3f)

write.table(out3_bm3, "../outputs/out3_bm3.txt", sep="\t")
