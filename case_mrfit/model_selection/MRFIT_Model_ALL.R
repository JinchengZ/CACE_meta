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

################# model 1 #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "alpha_z1", "alpha_z2") 
jags.m1 <- jags.model(file="../models/MRFIT_m1.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m1, n.iter=n.burnin) # burn in
samps.m1 <- coda.samples.dic(jags.m1, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m1$samples)[[1]], "out1_1.txt", sep="\t")
write.table(summary(samps.m1$samples)[[2]], "out2_1.txt", sep="\t")
#gelman.diag(samps.m1$samples)
#gelman.plot(samps.m1$samples)

out3_1 <- rep(NA, 3)
out3_1[1] <- samps.m1$dic[[1]] # mean deviance
out3_1[2] <- samps.m1$dic[[2]] # pD
out3_1[3] <- samps.m1$dic[[1]] + samps.m1$dic[[2]]  # DIC
write.table(out3_1, "../outputs/out3_1.txt", sep="\t")


################# model 2 a) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_s") 
jags.m2a <- jags.model(file="../models/MRFIT_m2a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2a, n.iter=n.burnin) # burn in
samps.m2a <- coda.samples.dic(jags.m2a, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2a$samples)[[1]], "out1_2a.txt", sep="\t")
write.table(summary(samps.m2a$samples)[[2]], "out2_2a.txt", sep="\t")
out3_2a <- rep(NA, 3)
out3_2a[1] <- samps.m2a$dic[[1]] # mean deviance
out3_2a[2] <- samps.m2a$dic[[2]] # pD
out3_2a[3] <- samps.m2a$dic[[1]] + samps.m2a$dic[[2]]  # DIC
write.table(out3_2a, "../outputs/out3_2a.txt", sep="\t")


################# Model 2 b) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_b") 
jags.m2b <- jags.model(file="../models/MRFIT_m2b.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2b, n.iter=n.burnin) # burn in
samps.m2b <- coda.samples.dic(jags.m2b, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2b$samples)[[1]], "out1_2b.txt", sep="\t")
write.table(summary(samps.m2b$samples)[[2]], "out2_2b.txt", sep="\t")

out3_2b <- rep(NA, 3)
out3_2b[1] <- samps.m2b$dic[[1]] # mean deviance
out3_2b[2] <- samps.m2b$dic[[2]] # pD
out3_2b[3] <- samps.m2b$dic[[1]] + samps.m2b$dic[[2]]  # DIC
write.table(out3_2b, "../outputs/out3_2b.txt", sep="\t")


################# Model 2 c) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_u") 
jags.m2c <- jags.model(file="../models/MRFIT_m2c.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2c, n.iter=n.burnin) # burn in
samps.m2c <- coda.samples.dic(jags.m2c, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2c$samples)[[1]], "out1_2c.txt", sep="\t")
write.table(summary(samps.m2c$samples)[[2]], "out2_2c.txt", sep="\t")

out3_2c <- rep(NA, 3)
out3_2c[1] <- samps.m2c$dic[[1]] # mean deviance
out3_2c[2] <- samps.m2c$dic[[2]] # pD
out3_2c[3] <- samps.m2c$dic[[1]] + samps.m2c$dic[[2]]  # DIC
write.table(out3_2c, "../outputs/out3_2c.txt", sep="\t")


################# Model 2 d) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_v") 
jags.m2d <- jags.model(file="../models/MRFIT_m2d.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2d, n.iter=n.burnin) # burn in
samps.m2d <- coda.samples.dic(jags.m2d, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2d$samples)[[1]], "out1_2d.txt", sep="\t")
write.table(summary(samps.m2d$samples)[[2]], "out2_2d.txt", sep="\t")

out3_2d <- rep(NA, 3)
out3_2d[1] <- samps.m2d$dic[[1]] # mean deviance
out3_2d[2] <- samps.m2d$dic[[2]] # pD
out3_2d[3] <- samps.m2d$dic[[1]] + samps.m2d$dic[[2]]  # DIC
write.table(out3_2d, "../outputs/out3_2d.txt", sep="\t")


################# Model 2 e) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "alpha_z1", "alpha_z2", "sigma_z1") 
jags.m2e <- jags.model(file="../models/MRFIT_m2e.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2e, n.iter=n.burnin) # burn in
samps.m2e <- coda.samples.dic(jags.m2e, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2e$samples)[[1]], "out1_2e.txt", sep="\t")
write.table(summary(samps.m2e$samples)[[2]], "out2_2e.txt", sep="\t")

out3_2e <- rep(NA, 3)
out3_2e[1] <- samps.m2e$dic[[1]] # mean deviance
out3_2e[2] <- samps.m2e$dic[[2]] # pD
out3_2e[3] <- samps.m2e$dic[[1]] + samps.m2e$dic[[2]]  # DIC
write.table(out3_2e, "../outputs/out3_2e.txt", sep="\t")


################# Model 2 f) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_z2") 
jags.m2f <- jags.model(file="../models/MRFIT_m2f.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2f, n.iter=n.burnin) # burn in
samps.m2f <- coda.samples.dic(jags.m2f, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2f$samples)[[1]], "out1_2f.txt", sep="\t")
write.table(summary(samps.m2f$samples)[[2]], "out2_2f.txt", sep="\t")

out3_2f <- rep(NA, 3)
out3_2f[1] <- samps.m2f$dic[[1]] # mean deviance
out3_2f[2] <- samps.m2f$dic[[2]] # pD
out3_2f[3] <- samps.m2f$dic[[1]] + samps.m2f$dic[[2]]  # DIC
write.table(out3_2f, "../outputs/out3_2f.txt", sep="\t")

################# Model 3 a) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_s") 
jags.m3a <- jags.model(file="../models/MRFIT_m3a.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3a, n.iter=n.burnin) # burn in
samps.m3a <- coda.samples.dic(jags.m3a, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m3a$samples)[[1]], "out1_3a.txt", sep="\t")
write.table(summary(samps.m3a$samples)[[2]], "out2_3a.txt", sep="\t")

out3_3a <- rep(NA, 3)
out3_3a[1] <- samps.m3a$dic[[1]] # mean deviance
out3_3a[2] <- samps.m3a$dic[[2]] # pD
out3_3a[3] <- samps.m3a$dic[[1]] + samps.m3a$dic[[2]]  # DIC
write.table(out3_3a, "../outputs/out3_3a.txt", sep="\t")


################# Model 3 b) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_b") 
jags.m3b <- jags.model(file="../models/MRFIT_m3b.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3b, n.iter=n.burnin) # burn in
samps.m3b <- coda.samples.dic(jags.m3b, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m3b$samples)[[1]], "out1_3b.txt", sep="\t")
write.table(summary(samps.m3b$samples)[[2]], "out2_3b.txt", sep="\t")

out3_3b <- rep(NA, 3)
out3_3b[1] <- samps.m3b$dic[[1]] # mean deviance
out3_3b[2] <- samps.m3b$dic[[2]] # pD
out3_3b[3] <- samps.m3b$dic[[1]] + samps.m3b$dic[[2]]  # DIC
write.table(out3_3b, "../outputs/out3_3b.txt", sep="\t")

################# Model 3 c) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_u") 
jags.m3c <- jags.model(file="../models/MRFIT_m3c.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3c, n.iter=n.burnin) # burn in
samps.m3c <- coda.samples.dic(jags.m3c, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m3c$samples)[[1]], "out1_3c.txt", sep="\t")
write.table(summary(samps.m3c$samples)[[2]], "out2_3c.txt", sep="\t")

out3_3c <- rep(NA, 3)
out3_3c[1] <- samps.m3c$dic[[1]] # mean deviance
out3_3c[2] <- samps.m3c$dic[[2]] # pD
out3_3c[3] <- samps.m3c$dic[[1]] + samps.m3c$dic[[2]]  # DIC
write.table(out3_3c, "../outputs/out3_3c.txt", sep="\t")


################# Model 3 d) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_v") 
jags.m3d <- jags.model(file="../models/MRFIT_m3d.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3d, n.iter=n.burnin) # burn in
samps.m3d <- coda.samples.dic(jags.m3d, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m3d$samples)[[1]], "out1_3d.txt", sep="\t")
write.table(summary(samps.m3d$samples)[[2]], "out2_3d.txt", sep="\t")

out3_3d <- rep(NA, 3)
out3_3d[1] <- samps.m3d$dic[[1]] # mean deviance
out3_3d[2] <- samps.m3d$dic[[2]] # pD
out3_3d[3] <- samps.m3d$dic[[1]] + samps.m3d$dic[[2]]  # DIC
write.table(out3_3d, "../outputs/out3_3d.txt", sep="\t")


################# Model 3 e) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "sigma_z1", "sigma_z2")  
jags.m3e <- jags.model(file="../models/MRFIT_m3e.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3e, n.iter=n.burnin) # burn in
samps.m3e <- coda.samples.dic(jags.m3e, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m3e$samples)[[1]], "out1_3e.txt", sep="\t")
write.table(summary(samps.m3e$samples)[[2]], "out2_3e.txt", sep="\t")

out3_3e <- rep(NA, 3)
out3_3e[1] <- samps.m3e$dic[[1]] # mean deviance
out3_3e[2] <- samps.m3e$dic[[2]] # pD
out3_3e[3] <- samps.m3e$dic[[1]] + samps.m3e$dic[[2]]  # DIC
write.table(out3_3e, "../outputs/out3_3e.txt", sep="\t")


################# Model 3 f) #########################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "Sigma_z")  
jags.m3f <- jags.model(file="../models/MRFIT_m3f.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m3f, n.iter=n.burnin) # burn in
samps.m3f <- coda.samples.dic(jags.m3f, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m3f$samples)[[1]], "out1_3f.txt", sep="\t")
write.table(summary(samps.m3f$samples)[[2]], "out2_3f.txt", sep="\t")

out3_3f <- rep(NA, 3)
out3_3f[1] <- samps.m3f$dic[[1]] # mean deviance
out3_3f[2] <- samps.m3f$dic[[2]] # pD
out3_3f[3] <- samps.m3f$dic[[1]] + samps.m3f$dic[[2]]  # DIC
write.table(out3_3f, "../outputs/out3_3f.txt", sep="\t")


################# Model 2 e) Sens 1 ####################################################
params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
            "pic", "pin", "pia", "alpha_z1", "alphz_z2", "sigma_z1") 
jags.m2e_s1 <- jags.model(file="../models/MRFIT_m2e_s1.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2e_s1, n.iter=n.burnin) # burn in
samps.m2e_s1 <- coda.samples.dic(jags.m2e_s1, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2e_s1$samples)[[1]], "out1_2e_s1.txt", sep="\t")
write.table(summary(samps.m2e_s1$samples)[[2]], "out2_2e_s1.txt", sep="\t")

out3_2e_s1 <- rep(NA, 3)
out3_2e_s1[1] <- samps.m2e_s1$dic[[1]] # mean deviance
out3_2e_s1[2] <- samps.m2e_s1$dic[[2]] # pD
out3_2e_s1[3] <- samps.m2e_s1$dic[[1]] + samps.m2e_s1$dic[[2]]  # DIC
write.table(out3_2e_s1, "../outputs/out3_2e_s1.txt", sep="\t")

################# Model 2 e) Sens 2 ######################################################
jags.m2e_s2 <- jags.model(file="../models/MRFIT_m2e_s2.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m2e_s2, n.iter=n.burnin) # burn in
samps.m2e_s2 <- coda.samples.dic(jags.m2e_s2, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m2e_s2$samples)[[1]], "out1_2e_s2.txt", sep="\t")
write.table(summary(samps.m2e_s2$samples)[[2]], "out2_2e_s2.txt", sep="\t")

out3_2e_s2 <- rep(NA, 3)
out3_2e_s2[1] <- samps.m2e_s2$dic[[1]] # mean deviance
out3_2e_s2[2] <- samps.m2e_s2$dic[[2]] # pD
out3_2e_s2[3] <- samps.m2e_s2$dic[[1]] + samps.m2e_s2$dic[[2]]  # DIC
write.table(out3_2e_s2, "../outputs/out3_2e_s2.txt", sep="\t")



############## save out3 tables ###################
setwd('../outputs/')
out3_1 <- read.table('out3_1.txt')
out3_2a <- read.table('out3_2a.txt')
out3_2b <- read.table('out3_2b.txt')
out3_2c <- read.table('out3_2c.txt')
out3_2d <- read.table('out3_2d.txt')
out3_2e <- read.table('out3_2e.txt')
out3_2f <- read.table('out3_2f.txt')
out3_3a <- read.table('out3_3a.txt')
out3_3b <- read.table('out3_3b.txt')
out3_3c <- read.table('out3_3c.txt')
out3_3d <- read.table('out3_3d.txt')
out3_3e <- read.table('out3_3e.txt')
out3_3f <- read.table('out3_3f.txt')
out3_2e_s1 <- read.table('out3_2e_s1.txt')
out3_2e_s2 <- read.table('out3_2e_s2.txt')

out3_12 <- cbind(out3_1, out3_2a, out3_2b, out3_2c, out3_2d, out3_2e, out3_2f)
out3_3 <- cbind(out3_3a, out3_3b, out3_3c, out3_3d, out3_3e, out3_3f)
out3_sens <- cbind(out3_2e, out3_2e_s1, out3_2e_s2)
write.table(out3_12, "out3_m12.txt", sep="\t")
write.table(out3_3, "out3_m3.txt", sep="\t")
write.table(out3_sens, "out3_sens.txt", sep="\t")
