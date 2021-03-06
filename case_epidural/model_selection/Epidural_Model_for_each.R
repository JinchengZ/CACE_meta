rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)
library(R2jags)

source('../fun_coda_dic.R')
obs <- read.table("../epidural.txt",header=T)
attach(obs)

Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
N0 <- n000+n001+n010+n011
N1 <- n100+n101+n110+n111
is.vector(Ntol)
R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
n <- length(Ntol)
r <- (n100+n101+n110+n111)/Ntol


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


params <- c("CACE", "alpha_n", "alpha_a", "alpha_s1", "alpha_b1", 
            "alpha_u", "alpha_v", "pi_c", "pi_n", "pi_a") 

################# rjags code #########################################################

for (i in 1:n) {
data <- list(N0=N0[i], N1=N1[i], R=R[i,], r=r[i])

jags.m_foreach <- jags.model(file="../models/epidural_m_foreach.txt", data=data, inits=init.jags, n.chains=n.chains, n.adapt=n.adapt)
update(jags.m_foreach, n.iter=n.burnin) # burn in
samps.m_foreach <- coda.samples.dic(jags.m_foreach, variable.names=params, n.iter=n.iter, thin=n.thin)

write.table(summary(samps.m_foreach$samples)[[1]], paste("../outputs/outeach1_",i,".txt",sep=""),sep="\t")
write.table(summary(samps.m_foreach$samples)[[2]], paste("../outputs/outeach2_",i,".txt",sep=""),sep="\t")

}



# put outcomes for the 10 studies together

foreach_CACE <- matrix(NA, n, 5)

for (i in 1:n){
  foreach_CACE[i, ] <- cbind(read.table(paste("../outputs/outeach2_",i,".txt",sep=""))[1, 1], # lower
                             read.table(paste("../outputs/outeach1_",i,".txt",sep=""))[1, 1], # mean
                             read.table(paste("../outputs/outeach2_",i,".txt",sep=""))[1, 3], # median
                             read.table(paste("../outputs/outeach2_",i,".txt",sep=""))[1, 5], # upper
                             read.table(paste("../outputs/outeach1_",i,".txt",sep=""))[1, 2]) # SD
}
data <- as.data.frame(foreach_CACE)
colnames(data)<-c("Lower", "Mean", "Median", "Upper", "SD")

# Save the data
write.table(data, "../outputs/CACE_foreach.txt", sep="\t", row.names = F)

meta1 <- rma(yi=Mean, sei=SD, data=data, method = "REML")
meta2 <- rma(yi=Mean, sei=SD, data=data, method = "DL")
meta3 <- rma(yi=Mean, sei=SD, data=data, method = "FE")

a <- cbind(meta1$beta, meta1$ci.lb, meta1$ci.ub)
out1 <- as.data.frame(a)
colnames(out1)<-c("Mean", "Lower", "Upper")
rownames(out1)<-c("CACE")
write.table(out1, "../outputs/CACEout.txt", sep="\t", row.names = F)
