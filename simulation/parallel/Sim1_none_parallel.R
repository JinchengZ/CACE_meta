rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)
library(VGAM)
library(parallel)

source('../fun_coda_dic.R')

I <- 20  # number of studies in each set
R <- 2000 #1000 # number of simulation, should be 1000

# assign true values
alpha_z1 <- -0.4
alpha_z2 <- -0.6
alpha_s <- 0.5
alpha_b <- -0.5
alpha_u <- -0.5
alpha_v <- 0.5
# alpha_z1 <- 0.5
# alpha_z2 <- -0.3
# alpha_s <- 0.5
# alpha_b <- -1
# alpha_u <- -1
# alpha_v <- 0.5
sigma_z1 <- sigma_u <- sigma_v <- 0.5

# assign empty
delta_z1 <- delta_u <- delta_v <- pi_c <- pi_n <- pi_a <- rep(NA,I)
z1 <- z2 <- u1 <- v1 <- s1 <- b1 <- rr <- rep(NA, I)
prob <- RR <- matrix(NA, nrow = I, ncol = 8)
prob0 <- prob1 <- R0_4 <- R1_4 <- matrix(NA, nrow = I, ncol = 4)

# assign empty matrix to store final estimates
sim1_pD <- sim1_DIC <- rep(NA, 8)
sim1_CACE <- rep(NA, 50)

seed <- 1234
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 1

wrap <- function(j)
{
  r <- rep(0.5, I)  # prob of Randomized trt=1
  Ntol <- rep(350, I) # number of subjects in each study
  N0 <- N1 <- rep(175, I)
  # Generate parameters and data, 1) no random
  set.seed(seed+j)
  for (i in 1:I) {
    set.seed(seed+j*R+i)
    
    z1[i] <- alpha_z1 
    z2[i] <- alpha_z2 
    
    pi_n[i] = exp(z1[i])/(1+exp(z1[i])+exp(z2[i]))
    pi_a[i] = exp(z2[i])/(1+exp(z1[i])+exp(z2[i]))
    pi_c[i] = 1-pi_a[i]-pi_n[i]
    
    u1[i] <- probit(alpha_u, inverse=T) 
    v1[i] <- probit(alpha_v, inverse=T)
    s1[i] <- logit(alpha_s, inverse=T) 
    b1[i] <- logit(alpha_b, inverse=T) 
    
    # generate multinomial data
    # prob[i, ] <- c((1-r[i])*(pi_n[i]*(1-s1[i]) + pi_c[i]*(1-v1[i])), 
    #                (1-r[i])*(pi_n[i]*s1[i] + pi_c[i]*v1[i]), 
    #                (1-r[i])*(pi_a[i]*(1-b1[i])), 
    #                (1-r[i])*(pi_a[i]*b1[i]), 
    #                r[i]*(pi_n[i]*(1-s1[i])), 
    #                r[i]*(pi_n[i]*s1[i]), 
    #                r[i]*(pi_c[i]*(1-u1[i])+pi_a[i]*(1-b1[i])), 
    #                r[i]*(pi_c[i]*u1[i]+pi_a[i]*b1[i])
    # )
    # RR[i, ] <- t(rmultinom(1, Ntol[i], prob[i, ]))
    # rr[i] <- (RR[i, 5]+RR[i, 6]+RR[i, 7]+RR[i, 8])/Ntol[i]
    
    prob0[i, ] <- c((pi_n[i]*(1-s1[i]) + pi_c[i]*(1-v1[i])), 
                    (pi_n[i]*s1[i] + pi_c[i]*v1[i]), 
                    (pi_a[i]*(1-b1[i])), 
                    (pi_a[i]*b1[i])    )
    prob1[i, ] <- c((pi_n[i]*(1-s1[i])), 
                    (pi_n[i]*s1[i]), 
                    (pi_c[i]*(1-u1[i])+pi_a[i]*(1-b1[i])), 
                    (pi_c[i]*u1[i]+pi_a[i]*b1[i])    )
    
    R0_4[i, ] <- t(rmultinom(1, N0[i], prob0[i, ]))
    R1_4[i, ] <- t(rmultinom(1, N1[i], prob1[i, ]))
  }
  
  sim1_CACE[1] <- probit(alpha_u, inverse=T) - probit(alpha_v, inverse=T)
  sim1_CACE[50] <- mean(u1-v1)
  
  #data = list(Ntol=Ntol, R=RR, I=I, r=rr)
  RR = cbind(R0_4, R1_4)
  data=list(N0=N0, N1=N1, R0_4=R0_4, R1_4=R1_4, I=I)
  
  write.table(RR, paste("../simdata/data1_", j,".txt",sep=""), sep="\t")
  
  ### set initials ##########################
  set.seed(seed)
  init.seeds <- sample(1:1000000, n.chains)
  
  Inits <- vector("list", n.chains)
  for(i in 1:n.chains){
    Inits[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i]+j)
  }
  
  ############# Model 1. none ############################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")  
  
  jags.m1_1 <- jags.model(file="../models/sim_m1.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_1, n.iter=n.burnin) # burn in
  samps.m1_1 <- coda.samples.dic(jags.m1_1, variable.names=params, n.iter=n.iter, thin=n.thin)
  out.m1_1 <- dic.samples(jags.m1_1, 10000, type='pD')
  
  sim1_CACE[2] <- summary(samps.m1_1$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[3:5] <- c(summary(samps.m1_1$samples)[[2]][1,1], 
                      summary(samps.m1_1$samples)[[2]][1,3], 
                      summary(samps.m1_1$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[6] <- ifelse((sim1_CACE[3]<=sim1_CACE[1]) & (sim1_CACE[5]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[7] <- sim1_CACE[5] - sim1_CACE[3]  # CACE length
  
  sim1_pD[1] <- samps.m1_1$dic[[2]]
  sim1_DIC[1] <- samps.m1_1$dic[[1]] + samps.m1_1$dic[[2]]
  
  
  ############# Model 2. delta_z1 ##########################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_z1") 
  
  jags.m1_2 <- jags.model(file="../models/sim_m2.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_2, n.iter=n.burnin) # burn in
  samps.m1_2 <- coda.samples.dic(jags.m1_2, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[8] <- summary(samps.m1_2$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[9:11] <- c(summary(samps.m1_2$samples)[[2]][1,1], 
                       summary(samps.m1_2$samples)[[2]][1,3], 
                       summary(samps.m1_2$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[12] <- ifelse((sim1_CACE[9]<=sim1_CACE[1]) & (sim1_CACE[11]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[13] <- sim1_CACE[11] - sim1_CACE[9]   # CACE length
  
  sim1_pD[2] <- samps.m1_2$dic[[2]]
  sim1_DIC[2] <- samps.m1_2$dic[[1]] + samps.m1_2$dic[[2]]
  
  
  ############# Model 3. delta_u  ########################################################  
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_u")  
  
  jags.m1_3 <- jags.model(file="../models/sim_m3.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_3, n.iter=n.burnin) # burn in
  samps.m1_3 <- coda.samples.dic(jags.m1_3, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[14] <- summary(samps.m1_3$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[15:17] <- c(summary(samps.m1_3$samples)[[2]][1,1], 
                        summary(samps.m1_3$samples)[[2]][1,3], 
                        summary(samps.m1_3$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[18] <- ifelse((sim1_CACE[15]<=sim1_CACE[1]) & (sim1_CACE[17]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[19] <- sim1_CACE[17] - sim1_CACE[15]  # CACE length
  
  sim1_pD[3] <- samps.m1_3$dic[[2]]
  sim1_DIC[3] <- samps.m1_3$dic[[1]] + samps.m1_3$dic[[2]]  
  
  
  ############# Model 4. delta_v ######################################################### 
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_v")  
  
  jags.m1_4 <- jags.model(file="../models/sim_m4.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_4, n.iter=n.burnin) # burn in
  samps.m1_4 <- coda.samples.dic(jags.m1_4, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[20] <- summary(samps.m1_4$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[21:23] <- c(summary(samps.m1_4$samples)[[2]][1,1], 
                        summary(samps.m1_4$samples)[[2]][1,3], 
                        summary(samps.m1_4$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[24] <- ifelse((sim1_CACE[21]<=sim1_CACE[1]) & (sim1_CACE[23]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[25] <- sim1_CACE[23] - sim1_CACE[21]  # CACE length
  
  sim1_pD[4] <- samps.m1_4$dic[[2]]
  sim1_DIC[4] <- samps.m1_4$dic[[1]] + samps.m1_4$dic[[2]] 
  
  
  ############# Model 5. delta_z1, delta_u ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_z1", "sigma_u")  
  
  jags.m1_5 <- jags.model(file="../models/sim_m5.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_5, n.iter=n.burnin) # burn in
  samps.m1_5 <- coda.samples.dic(jags.m1_5, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[26] <- summary(samps.m1_5$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[27:29] <- c(summary(samps.m1_5$samples)[[2]][1,1], 
                        summary(samps.m1_5$samples)[[2]][1,3], 
                        summary(samps.m1_5$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[30] <- ifelse((sim1_CACE[27]<=sim1_CACE[1]) & (sim1_CACE[29]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[31] <- sim1_CACE[29] - sim1_CACE[27]  # CACE length
  
  sim1_pD[5] <- samps.m1_5$dic[[2]]
  sim1_DIC[5] <- samps.m1_5$dic[[1]] + samps.m1_5$dic[[2]]
  
  
  ############# Model 6. delta_z1, delta_v ###################################################   
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_z1", "sigma_v") 
  
  jags.m1_6 <- jags.model(file="../models/sim_m6.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_6, n.iter=n.burnin) # burn in
  samps.m1_6 <- coda.samples.dic(jags.m1_6, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[32] <- summary(samps.m1_6$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[33:35] <- c(summary(samps.m1_6$samples)[[2]][1,1], 
                        summary(samps.m1_6$samples)[[2]][1,3], 
                        summary(samps.m1_6$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[36] <- ifelse((sim1_CACE[33]<=sim1_CACE[1]) & (sim1_CACE[35]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[37] <- sim1_CACE[35] - sim1_CACE[33]  # CACE length
  
  sim1_pD[6] <- samps.m1_6$dic[[2]]
  sim1_DIC[6] <- samps.m1_6$dic[[1]] + samps.m1_6$dic[[2]]  
  
  
  ############# Model 7. delta_u, delta_v ################################################### 
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_u", "sigma_v") 
  
  jags.m1_7 <- jags.model(file="../models/sim_m7.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_7, n.iter=n.burnin) # burn in
  samps.m1_7 <- coda.samples.dic(jags.m1_7, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[38] <- summary(samps.m1_7$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[39:41] <- c(summary(samps.m1_7$samples)[[2]][1,1], 
                        summary(samps.m1_7$samples)[[2]][1,3], 
                        summary(samps.m1_7$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[42] <- ifelse((sim1_CACE[39]<=sim1_CACE[1]) & (sim1_CACE[41]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[43] <- sim1_CACE[41] - sim1_CACE[39]  # CACE length
  
  sim1_pD[7] <- samps.m1_7$dic[[2]]
  sim1_DIC[7] <- samps.m1_7$dic[[1]] + samps.m1_7$dic[[2]]   
  
  
  ############# Model 8. delta_z1, delta_u, delta_v ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_z1", "sigma_u", "sigma_v") 
  
  jags.m1_8 <- jags.model(file="../models/sim_m8.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m1_8, n.iter=n.burnin) # burn in
  samps.m1_8 <- coda.samples.dic(jags.m1_8, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim1_CACE[44] <- summary(samps.m1_8$samples)[[1]][1,1]  # mean CACE
  sim1_CACE[45:47] <- c(summary(samps.m1_8$samples)[[2]][1,1], 
                        summary(samps.m1_8$samples)[[2]][1,3], 
                        summary(samps.m1_8$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim1_CACE[48] <- ifelse((sim1_CACE[45]<=sim1_CACE[1]) & (sim1_CACE[47]>=sim1_CACE[1]), 1, 0)
  sim1_CACE[49] <- sim1_CACE[47] - sim1_CACE[45]  # CACE length
  
  sim1_pD[8] <- samps.m1_8$dic[[2]]
  sim1_DIC[8] <- samps.m1_8$dic[[1]] + samps.m1_8$dic[[2]]  
  
  write.table(sim1_CACE, paste("../sim_CACE/sim1_CACE_", j, ".txt", sep=""), sep="\t")
  write.table(sim1_pD, paste("../sim_pD/sim1_pD_", j, ".txt", sep=""), sep="\t")
  write.table(sim1_DIC, paste("../sim_DIC/sim1_DIC_", j, ".txt", sep=""), sep="\t")
  
  output <- list(sim1_CACE, sim1_pD, sim1_DIC)
  return(output)
}

runsim <- mclapply(1:R, wrap, mc.cores=24)
