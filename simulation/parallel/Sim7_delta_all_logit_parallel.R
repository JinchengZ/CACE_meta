rm(list=ls())
library(MASS)
#install.packages("mvtnorm")
library(mvtnorm)
library(coda)
library(rjags)
library(VGAM)
library(parallel)

source('../fun_coda_dic.R')

I <- 20  # number of studies in each set
R <- 2000 # number of simulation, should be 1000

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
sigma_z1 <- sigma_z2 <- sigma_u <- sigma_v <-  sigma_s <- sigma_b <- 0.5

# assign empty
delta_z1 <- delta_z2 <- delta_u <- delta_v <-  delta_s <- delta_b <- pi_c <- pi_n <- pi_a <- rep(NA,I)
z1 <- z2 <- u1 <- v1 <- s1 <- b1 <- rr <- rep(NA, I)
prob <- RR <- matrix(NA, nrow = I, ncol = 8)
prob0 <- prob1 <- R0_4 <- R1_4 <- matrix(NA, nrow = I, ncol = 4)
delta_Z <- matrix(NA, nrow = I, ncol = 2)

# assign empty matrix to store final estimates
sim7_pD <- sim7_DIC <- rep(NA, 9)
sim7_CACE <- rep(NA, 56)

seed <- 1234
n.chains <- 3
n.adapt <- 1000
n.burnin <- 10000
n.iter <- 100000
n.thin <- 1

wrap=function(j)
{
  r <- rep(0.5, I)  # prob of Randomized trt=1
  Ntol <- rep(350, I) # number of subjects in each study
  N0 <- N1 <- rep(175, I)  
  # Generate parameters and data, 5) random all
  set.seed(123+j)
  for (i in 1:I) {
    set.seed(123+j*R+i)
    delta_Z[i, 1:2] <- mvrnorm(n = 1, mu = c(0, 0), 
                       Sigma = matrix(c(sigma_z1^2, 0.5*sigma_z1*sigma_z2,
                                        0.5*sigma_z1*sigma_z2, sigma_z2^2), 2, 2))
    z1[i] <- alpha_z1 + delta_Z[i, 1]
    z2[i] <- alpha_z2 + delta_Z[i, 2]
    # delta_z1[i] <- rnorm(1, 0, sigma_z1)
    # delta_z2[i] <- rnorm(1, 0, sigma_z2)
    # z1[i] <- alpha_z1 + delta_z1[i]
    # z2[i] <- alpha_z2 + delta_z2[i]
    
    pi_n[i] = exp(z1[i])/(1+exp(z1[i])+exp(z2[i]))
    pi_a[i] = exp(z2[i])/(1+exp(z1[i])+exp(z2[i]))
    pi_c[i] = 1-pi_a[i]-pi_n[i]
    
    set.seed(1234+j*R+i)
    delta_u[i] <- rnorm(1, 0, sigma_u)
    delta_v[i] <- rnorm(1, 0, sigma_v)
    delta_s[i] <- rnorm(1, 0, sigma_s)
    delta_b[i] <- rnorm(1, 0, sigma_b)
    
    u1[i] <- probit(alpha_u + delta_u[i], inverse=T) 
    v1[i] <- probit(alpha_v + delta_v[i], inverse=T)
    s1[i] <- logit(alpha_s + delta_s[i], inverse=T) 
    b1[i] <- logit(alpha_b + delta_b[i], inverse=T) 
    
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
  
  sim7_CACE[1] <- probit(alpha_u/sqrt(1+sigma_u^2), inverse=T) - probit(alpha_v/sqrt(1+sigma_v^2), inverse=T)
  sim7_CACE[56] <- mean(u1-v1)
  
  #data = list(Ntol=Ntol, R=RR, I=I, r=rr)
  RR = cbind(R0_4, R1_4)
  Csq <- (16*16*3)/(15*15*pi*pi)
  data=list(N0=N0, N1=N1, R0_4=R0_4, R1_4=R1_4, I=I, Csq=Csq)
  
  write.table(RR, paste("../simdata/data7_", j,".txt",sep=""), sep="\t")

  
  ### set initials ##########################
  set.seed(seed+j)
  init.seeds <- sample(1:1000000, n.chains)
  
  Inits <- vector("list", n.chains)
  for(i in 1:n.chains){
    Inits[[i]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = init.seeds[i])
  }
  
  ############# Model 1. none ############################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia")   
  
  jags.m7_1 <- jags.model(file="../models/sim_logit_m1.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_1, n.iter=n.burnin) # burn in
  samps.m7_1 <- coda.samples.dic(jags.m7_1, variable.names=params, n.iter=n.iter, thin=n.thin)
  out.m7_1 <- dic.samples(jags.m7_1, 10000, type='pD')
  
  sim7_CACE[2] <- summary(samps.m7_1$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[3:5] <- c(summary(samps.m7_1$samples)[[2]][1,1], 
                      summary(samps.m7_1$samples)[[2]][1,3], 
                      summary(samps.m7_1$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[6] <- ifelse((sim7_CACE[3]<=sim7_CACE[1]) & (sim7_CACE[5]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[7] <- sim7_CACE[5] - sim7_CACE[3]  # CACE length
  
  sim7_pD[1] <- samps.m7_1$dic[[2]]
  sim7_DIC[1] <- samps.m7_1$dic[[1]] + samps.m7_1$dic[[2]]
  
  
  ############# Model 2. delta_z1 ##########################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_z1") 
  
  jags.m7_2 <- jags.model(file="../models/sim_logit_m2.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_2, n.iter=n.burnin) # burn in
  samps.m7_2 <- coda.samples.dic(jags.m7_2, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[8] <- summary(samps.m7_2$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[9:11] <- c(summary(samps.m7_2$samples)[[2]][1,1], 
                       summary(samps.m7_2$samples)[[2]][1,3], 
                       summary(samps.m7_2$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[12] <- ifelse((sim7_CACE[9]<=sim7_CACE[1]) & (sim7_CACE[11]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[13] <- sim7_CACE[11] - sim7_CACE[9]   # CACE length
  
  sim7_pD[2] <- samps.m7_2$dic[[2]]
  sim7_DIC[2] <- samps.m7_2$dic[[1]] + samps.m7_2$dic[[2]]
  
  
  ############# Model 3. delta_u  ########################################################  
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_u")  
  
  jags.m7_3 <- jags.model(file="../models/sim_logit_m3.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_3, n.iter=n.burnin) # burn in
  samps.m7_3 <- coda.samples.dic(jags.m7_3, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[14] <- summary(samps.m7_3$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[15:17] <- c(summary(samps.m7_3$samples)[[2]][1,1], 
                        summary(samps.m7_3$samples)[[2]][1,3], 
                        summary(samps.m7_3$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[18] <- ifelse((sim7_CACE[15]<=sim7_CACE[1]) & (sim7_CACE[17]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[19] <- sim7_CACE[17] - sim7_CACE[15]  # CACE length
  
  sim7_pD[3] <- samps.m7_3$dic[[2]]
  sim7_DIC[3] <- samps.m7_3$dic[[1]] + samps.m7_3$dic[[2]]  
  
  
  ############# Model 4. delta_v ######################################################### 
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", "pic", "pin", "pia", "sigma_v")  
  
  jags.m7_4 <- jags.model(file="../models/sim_logit_m4.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_4, n.iter=n.burnin) # burn in
  samps.m7_4 <- coda.samples.dic(jags.m7_4, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[20] <- summary(samps.m7_4$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[21:23] <- c(summary(samps.m7_4$samples)[[2]][1,1], 
                        summary(samps.m7_4$samples)[[2]][1,3], 
                        summary(samps.m7_4$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[24] <- ifelse((sim7_CACE[21]<=sim7_CACE[1]) & (sim7_CACE[23]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[25] <- sim7_CACE[23] - sim7_CACE[21]  # CACE length
  
  sim7_pD[4] <- samps.m7_4$dic[[2]]
  sim7_DIC[4] <- samps.m7_4$dic[[1]] + samps.m7_4$dic[[2]] 
  
  
  ############# Model 5. delta_z1, delta_u ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_z1", "sigma_u")  
  
  jags.m7_5 <- jags.model(file="../models/sim_logit_m5.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_5, n.iter=n.burnin) # burn in
  samps.m7_5 <- coda.samples.dic(jags.m7_5, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[26] <- summary(samps.m7_5$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[27:29] <- c(summary(samps.m7_5$samples)[[2]][1,1], 
                        summary(samps.m7_5$samples)[[2]][1,3], 
                        summary(samps.m7_5$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[30] <- ifelse((sim7_CACE[27]<=sim7_CACE[1]) & (sim7_CACE[29]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[31] <- sim7_CACE[29] - sim7_CACE[27]  # CACE length
  
  sim7_pD[5] <- samps.m7_5$dic[[2]]
  sim7_DIC[5] <- samps.m7_5$dic[[1]] + samps.m7_5$dic[[2]]
  
  
  ############# Model 6. delta_z1, delta_v ###################################################   
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_z1", "sigma_v") 
  
  jags.m7_6 <- jags.model(file="../models/sim_logit_m6.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_6, n.iter=n.burnin) # burn in
  samps.m7_6 <- coda.samples.dic(jags.m7_6, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[32] <- summary(samps.m7_6$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[33:35] <- c(summary(samps.m7_6$samples)[[2]][1,1], 
                        summary(samps.m7_6$samples)[[2]][1,3], 
                        summary(samps.m7_6$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[36] <- ifelse((sim7_CACE[33]<=sim7_CACE[1]) & (sim7_CACE[35]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[37] <- sim7_CACE[35] - sim7_CACE[33]  # CACE length
  
  sim7_pD[6] <- samps.m7_6$dic[[2]]
  sim7_DIC[6] <- samps.m7_6$dic[[1]] + samps.m7_6$dic[[2]]  
  
  
  ############# Model 7. delta_u, delta_v ################################################### 
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_u", "sigma_v") 
  
  jags.m7_7 <- jags.model(file="../models/sim_logit_m7.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_7, n.iter=n.burnin) # burn in
  samps.m7_7 <- coda.samples.dic(jags.m7_7, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[38] <- summary(samps.m7_7$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[39:41] <- c(summary(samps.m7_7$samples)[[2]][1,1], 
                        summary(samps.m7_7$samples)[[2]][1,3], 
                        summary(samps.m7_7$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[42] <- ifelse((sim7_CACE[39]<=sim7_CACE[1]) & (sim7_CACE[41]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[43] <- sim7_CACE[41] - sim7_CACE[39]  # CACE length
  
  sim7_pD[7] <- samps.m7_7$dic[[2]]
  sim7_DIC[7] <- samps.m7_7$dic[[1]] + samps.m7_7$dic[[2]]   
  
  
  ############# Model 8. delta_z1, delta_u, delta_v ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "sigma_z1", "sigma_u", "sigma_v") 
  
  jags.m7_8 <- jags.model(file="../models/sim_logit_m8.txt", data=data, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_8, n.iter=n.burnin) # burn in
  samps.m7_8 <- coda.samples.dic(jags.m7_8, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[44] <- summary(samps.m7_8$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[45:47] <- c(summary(samps.m7_8$samples)[[2]][1,1], 
                        summary(samps.m7_8$samples)[[2]][1,3], 
                        summary(samps.m7_8$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[48] <- ifelse((sim7_CACE[45]<=sim7_CACE[1]) & (sim7_CACE[47]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[49] <- sim7_CACE[47] - sim7_CACE[45]  # CACE length
  
  sim7_pD[8] <- samps.m7_8$dic[[2]]
  sim7_DIC[8] <- samps.m7_8$dic[[1]] + samps.m7_8$dic[[2]]  
  
  ############# Model 9. all deltas ###################################################
  params <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
              "pic", "pin", "pia", "Sigma_Z", "sigma_u", "sigma_v",
              "sigma_s", "sigma_b") 
  data_9 <- list(N0=N0, N1=N1, R0_4=R0_4, R1_4=R1_4, I=I, Csq=Csq,
                 II = structure(.Data = c(1, 0, 0, 1), .Dim = c(2, 2))  )
  jags.m7_9 <- jags.model(file="../models/sim_logit_m9.txt", data=data_9, inits=Inits, n.chains=n.chains, n.adapt=n.adapt)
  update(jags.m7_9, n.iter=n.burnin) # burn in
  samps.m7_9 <- coda.samples.dic(jags.m7_9, variable.names=params, n.iter=n.iter, thin=n.thin)
  
  sim7_CACE[44+6] <- summary(samps.m7_9$samples)[[1]][1,1]  # mean CACE
  sim7_CACE[45+6:47+6] <- c(summary(samps.m7_9$samples)[[2]][1,1], 
                        summary(samps.m7_9$samples)[[2]][1,3], 
                        summary(samps.m7_9$samples)[[2]][1,5])  # 2.5% 50% 97.5% CACE
  sim7_CACE[48+6] <- ifelse((sim7_CACE[45+6]<=sim7_CACE[1]) & (sim7_CACE[47+6]>=sim7_CACE[1]), 1, 0)
  sim7_CACE[49+6] <- sim7_CACE[47+6] - sim7_CACE[45+6]  # CACE length
  
  sim7_pD[9] <- samps.m7_9$dic[[2]]
  sim7_DIC[9] <- samps.m7_9$dic[[1]] + samps.m7_9$dic[[2]]  
  
  write.table(sim7_CACE, paste("../sim_CACE/sim7_CACE_", j, ".txt", sep=""), sep="\t")
  write.table(sim7_pD, paste("../sim_pD/sim7_pD_", j, ".txt", sep=""), sep="\t")
  write.table(sim7_DIC, paste("../sim_DIC/sim7_DIC_", j, ".txt", sep=""), sep="\t")
  
  output <- list(sim7_CACE, sim7_pD, sim7_DIC)
  return(output)
}

runsim <- mclapply(1:R, wrap, mc.cores=24)