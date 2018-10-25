R<- 200
out1_CACE <- out2_CACE <- out3_CACE <- out4_CACE <- matrix(NA, R, 6)

library(VGAM)
# assign true values
alpha_z1 <- -0.4
alpha_z2 <- -0.6
alpha_s <- 0.5
alpha_b <- -0.5
alpha_u <- -0.5
alpha_v <- 0.5
sigma_z1 <- sigma_u <- sigma_v <- 0.5

true12 <- probit(alpha_u, inverse=T) - probit(alpha_v, inverse=T)
true34 <- probit(alpha_u/sqrt(1+sigma_u^2), inverse=T) - probit(alpha_v, inverse=T)

for (i in 1:R){
  setwd("../foreach")
  out1_CACE[i, 1] <- true12
  out2_CACE[i, 1] <- true12
  out3_CACE[i, 1] <- true34
  out4_CACE[i, 1] <- true34
  
  out1_CACE[i, 2:4] <- as.matrix(read.table(paste("CACEout1_",i,".txt",sep=""), header = T) )[1, ]
  out2_CACE[i, 2:4] <- as.matrix(read.table(paste("CACEout2_",i,".txt",sep=""), header = T) )[1, ]
  out3_CACE[i, 2:4] <- as.matrix(read.table(paste("CACEout3_",i,".txt",sep=""), header = T) )[1, ]
  out4_CACE[i, 2:4] <- as.matrix(read.table(paste("CACEout4_",i,".txt",sep=""), header = T) )[1, ]
  
  out1_CACE[i, 5] <- out1_CACE[i, 2] - out1_CACE[i, 1]
  out2_CACE[i, 5] <- out2_CACE[i, 2] - out2_CACE[i, 1]
  out3_CACE[i, 5] <- out3_CACE[i, 2] - out3_CACE[i, 1]
  out4_CACE[i, 5] <- out4_CACE[i, 2] - out4_CACE[i, 1]
  
  out1_CACE[i, 6] <- ifelse((out1_CACE[3]<=out1_CACE[1]) & (out1_CACE[4]>=out1_CACE[1]), 1, 0)
  out2_CACE[i, 6] <- ifelse((out2_CACE[3]<=out2_CACE[1]) & (out2_CACE[4]>=out2_CACE[1]), 1, 0)
  out3_CACE[i, 6] <- ifelse((out3_CACE[3]<=out3_CACE[1]) & (out3_CACE[4]>=out3_CACE[1]), 1, 0)
  out4_CACE[i, 6] <- ifelse((out4_CACE[3]<=out4_CACE[1]) & (out4_CACE[4]>=out4_CACE[1]), 1, 0)
}

CACE_out <- matrix(NA, 4, 6)
  setwd("../foreach")
  CACE_out <- rbind(colMeans(out1_CACE), colMeans(out2_CACE),  
                    colMeans(out3_CACE), colMeans(out4_CACE))


# Save the data
write.table(out1_CACE, "../outputs/out1_CACE.txt", sep="\t", row.names = FALSE)
write.table(out2_CACE, "../outputs/out2_CACE.txt", sep="\t", row.names = FALSE)
write.table(out3_CACE, "../outputs/out3_CACE.txt", sep="\t", row.names = FALSE)
write.table(out4_CACE, "../outputs/out4_CACE.txt", sep="\t", row.names = FALSE)

write.table(CACE_out, "../outputs/CACE_out.txt", sep="\t", row.names = FALSE)
