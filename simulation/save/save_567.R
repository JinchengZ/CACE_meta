R <- 2000
sim5_CACE <-sim6_CACE <- sim7_CACE <- matrix(NA, R, 98)
sim5_DIC <- sim6_DIC <- sim7_DIC <- matrix(NA, R, 9)
sim5_pD <- sim6_pD <- sim7_pD  <- matrix(NA, R, 9)
sim5_loc <- sim6_loc <- sim7_loc <- rep(NA, R)
for (i in 1:R){
  setwd("../sim_CACE")
  #sim5_CACE[i, ] <- read.table(paste("sim5_CACE_",i,".txt",sep=""))$x
  sim6_CACE[i, ] <- read.table(paste("sim6_CACE_",i,".txt",sep=""))$x
  sim6_CACE[i, 51] <- read.table(paste("sim6_CACE_",i,".txt",sep=""))$x[57]
  sim6_CACE[i, 52] <- read.table(paste("sim6_CACE_",i,".txt",sep=""))$x[58]
  sim6_CACE[i, 53] <- read.table(paste("sim6_CACE_",i,".txt",sep=""))$x[59]
  sim6_CACE[i, 54] <- ifelse((sim6_CACE[i, 51]<=sim6_CACE[i, 1]) & (sim6_CACE[i, 53]>=sim6_CACE[i, 1]), 1, 0)
  sim6_CACE[i, 55] <- sim6_CACE[i, 53] - sim6_CACE[i, 51]
  
  sim7_CACE[i, ] <- read.table(paste("sim7_CACE_",i,".txt",sep=""))$x  
  sim7_CACE[i, 51] <- read.table(paste("sim7_CACE_",i,".txt",sep=""))$x[57]
  sim7_CACE[i, 52] <- read.table(paste("sim7_CACE_",i,".txt",sep=""))$x[58]
  sim7_CACE[i, 53] <- read.table(paste("sim7_CACE_",i,".txt",sep=""))$x[59]
  sim7_CACE[i, 54] <- ifelse((sim7_CACE[i, 51]<=sim7_CACE[i, 1]) & (sim7_CACE[i, 53]>=sim7_CACE[i, 1]), 1, 0)
  sim7_CACE[i, 55] <- sim7_CACE[i, 53] - sim7_CACE[i, 51]

  setwd("../sim_DIC")
  #sim5_DIC[i, ] <- read.table(paste("sim5_DIC_",i,".txt",sep=""))$x
  #sim5_loc[i] <- which.min(sim5_DIC[i,])
  sim6_DIC[i, ] <- read.table(paste("sim6_DIC_",i,".txt",sep=""))$x
  sim6_loc[i] <- which.min(sim6_DIC[i,])
  sim7_DIC[i, ] <- read.table(paste("sim7_DIC_",i,".txt",sep=""))$x
  sim7_loc[i] <- which.min(sim7_DIC[i,])

  setwd("../sim_pD")
  #sim5_pD[i, ] <- read.table(paste("sim5_pD_",i,".txt",sep=""))$x
  sim6_pD[i, ] <- read.table(paste("sim6_pD_",i,".txt",sep=""))$x
  sim7_pD[i, ] <- read.table(paste("sim7_pD_",i,".txt",sep=""))$x
}

table1_5 <- table1_6 <- table1_7  <- matrix(NA, 1, 9)
for (i in 1:9) {
  #table1_5[i] <- length(which(sim5_loc == i)) 
  table1_6[i] <- length(which(sim6_loc == i)) 
  table1_7[i] <- length(which(sim7_loc == i)) 
}

table2_5 <- table2_6 <- table2_7 <- table2_4 <- matrix(NA, 7, 9)
for (i in 1:9) {
  #table2_5[1, i] <- mean(sim5_CACE[, 1])      # true value
  #table2_5[2, i] <- mean(sim5_CACE[, 6*i-4])  # mean
  #table2_5[3, i] <- mean(sim5_CACE[, 6*i-2])  # median
  #table2_5[4, i] <- mean(sim5_CACE[, 6*i-4]-sim5_CACE[, 1])  # bias, mean-true
  #table2_5[5, i] <- mean(sim5_CACE[, 6*i+1])  # CI length
  #table2_5[6, i] <- sum(sim5_CACE[, 6*i])/R   # Coverage prob
  #table2_5[7, i] <- mean(sim5_CACE[, 50])      # mean true
  
  table2_6[1, i] <- mean(sim6_CACE[, 1])      # true value
  table2_6[2, i] <- mean(sim6_CACE[, 6*i-4])  # mean
  table2_6[3, i] <- mean(sim6_CACE[, 6*i-2])  # median
  table2_6[4, i] <- mean(sim6_CACE[, 6*i-4]-sim6_CACE[, 1])  # bias
  table2_6[5, i] <- mean(sim6_CACE[, 6*i+1])  # CI length
  table2_6[6, i] <- sum(sim6_CACE[, 6*i])/R   # Coverage prob
  table2_6[7, i] <- mean(sim6_CACE[, 56])      # mean true
  
  table2_7[1, i] <- mean(sim7_CACE[, 1])      # true value
  table2_7[2, i] <- mean(sim7_CACE[, 6*i-4])  # mean
  table2_7[3, i] <- mean(sim7_CACE[, 6*i-2])  # median
  table2_7[4, i] <- mean(sim7_CACE[, 6*i-4]-sim7_CACE[, 1])  # bias
  table2_7[5, i] <- mean(sim7_CACE[, 6*i+1])  # CI length
  table2_7[6, i] <- sum(sim7_CACE[, 6*i])/R   # Coverage prob
  table2_7[7, i] <- mean(sim7_CACE[, 56])      # mean true
  
}

table1  <- rbind(table1_6, table1_7)
table2 <- rbind(table2_6[c(2, 4, 5, 6), ], 
                table2_7[c(2, 4, 5, 6), ])


# Save the data
setwd("../")
#write.table(table1_5, "table1_5.txt", sep="\t", row.names = FALSE)
#write.table(table2_5, "table2_5.txt", sep="\t", row.names = FALSE)
write.table(table1_6, "table1_6.txt", sep="\t", row.names = FALSE)
write.table(table2_6, "table2_6.txt", sep="\t", row.names = FALSE)
write.table(table1_7, "table1_7.txt", sep="\t", row.names = FALSE)
write.table(table2_7, "table2_7.txt", sep="\t", row.names = FALSE)

write.table(table1, "table1_67.txt", sep="\t", row.names = FALSE)
write.table(table2, "table2_67.txt", sep="\t", row.names = FALSE)
