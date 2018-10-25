R <- 2000
sim1_CACE <-sim2_CACE <- sim3_CACE <- sim4_CACE <- matrix(NA, R, 50)
sim1_DIC <- sim2_DIC <- sim3_DIC <- sim4_DIC <- matrix(NA, R, 8)
sim1_pD <- sim2_pD <- sim3_pD <- sim4_pD <- matrix(NA, R, 8)
sim1_loc <- sim2_loc <- sim3_loc <- sim4_loc <- rep(NA, R)
for (i in 1:R){
  setwd("../sim_CACE")
  sim1_CACE[i, ] <- read.table(paste("sim1_CACE_",i,".txt",sep=""))$x
  sim2_CACE[i, ] <- read.table(paste("sim2_CACE_",i,".txt",sep=""))$x
  sim3_CACE[i, ] <- read.table(paste("sim3_CACE_",i,".txt",sep=""))$x
  sim4_CACE[i, ] <- read.table(paste("sim4_CACE_",i,".txt",sep=""))$x
  setwd("../sim_DIC")
  sim1_DIC[i, ] <- read.table(paste("sim1_DIC_",i,".txt",sep=""))$x
  sim1_loc[i] <- which.min(sim1_DIC[i,])
  sim2_DIC[i, ] <- read.table(paste("sim2_DIC_",i,".txt",sep=""))$x
  sim2_loc[i] <- which.min(sim2_DIC[i,])
  sim3_DIC[i, ] <- read.table(paste("sim3_DIC_",i,".txt",sep=""))$x
  sim3_loc[i] <- which.min(sim3_DIC[i,])
  sim4_DIC[i, ] <- read.table(paste("sim4_DIC_",i,".txt",sep=""))$x
  sim4_loc[i] <- which.min(sim4_DIC[i,])
  setwd("../sim_pD")
  sim1_pD[i, ] <- read.table(paste("sim1_pD_",i,".txt",sep=""))$x
  sim2_pD[i, ] <- read.table(paste("sim2_pD_",i,".txt",sep=""))$x
  sim3_pD[i, ] <- read.table(paste("sim3_pD_",i,".txt",sep=""))$x
  sim4_pD[i, ] <- read.table(paste("sim4_pD_",i,".txt",sep=""))$x
}

table1_1 <- table1_2 <- table1_3 <- table1_4 <- matrix(NA, 1, 8)
for (i in 1:8) {
  table1_1[i] <- length(which(sim1_loc == i)) 
  table1_2[i] <- length(which(sim2_loc == i)) 
  table1_3[i] <- length(which(sim3_loc == i)) 
  table1_4[i] <- length(which(sim4_loc == i)) 
}

table2_1 <- table2_2 <- table2_3 <- table2_4 <- matrix(NA, 7, 8)
for (i in 1:8) {
  table2_1[1, i] <- mean(sim1_CACE[, 1])      # true value
  table2_1[2, i] <- mean(sim1_CACE[, 6*i-4])  # mean
  table2_1[3, i] <- mean(sim1_CACE[, 6*i-2])  # median
  table2_1[4, i] <- mean(sim1_CACE[, 6*i-4]-sim1_CACE[, 1])  # bias, mean-true
  table2_1[5, i] <- mean(sim1_CACE[, 6*i+1])  # CI length
  table2_1[6, i] <- sum(sim1_CACE[, 6*i])/R   # Coverage prob
  table2_1[7, i] <- mean(sim1_CACE[, 50])      # mean true
  
  table2_2[1, i] <- mean(sim2_CACE[, 1])      # true value
  table2_2[2, i] <- mean(sim2_CACE[, 6*i-4])  # mean
  table2_2[3, i] <- mean(sim2_CACE[, 6*i-2])  # median
  table2_2[4, i] <- mean(sim2_CACE[, 6*i-4]-sim2_CACE[, 1])  # bias
  table2_2[5, i] <- mean(sim2_CACE[, 6*i+1])  # CI length
  table2_2[6, i] <- sum(sim2_CACE[, 6*i])/R   # Coverage prob
  table2_2[7, i] <- mean(sim2_CACE[, 50])      # mean true
  
  table2_3[1, i] <- mean(sim3_CACE[, 1])      # true value
  table2_3[2, i] <- mean(sim3_CACE[, 6*i-4])  # mean
  table2_3[3, i] <- mean(sim3_CACE[, 6*i-2])  # median
  table2_3[4, i] <- mean(sim3_CACE[, 6*i-4]-sim3_CACE[, 1])  # bias
  table2_3[5, i] <- mean(sim3_CACE[, 6*i+1])  # CI length
  table2_3[6, i] <- sum(sim3_CACE[, 6*i])/R   # Coverage prob
  table2_3[7, i] <- mean(sim3_CACE[, 50])      # mean true
  
  table2_4[1, i] <- mean(sim4_CACE[, 1])      # true value
  table2_4[2, i] <- mean(sim4_CACE[, 6*i-4])  # mean
  table2_4[3, i] <- mean(sim4_CACE[, 6*i-2])  # median
  table2_4[4, i] <- mean(sim4_CACE[, 6*i-4]-sim4_CACE[, 1])  # bias
  table2_4[5, i] <- mean(sim4_CACE[, 6*i+1])  # CI length
  table2_4[6, i] <- sum(sim4_CACE[, 6*i])/R   # Coverage prob
  table2_4[7, i] <- mean(sim4_CACE[, 50])      # mean true
}
table2_1
table1  <- rbind(table1_1, table1_2, table1_3, table1_4)
table2 <- rbind(table2_1[c(2, 4, 5, 6), ], table2_2[c(2, 4, 5, 6), ], 
                table2_3[c(2, 4, 5, 6), ], table2_4[c(2, 4, 5, 6), ])


# Save the data
setwd("../")
write.table(table1_1, "table1_1.txt", sep="\t", row.names = FALSE)
write.table(table2_1, "table2_1.txt", sep="\t", row.names = FALSE)
write.table(table1_2, "table1_2.txt", sep="\t", row.names = FALSE)
write.table(table2_2, "table2_2.txt", sep="\t", row.names = FALSE)
write.table(table1_3, "table1_3.txt", sep="\t", row.names = FALSE)
write.table(table2_3, "table2_3.txt", sep="\t", row.names = FALSE)
write.table(table1_4, "table1_4.txt", sep="\t", row.names = FALSE)
write.table(table2_4, "table2_4.txt", sep="\t", row.names = FALSE)

write.table(table1, "table1.txt", sep="\t", row.names = FALSE)
write.table(table2, "table2.txt", sep="\t", row.names = FALSE)
