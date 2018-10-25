R <- 2000
sim8_CACE <- matrix(NA, R, 56)
sim8_DIC <- matrix(NA, R, 9)
sim8_pD <- matrix(NA, R, 9)
sim8_loc <- rep(NA, R)
for (i in 1:R){
  setwd("../sim_CACE")
  sim8_CACE[i, ] <- read.table(paste("sim8_CACE_",i,".txt",sep=""))$x

  setwd("../sim_DIC")
  sim8_DIC[i, ] <- read.table(paste("sim8_DIC_",i,".txt",sep=""))$x
  sim8_loc[i] <- which.min(sim8_DIC[i,])

  setwd("../sim_pD")
  sim8_pD[i, ] <- read.table(paste("sim8_pD_",i,".txt",sep=""))$x
}

table1_8  <- matrix(NA, 1, 9)
for (i in 1:9) {
  table1_8[i] <- length(which(sim8_loc == i)) 
}

table2_8 <- matrix(NA, 7, 9)
for (i in 1:9) {
  
  table2_8[1, i] <- mean(sim8_CACE[, 1])      # true value
  table2_8[2, i] <- mean(sim8_CACE[, 6*i-4])  # mean
  table2_8[3, i] <- mean(sim8_CACE[, 6*i-2])  # median
  table2_8[4, i] <- mean(sim8_CACE[, 6*i-4]-sim8_CACE[, 1])  # bias
  table2_8[5, i] <- mean(sim8_CACE[, 6*i+1])  # CI length
  table2_8[6, i] <- sum(sim8_CACE[, 6*i])/R   # Coverage prob
  table2_8[7, i] <- mean(sim8_CACE[, 56])      # mean true
}

#table1  <- rbind(table1_6, table1_7)
#table2 <- rbind(table2_6[c(2, 4, 5, 6), ], 
#                table2_7[c(2, 4, 5, 6), ])


# Save the data
setwd("../")
write.table(table1_8, "table1_8.txt", sep="\t", row.names = FALSE)
write.table(table2_8, "table2_8.txt", sep="\t", row.names = FALSE)

#write.table(table1, "table1_67.txt", sep="\t", row.names = FALSE)
#write.table(table2, "table2_67.txt", sep="\t", row.names = FALSE)
