R <- 100

DIC_1_set_1 <- read.table('DIC_1_set_1.txt')
DIC_1_set_2 <- read.table('DIC_1_set_2.txt')
DIC_1_set_3 <- read.table('DIC_1_set_3.txt')
DIC_1_set_4 <- read.table('DIC_1_set_4.txt')

DIC_1 <- rbind(DIC_1_set_1, DIC_1_set_2, DIC_1_set_3, DIC_1_set_4)

R <- 100
sim1_loc <- rep(NA, R)
for (i in 1:R) {
  sim1_loc[i] <- which.min(DIC_1[i,])
}
outtable1_1 <- matrix(NA, 1, 8)
for (i in 1:8) {
  outtable1_1[i] <- length(which(sim1_loc == i)) 
}
outtable1_1


CACE_1_set_1 <- read.table('CACE_1_set_1.txt')
CACE_1_set_2 <- read.table('CACE_1_set_2.txt')
CACE_1_set_3 <- read.table('CACE_1_set_3.txt')
CACE_1_set_4 <- read.table('CACE_1_set_4.txt')

CACE_1 <- rbind(CACE_1_set_1, CACE_1_set_2, CACE_1_set_3, CACE_1_set_4)

outtable2_1 <- matrix(NA, 4, 8)
for (i in 1:8) {
  outtable2_1[1, i] <- mean(CACE_1[, 6*i-4])  # mean
  outtable2_1[2, i] <- mean(CACE_1[, 6*i-2])  # median
  outtable2_1[3, i] <- mean(CACE_1[, 6*i+1])  # CI length
  outtable2_1[4, i] <- sum(CACE_1[, 6*i])/R   # Coverage prob
}
outtable2_1


CACE_2_set_1 <- read.table('CACE_2_set_1.txt')
CACE_2_set_2 <- read.table('CACE_2_set_2.txt')
CACE_2_set_3 <- read.table('CACE_2_set_3.txt')
CACE_2_set_4 <- read.table('CACE_2_set_4.txt')

CACE_2 <- rbind(CACE_2_set_1, CACE_2_set_2, CACE_2_set_3, CACE_2_set_4)

outtable2_2 <- matrix(NA, 4, 8)
for (i in 1:8) {
  outtable2_2[1, i] <- mean(CACE_2[, 6*i-4])  # mean
  outtable2_2[2, i] <- mean(CACE_2[, 6*i-2])  # median
  outtable2_2[3, i] <- mean(CACE_2[, 6*i+1])  # CI length
  outtable2_2[4, i] <- sum(CACE_2[, 6*i])/R   # Coverage prob
}
outtable2_2


CACE_1_set_1 <- read.table('CACE_1_set_1.txt')
CACE_1_set_2 <- read.table('CACE_1_set_2.txt')
CACE_1_set_3 <- read.table('CACE_1_set_3.txt')
CACE_1_set_4 <- read.table('CACE_1_set_4.txt')

CACE_1 <- rbind(CACE_1_set_1, CACE_1_set_2, CACE_1_set_3, CACE_1_set_4)

outtable2_1 <- matrix(NA, 4, 8)
for (i in 1:8) {
  outtable2_1[1, i] <- mean(CACE_1[, 6*i-4])  # mean
  outtable2_1[2, i] <- mean(CACE_1[, 6*i-2])  # median
  outtable2_1[3, i] <- mean(CACE_1[, 6*i+1])  # CI length
  outtable2_1[4, i] <- sum(CACE_1[, 6*i])/R   # Coverage prob
}
outtable2_1


CACE_3_set_1 <- read.table('CACE_3_set_1.txt')
CACE_3_set_2 <- read.table('CACE_3_set_2.txt')
CACE_3_set_3 <- read.table('CACE_3_set_3.txt')
CACE_3_set_4 <- read.table('CACE_3_set_4.txt')

CACE_3 <- rbind(CACE_3_set_1, CACE_3_set_2, CACE_3_set_3, CACE_3_set_4)

outtable2_3 <- matrix(NA, 4, 8)
for (i in 1:8) {
  outtable2_3[1, i] <- mean(CACE_3[, 6*i-4])  # mean
  outtable2_3[2, i] <- mean(CACE_3[, 6*i-2])  # median
  outtable2_3[3, i] <- mean(CACE_3[, 6*i+1])  # CI length
  outtable2_3[4, i] <- sum(CACE_3[, 6*i])/R   # Coverage prob
}
outtable2_3


CACE_4_set_1 <- read.table('CACE_4_set_1.txt')
CACE_4_set_2 <- read.table('CACE_4_set_2.txt')
CACE_4_set_3 <- read.table('CACE_4_set_3.txt')
CACE_4_set_4 <- read.table('CACE_4_set_4.txt')

CACE_4 <- rbind(CACE_4_set_1, CACE_4_set_2, CACE_4_set_3, CACE_4_set_4)

outtable2_4 <- matrix(NA, 4, 8)
for (i in 1:8) {
  outtable2_4[1, i] <- mean(CACE_4[, 6*i-4])  # mean
  outtable2_4[2, i] <- mean(CACE_4[, 6*i-2])  # median
  outtable2_4[3, i] <- mean(CACE_4[, 6*i+1])  # CI length
  outtable2_4[4, i] <- sum(CACE_4[, 6*i])/R   # Coverage prob
}
outtable2_4