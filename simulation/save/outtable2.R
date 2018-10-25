R <- 100

DIC_1_set_1 <- read.table('DIC_1_set_1.txt')
DIC_1_set_2 <- read.table('DIC_1_set_2.txt')
DIC_1_set_3 <- read.table('DIC_1_set_3.txt')
DIC_1_set_4 <- read.table('DIC_1_set_4.txt')

DIC_1 <- rbind(DIC_1_set_1, DIC_1_set_2, DIC_1_set_3, DIC_1_set_4)

sim1_loc <- rep(NA, R)
for (i in 1:R) {
  sim1_loc[i] <- which.min(DIC_1[i,])
}
outtable1_1 <- matrix(NA, 1, 8)
for (i in 1:8) {
  outtable1_1[i] <- length(which(sim1_loc == i)) 
}
outtable1_1


DIC_2_set_1 <- read.table('DIC_2_set_1.txt')
DIC_2_set_2 <- read.table('DIC_2_set_2.txt')
DIC_2_set_3 <- read.table('DIC_2_set_3.txt')
DIC_2_set_4 <- read.table('DIC_2_set_4.txt')

DIC_2 <- rbind(DIC_2_set_1, DIC_2_set_2, DIC_2_set_3, DIC_2_set_4)

sim2_loc <- rep(NA, R)
for (i in 1:R) {
  sim2_loc[i] <- which.min(DIC_2[i,])
}
outtable1_2 <- matrix(NA, 1, 8)
for (i in 1:8) {
  outtable1_2[i] <- length(which(sim2_loc == i)) 
}
outtable1_2


DIC_3_set_1 <- read.table('DIC_3_set_1.txt')
DIC_3_set_2 <- read.table('DIC_3_set_2.txt')
DIC_3_set_3 <- read.table('DIC_3_set_3.txt')
DIC_3_set_4 <- read.table('DIC_3_set_4.txt')

DIC_3 <- rbind(DIC_3_set_1, DIC_3_set_2, DIC_3_set_3, DIC_3_set_4)

sim3_loc <- rep(NA, R)
for (i in 1:R) {
  sim3_loc[i] <- which.min(DIC_3[i,])
}
outtable1_3 <- matrix(NA, 1, 8)
for (i in 1:8) {
  outtable1_3[i] <- length(which(sim3_loc == i)) 
}
outtable1_3


DIC_4_set_1 <- read.table('DIC_4_set_1.txt')
DIC_4_set_2 <- read.table('DIC_4_set_2.txt')
DIC_4_set_3 <- read.table('DIC_4_set_3.txt')
DIC_4_set_4 <- read.table('DIC_4_set_4.txt')

DIC_4 <- rbind(DIC_4_set_1, DIC_4_set_2, DIC_4_set_3, DIC_4_set_4)

sim4_loc <- rep(NA, R)
for (i in 1:R) {
  sim4_loc[i] <- which.min(DIC_4[i,])
}
outtable1_4 <- matrix(NA, 1, 8)
for (i in 1:8) {
  outtable1_4[i] <- length(which(sim4_loc == i)) 
}
outtable1_4

outtable1 <- rbind(outtable1_1, outtable1_2, outtable1_3, outtable1_4)
write.table(outtable1, "outtable1.txt", sep="\t")

outtable2 <- rbind(outtable2_1, outtable2_2, outtable2_3, outtable2_4)
write.table(outtable2, "outtable2.txt", sep="\t")

