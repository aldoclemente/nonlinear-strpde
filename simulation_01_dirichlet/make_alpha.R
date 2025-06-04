alpha = 3.0
set.seed(0)
nsim = 30
alpha_0.02 = alpha + rnorm(nsim, sd = 0.02)
alpha_0.05 = alpha + rnorm(nsim, sd = 0.05)
alpha_0.1 = alpha + rnorm(nsim, sd = 0.1)

par(mfrow=c(1,3))
hist(alpha_0.02, probability = T)
hist(alpha_0.05, probability = T)
hist(alpha_0.1, probability = T)
if(!dir.exists("data")) dir.create("data")
write.table(alpha_0.02, file = "data/alpha_0.02.txt", row.names = F, col.names = F)
write.table(alpha_0.05, file = "data/alpha_0.05.txt", row.names = F, col.names = F)
write.table(alpha_0.1, file = "data/alpha_0.1.txt", row.names = F, col.names = F)
