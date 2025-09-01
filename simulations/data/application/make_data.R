
table_0 = read.table("table_0.txt", header = T)
table_1 = read.table("table_1.txt", header = T)
table_2 = read.table("table_2.txt", header = T)

obs_0 = as.matrix(table_0$Mean / table_0$Volume_mm3)
obs_1 = as.matrix(table_1$Mean / table_1$Volume_mm3)
obs_2 = as.matrix(table_2$Mean / table_2$Volume_mm3)

write.csv(format(obs_0, digits=16), file = "obs0.csv")
write.csv(format(c(obs_1, obs_2), digits=16), file = "obs.csv")
