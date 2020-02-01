args = commandArgs(trailingOnly = TRUE)
paintor <- read.table(args[1], header = T)
causal <- read.table(args[2])
snps <- dim(paintor)[1]
paintor$Posterior_Prob <- paintor$Posterior_Prob/sum(paintor$Posterior_Prob)
x = max(paintor$Posterior_Prob)
summ =  x
set <- paintor[which(paintor$Posterior_Prob == x)[1],"pos"]
paintor <- paintor[-which(paintor$Posterior_Prob == x)[1],]
count = 1
while( summ < 0.95){
  x = max(paintor$Posterior_Prob)
  summ = summ + x
  set <-c(set, paintor[which(paintor$Posterior_Prob == x)[1],"pos"])
  paintor <- paintor[-which(paintor$Posterior_Prob == x)[1],]
  count = count + 1
}
set <- as.data.frame(set)
write.table(set, args[3], append = T, col.names = F, row.names = F, quote = F)
write.table(length(which(causal$V1 %in% set$set))/length(causal$V1), args[4], append = T, col.names = F, row.names = F, quote = F)
write.table(count, args[5], append = T, col.names = F, row.names = F, quote = F)
