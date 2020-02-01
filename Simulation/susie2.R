library(susieR)
set.seed(1)
args = commandArgs(trailingOnly = TRUE)
# [1] is summary statistics,[2] is LD matrix, [3] is the true causal file, [4] is number of causal snps,
# [5] is output subsets file name, [6] is output set file name, [7] is accuracy (sensitivity) file name, 
# [8] is set size file name, [9] is the number of credible sets file name.
# z reports line 1 did not have 10 elements
z <- read.table(args[1], header = FALSE)
R <- read.table(args[2], header = FALSE)
causal <- read.table(args[3], header = FALSE)
num_causal <- as.numeric(args[4])

R <- data.matrix(R)

fitted <- susie_rss(z[,2], R,
                L = 10,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                verbose = TRUE, check_R = FALSE)

# number of causal sets (CS)
num_cs <- length(fitted$sets$cs)

if (num_cs == 0) {
  # when susie does not converge
  write.table(0, args[6], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(0, args[7], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(0, args[8], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
} else {
  for (each in 1:num_cs) {
    subset <- fitted$sets$cs[[each]] # the n-th causal set among all causal sets
    subset <- c(subset)
    union <- which(causal$V1 %in% subset)
    # subname <- paste(args[5], each, sep = "")
    write.table(paste("cs", each, sep=""), args[5], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(subset, args[5], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    # sensitivity
  }
  set <- do.call(c, fitted$sets$cs)
  set <- as.data.frame(set)

  # adjust the set index
  for (loci in 1:nrow(set)) {
    set[loci,] = set[loci,] - 1
  }
  
  count <- length(set[,1])
  
  write.table(set, args[6], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(length(which(causal$V1 %in% set$set))/length(causal$V1), args[7], append = TRUE, 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(count, args[8], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
}

write.table(num_cs, args[9], append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)