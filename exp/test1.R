library(igraph)
library(bnlearn)
addDiffNoiseLevelsToGraphData <- function(data, noiseType = 'gaussian',
                                          noiseLevel = seq(0,1,1/(ncol(data)-1))) {
  #noiseLevel is the standard deviation of the noise of variables
  stopifnot(noiseType == 'gaussian')
  for (ci in 1:ncol(data)) {
    if(noiseLevel[ci] != 0) {
      data[, ci] = data[, ci] + rnorm(nrow(data), 0, noiseLevel[ci])
    }
  }
  return(data)
}
ham_dist_btw_matrices <- function(matrix_one, matrix_two) {
  if(nrow(matrix_one) != nrow(matrix_two)) {
    stop("The number of rows should be the same.")
  }
  if(ncol(matrix_one) != ncol(matrix_two)) {
    stop("The number of columns should be the same.")
  }
  sum = 0
  for(i in 1:nrow(matrix_one)) {
    for(j in 1:ncol(matrix_one)) {
      sum = sum + abs(matrix_one[i,j] - matrix_two[i,j])
    }
  }
  return(sum)
}

# Add different noise levels to nodes in the network
set.seed(1)
data <- readRDS("C:\\Users/Ehsan/Desktop/mmpc/data2.RData")
g <- readRDS("C:\\Users/Ehsan/Desktop/mmpc/g2.RData")
mat_true <- as.matrix(as_adj(g)) + t(as.matrix(as_adj(g))) # true network adjacency matrix
num_iterations = 10
t16_noise_allNodes <- data.frame("beforeCancel" = rep(-1,num_iterations),
                                 "afterCancel" = rep(-1,num_iterations))
for (iteration in 1:num_iterations) {
  noiseLevels <- sample(seq(0,1,0.1),ncol(data),replace = T)
  gData_noisy <- addDiffNoiseLevelsToGraphData(data,noiseLevel = noiseLevels)
  noiseLevels <- data.frame(t(noiseLevels))
  colnames(noiseLevels) <- colnames(gData_noisy)
  learned_noisy <- mmhc(data.frame(gData_noisy))
  t16_noise_allNodes[iteration,"beforeCancel"] <- ham_dist_btw_matrices(matrix_one = mat_true,
                                                                        matrix_two = amat(learned_noisy))/2
  # Cancel the noise
  learned_noise_cancel <- mmhc(data.frame(gData_noisy), noise.levels = noiseLevels)
  adj_mat_noise_cancel <- amat(learned_noise_cancel)
  t16_noise_allNodes[iteration,"afterCancel"] <- ham_dist_btw_matrices(matrix_one = mat_true,
                                                                       matrix_two = adj_mat_noise_cancel)/2
}

print(t16_noise_allNodes)
