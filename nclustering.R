##### N-Clustering sample
source('pso.R')

##### global parameters
# iterations <- 1;
# individuals <- 20;
# clusters <- 1;
# data <- array(runif(1000)*2 -1 , c(500,2));


##### general parameters
clusters_dimension <- 2;
dimensions <- clusters_dimension*clusters;
c1 <- 2;
c2 <- 1;
inertia <- .8;

# fitness function in clustering process
clustering_fitness <- function(ind) {
  val <- array(0, c(length(ind[,1]),1));
  for(individual in 1:length(ind[,1])) {
    dist <- 0;
    for(i in 1:length(data[,1])) {
      dist_vector <- (data[i,]-ind[individual,])^2;
      dist_min <- array(0, clusters);
      for (j in 1:clusters) {
        dist_min[j] <-  sqrt(sum(dist_vector[((j-1)*clusters_dimension + 1):(j*clusters_dimension)]));
      }
      dist <- dist + min(dist_min)/length(data[,1]);
    }
    val[individual] <- 1 / (1 + dist);
  }
  return(val)
}

res <- pso(FALSE, individuals, dimensions, iterations, c1, c2, inertia, clustering_fitness);

# setup margins and title size
par(mar=c(0.2,0.2,2,0.2),bg=rgb(0.9,0.9,0.9), cex.main=1.5);

# plot data
plot(data[,1], data[,2], cex= 0.1, xlab='', ylab='',lwd =1, axes=FALSE);
title(paste('#Clusters:', clusters, '#Population: ', individuals, '#Iterations: ', iterations));

distribution <- array(0,c(1,clusters), dimnames = list('Number of elements:',paste('Cluster',1:clusters)));

for(i in 1:length(data[,1])) {
  dist_vector <- (data[i,] - res$best[1:dimensions])^2;
  dist_min <- array(0, clusters);
  for (j in 1:clusters) {
    dist_min[j] <-  sqrt(sum(dist_vector[((j-1)*clusters_dimension + 1):(j*clusters_dimension)]));
  }
  dist_min <- 1 / (1 + dist_min);
  winner <-  max.col(t(dist_min));
  distribution[winner] <- distribution[winner] + 1;
  points(data[i,1], data[i,2],lwd=2,cex=(4*max(dist_min)^2), pch=19, col = winner+2);
}

for (i in 1:clusters)
  points(res$best[((i-1)*clusters_dimension)+1], res$best[(i)*clusters_dimension], col='red', pch = 3,lwd=5, cex=2);
