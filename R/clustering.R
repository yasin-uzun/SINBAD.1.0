library(densityClust)

dens_clus <- function(dim_red, rho_threshold=1, delta_threshold=1){
  dataDist <- dist(dim_red)
  dataClust <- densityClust::densityClust(dataDist, gaussian = T)
  dataClust <- densityClust::findClusters(dataClust, 
                                          rho = rho_threshold,
                                          delta = delta_threshold)
  clusters = factor(dataClust$clusters)
  head(clusters)
  return(clusters)
}
