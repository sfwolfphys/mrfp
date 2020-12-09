#############################################################
##
## mrfp.R
## Author: SFW
##
## This function will make a reduced network from an existing partition
## (e.g., edge_betweenness) rather than generating it from concor.
#############################################################


library(concorR)
library(igraph)

## Exists in the repo
.edge_dens <- function(adj_mat) {
    adj_mat[adj_mat > 0] <- 1

    a <- sum(adj_mat)
    m <- length(adj_mat) - sqrt(length(adj_mat))
    d <- a / m
    return(d)
}

.scaledDegree <-  function(adj_mat){
    adj_mat[adj_mat > 0] <- 1

    avgOutDegree = sum(adj_mat)/nrow(adj_mat)
    maxOutDegree = max(rowSums(adj_mat))
    scaledDegree = ifelse(maxOutDegree > 0,avgOutDegree/maxOutDegree, 0)
    return(scaledDegree)
}


###############################################################
## New functions to add to repo.

.block_edge_dens  <- function(adj_mat){
    adj_mat[adj_mat > 0]  <- 1

    a <- sum(adj_mat)
    m <- length(adj_mat)
    d <- a / m
    return(d)
}    

make_reduced_from_partition <- function(adj_mat, partition, stat='density') {
  if(stat=='density'){  
      dens <- .edge_dens(adj_mat)

      nb = max(partition)
      reduced_den = matrix(0, nrow = nb, ncol = nb)
      rownames(reduced_den) = paste("Block",1:nb)
      colnames(reduced_den) = paste("Block",1:nb)
      for(j in 1:nb){
          nRows = sum(j==partition)
          for(k in 1:nb){
              nCols = sum(k==partition)
              if(nRows==1){
                  if(nCols==1){
                      blk_adj_mat = adj_mat[j==partition, k==partition] 
                      d = ifelse(blk_adj_mat>0,1,0) 
                  }else{
                      blk_adj_mat = adj_mat[j==partition, k==partition] 
                      blk_adj_mat = matrix(blk_adj_mat,nrow=1)
                      d = .block_edge_dens(blk_adj_mat)
                  }
              }else{
                  if(nCols==1){
                      blk_adj_mat = adj_mat[j==partition, k==partition]
                      blk_adj_mat = matrix(blk_adj_mat,ncol=1)
                  }else{
                      blk_adj_mat = adj_mat[j==partition, k==partition]
                  }
                  d = ifelse(j==k,.edge_dens(blk_adj_mat),
                             .block_edge_dens(blk_adj_mat))
              }
              reduced_den[j,k] = d
          }
      }
      reduced_den[is.nan(reduced_den)] <- 0
      reduced_den[reduced_den < dens] <- 0
      reduced_den[reduced_den > 0] <- 1
 
      return_list <- list()
      return_list$reduced_mat <- reduced_den
      return_list$dens <- dens
      return(return_list)
  }else if(stat=='degree'){
      outdegree = .scaledDegree(adj_mat)

      nb = max(partition)
      reduced_degree = matrix(0, nrow = nb, ncol = nb)
      rownames(reduced_degree) = paste("Block",1:nb)
      colnames(reduced_degree) = paste("Block",1:nb)
      for(j in 1:nb){
          nRows = sum(j==partition)
          for(k in 1:nb){
              nCols = sum(k==partition)
              if(nRows==1){
                  if(nCols==1){
                      blk_adj_mat = adj_mat[j==partition, k==partition] 
                      outDeg = ifelse(blk_adj_mat>0,1,0) 
                  }else{
                      blk_adj_mat = adj_mat[j==partition, k==partition] 
                      blk_adj_mat = matrix(blk_adj_mat,nrow=1)
                      outDeg = .scaledDegree(blk_adj_mat)
                  }
              }else{
                  if(nCols==1){
                      blk_adj_mat = adj_mat[j==partition, k==partition]
                      blk_adj_mat = matrix(blk_adj_mat,ncol=1)
                  }else{
                      blk_adj_mat = adj_mat[j==partition, k==partition]
                  }
                  outDeg = .scaledDegree(blk_adj_mat)
              }
              reduced_degree[j,k] = outDeg
          }
      }
      reduced_degree[is.nan(reduced_degree)] <- 0
      reduced_degree[reduced_degree < outdegree] <- 0
      reduced_degree[reduced_degree > 0] <- 1
      
      return_list <- list()
      return_list$reduced_mat <- reduced_degree
      return_list$deg <- outdegree
      return(return_list)
  }else{
      stop('Statistics implemented for determining edges in reduced networks are only 
         density and degree.')
  }
}

