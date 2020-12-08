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

## set.seed(1234)
## g  <-  erdos.renyi.game(50,p=0.2)
## g_adj  <- as.matrix(as_adjacency_matrix(g))
## ebc.g  <- edge.betweenness.community(g)
## ebPart  <- list(ebc.g$membership)
## g.red  <- make_reduced_from_partition(list(g_adj), ebPart, stat='degree')
## plot_reduced(make_reduced_igraph(g.red$reduced_mat[[1]]))

## g.red.den  <- make_reduced_from_partition(list(g_adj), ebPart, stat='density')
## plot_reduced(make_reduced_igraph(g.red.den$reduced_mat[[1]]))


make_reduced_from_partition <- function(adj_list, partition_list, stat='density') {
  if(stat=='density'){  
      dens_vec <- sapply(adj_list, function(x) .edge_dens(x))
      mat_return <- vector("list", length = length(dens_vec))

      for(i in 1:length(dens_vec)){
          this_adj_mat = adj_list[[i]]
          thisBlk = partition_list[[i]]
          nb = max(thisBlk)
          reduced_den = matrix(0, nrow = nb, ncol = nb)
          rownames(reduced_den) = paste("Block",1:nb)
          colnames(reduced_den) = paste("Block",1:nb)
          for(j in 1:nb){
              nRows = sum(j==thisBlk)
              for(k in 1:nb){
                  nCols = sum(k==thisBlk)
                  if(nRows==1){
                      if(nCols==1){
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk] 
                          d = ifelse(blk_adj_mat>0,1,0) 
                      }else{
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk] 
                          blk_adj_mat = matrix(blk_adj_mat,nrow=1)
                          d = .block_edge_dens(blk_adj_mat)
                      }
                  }else{
                      if(nCols==1){
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk]
                          blk_adj_mat = matrix(blk_adj_mat,ncol=1)
                      }else{
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk]
                      }
                      d = ifelse(i==j,.edge_dens(blk_adj_mat),
                                      .block_edge_dens(blk_adj_mat))
                  }
                  reduced_den[j,k] = d
              }
          }
          temp1 <- reduced_den
          temp1[is.nan(temp1)] <- 0
          temp1[temp1 < dens_vec[[i]]] <- 0
          temp1[temp1 > 0] <- 1
          mat_return[[i]] <- temp1
      }
 
      return_list <- list()
      return_list$reduced_mat <- mat_return
      return_list$dens <- dens_vec
      return(return_list)
  }else if(stat=='degree'){
      outdegree = lapply(adj_list, function(x) .scaledDegree(x))
      mat_return <- vector("list", length = length(outdegree))
      
      for(i in 1:length(outdegree)){ 
          this_adj_mat = adj_list[[i]]
          thisBlk = partition_list[[i]]
          nb = max(thisBlk)
          reduced_degree = matrix(0, nrow = nb, ncol = nb)
          rownames(reduced_degree) = paste("Block",1:nb)
          colnames(reduced_degree) = paste("Block",1:nb)
          for(j in 1:nb){
              nRows = sum(j==thisBlk)
              for(k in 1:nb){
                  nCols = sum(k==thisBlk)
                  if(nRows==1){
                      if(nCols==1){
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk] 
                          outDeg = ifelse(blk_adj_mat>0,1,0) 
                      }else{
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk] 
                          blk_adj_mat = matrix(blk_adj_mat,nrow=1)
                          outDeg = .scaledDegree(blk_adj_mat)
                      }
                  }else{
                      if(nCols==1){
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk]
                          blk_adj_mat = matrix(blk_adj_mat,ncol=1)
                      }else{
                          blk_adj_mat = this_adj_mat[j==thisBlk, k==thisBlk]
                      }
                      outDeg = .scaledDegree(blk_adj_mat)
                  }
                  reduced_degree[j,k] = outDeg
              }
          }
          temp1 <- reduced_degree
          temp1[is.nan(temp1)] <- 0
          temp1[temp1 < outdegree[[i]]] <- 0
          temp1[temp1 > 0] <- 1
          mat_return[[i]] <- temp1
      }
                                                                
      return_list <- list()
      return_list$reduced_mat <- mat_return
      return_list$deg <- outdegree
      return(return_list)
  }else{
      stop('Statistics implemented for determining edges in reduced networks are only 
         density and degree.')
  }
}

.scaledDegree <-  function(adj_mat){
    adj_mat[adj_mat > 0] <- 1

    avgOutDegree = sum(adj_mat)/nrow(adj_mat)
    maxOutDegree = max(rowSums(adj_mat))
    scaledDegree = ifelse(maxOutDegree > 0,avgOutDegree/maxOutDegree, 0)
    return(scaledDegree)
}

.block_edge_dens  <- function(adj_mat){
    adj_mat[adj_mat > 0]  <- 1

    a <- sum(adj_mat)
    m <- length(adj_mat)
    d <- a / m
    return(d)
}    

.edge_dens <- function(adj_mat) {
    adj_mat[adj_mat > 0] <- 1

    a <- sum(adj_mat)
    m <- length(adj_mat) - sqrt(length(adj_mat))
    d <- a / m
    return(d)
}
