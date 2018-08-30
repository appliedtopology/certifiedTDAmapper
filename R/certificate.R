# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'





### Try to generate a certificate for a TDAMapper analysis

#' @importFrom magrittr %>%
#' @importFrom tibble as.tibble
#' @importFrom dplyr summarise_all
#' @importFrom TDA ripsDiag
#' @importFrom igraph E V ends graph_from_adjacency_matrix
NULL

calculate.persistence = function(pts) {
  max.eps = (pts %>% dist %>% max) * 2
  dgm = ripsDiag(pts, NCOL(pts), max.eps)
  dgm$diagram[dgm$diagram == max.eps] = Inf
  diagram = dgm$diagram[-1,]
  if(NCOL(diagram) != 3)
    data.frame(max.death=diagram[3],
               max.eps=diagram[3]-diagram[2],
               max.log.eps=log(diagram[3]-diagram[2]))
  else
    data.frame(max.death=max(diagram[,3]),
               max.eps=max(diagram[,3]-diagram[,2]),
               max.log.eps=max(log(diagram[,3]-diagram[,2])))
}

random.like = function(pts) {
  do.call(cbind, lapply(1:NCOL(pts), function(j) {
    lo = min(pts[,j])
    hi = max(pts[,j])
    N = NROW(pts)
    runif(N, (lo*N-hi)/(N-1), (hi*N-lo)/(N-1))
  }))
}


#' Mapper graph certificates
#'
#' Calculates certificates of non-obstruction for \code{TDAmapper::mapper1d} graphs
#'
#' Based on Vejdemo-Johansson and Leshchenko, this is an implementation of computing
#' certificates for 1-dimensional mapper graphs. Here, a certificate is a quantification
#' of the absence or indication of the presence of an obstruction to a nerve lemma for
#' the graph, which would imply that the graph is close in persistent homology to the
#' true shape of the point cloud.
#'
#' This function implements the global method from Vejdemo-Johansson and Leshchenko for
#' probabilistic certificates.
#'
#' @param mapper A \code{TDAmapper::mapper1d} output object.
#' @param data The underlying data for the Mapper graph.
#' @param graph The adjacency matrix \code{igraph::graph} of the Mapper graph.
#'              Will be calculated if missing.
#' @param type Which kind(s) of certificate to calculate: interleaving distance or
#'             probabilistic test. Takes a character vector and calculates either
#'             or both depending on the presence or absense of \code{interleaving}
#'             or \code{probabilistic} as strings in the vector.
#' @param cl A \code{parallel} cluster object. If present, will use \code{parLapply}
#'           to speed up calculations.
#' @param N.sim The number of entries in the simulation testing for the probabilistic
#'              certificate.
#' @return A list with entries depending on the type of certificate chosen:
#'   If calculating an \code{interleaving} certificate
#'   \describe{
#'     \item{\code{interleaving.epsilon}}{The interleaving distance according to the
#'       persistent interleaving results by Govc and Skraba.}
#'     \item{\code{persistences}}{The actual persistence values used to calculate
#'       \code{interleaving.epsilon}}
#'   }
#'
#'   If calculating an \code{probabilistic} certificate
#'   \describe{
#'     \item{\code{z.scores}}{The Z-scores calculated.}
#'     \item{\code{p.value}}{The estimated p-value derived.}
#'   }
#' @export
certificate = function(mapper, data, graph=NULL, type=c("interleaving", "probabilistic"), cl=NULL,
                       N.sim=100) {
  if(is.null(graph)) {
    graph = graph_from_adjacency_matrix(mapper$adjacency, mode="undirected")
  }

  tasks = c(
    lapply(E(graph), function(e) {
      vs = ends(graph, e)
      pts = intersect(mapper$points_in_vertex[[vs[1]]], mapper$points_in_vertex[[vs[2]]])
      pts
    }),
    mapper$points_in_vertex
  )

  ret.struct = list(tasks=tasks)

  if("interleaving" %in% type) {
    persistences = do.call(rbind, lapply(tasks, function(task) {
      points = data[task,]
      if(dim(points)[1] <= 1) {
        data.frame(max.death=0, max.eps=0, max.log.eps=-Inf)
      } else {
        calculate.persistence(points)
      }
    }))
    # Govc-Skraba: eps is bounded by 2*(d+1)*max(persistences) where d is the dimension
    # of the Mapper complex. Since we are working with graphs, d+1=2.
    ret.struct$interleaving.eps = 2*2*max(persistences %>% select(max.death, max.eps))
    ret.struct$persistences = persistences
  }

  if("probabilistic" %in% type) {
    if(!is.null(cl)) {
      parallel::clusterEvalQ(cl, {
        library(tidyverse)
        library(TDA)
        library(TDAmapper)
        library(igraph)
        calculate.persistence = function(pts) {
          max.eps = (pts %>% dist %>% max) * 2
          dgm = ripsDiag(pts, NCOL(pts), max.eps)
          dgm$diagram[dgm$diagram == max.eps] = Inf
          diagram = dgm$diagram[-1,]
          if(NCOL(diagram) != 3)
            data.frame(max.death=diagram[3],
                       max.eps=diagram[3]-diagram[2],
                       max.log.eps=log(diagram[3]-diagram[2]))
          else
            data.frame(max.death=max(diagram[,3]),
                       max.eps=max(diagram[,3]-diagram[,2]),
                       max.log.eps=max(log(diagram[,3]-diagram[,2])))
        }
      })
    }
    z.scores = do.call(rbind, lapply(tasks, function(task) {
      if(length(task) == 1)
        return(numeric(N.sim))
      points = data[task,]
      simulations = lapply(1:(N.sim-1), function(j) random.like(points))
      if(is.null(cl))
        sim.persistences = do.call(rbind,lapply(cl, simulations, calculate.persistence))
      else
        sim.persistences = do.call(rbind,parallel::parLapply(cl, simulations, calculate.persistence))
      persistences = rbind(
        calculate.persistence(points),
        sim.persistences
      )
      z.scores = (persistences$max.log.eps - mean(persistences$max.log.eps[-1], na.rm=T))/sd(persistences$max.log.eps[-1],na.rm=T)
      z.scores
    }))
    p.rank = (z.scores %>% as.tibble %>% summarise_all(max) %>% t %>% rank)[1]
    ret.struct$p.value = 1-p.rank/N.sim
    ret.struct$z.scores = z.scores
  }
  return(ret.struct)
}

