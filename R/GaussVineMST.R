#' Sequential MST based on partial correlations
#'
#' @description
#' Sequential minimum spanning tree (MST) based on partial correlations for best truncated Gaussian vine.
#' This depends on the R package igraph for the minimum spanning tree algorithm
#' and operations on graphs.
#'
#' @param rmat correlation matrix
#' @param n sample size
#' @param CFIbd lower bound of CFI (comparative fit index) to stop, default 0.95,
#' set as 1.01 for no truncation
#' @param iprint print flag for intermediate results on truncated vine
#' @param iprint2 print flag for intermediate results on graph objects
#'
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{VineA} -  vine array of dimension dxd
#'  \item \code{pcmat} -  matrix of partial correlations by tree
#'  \item \code{treeweight} -  vector of length d-1 with
#'                 sum_edge log(1-rho[edge]^2) for trees 1,...d-1
#'  \item \code{edgemat} -  matrix with columns node1 node2 level vector-of-conditioning
#'  \item \code{fitval} - cumsum(treeweight)/\eqn{sum_{1:(d-1)}} treeweight
#'  \item \code{CFIv} - vector of CFI values
#'  \item \code{ntrunc} - truncation level to reach CFIbd
#' }
#'
#' @details
#' CFI = 1- numerator/denominator, \cr
#' numerator = \eqn{\max(0,D_t-\nu_t)},\cr
#' denominator = \eqn{\max(0,D_t-\nu_t,D_0-\nu_0)},\cr
#' \eqn{D_0=-n\log(\det(R))},\cr
#' \eqn{D_t = n[-L_t(V) - \log(\det(R))]} is decreasing as t increases,\cr
#' \eqn{L_t(V) = -\sum_{1:t} \sum_e -\log(1-r_e^2)} (sum pcor up to tree t),\cr
#' \eqn{\nu_t = (d-t)(d-t-1)/2; \nu_0=d(d-1)/2},\cr
#' \eqn{L_t(V)} is increasing in t.
#'
#' @references
#' For use of CFI in the context of partial correlation vines, please see: \cr
#' Brechmann E C and Joe H (2015), Truncation of vine copulas using fit indices.
#' Journal of Multivariate Analysis, 138, 19-33.
#'
#' @examples
#' n=400
#' rmat=matrix(c(
#'   1.000,  0.172, -0.062,  0.385,  0.499, 0.267,  0.578,
#'   0.172,  1.000, -0.047,  0.493,  0.602, 0.396,  0.600,
#'  -0.062, -0.062,  1.000, -0.120, -0.092, 0.070, -0.070,
#'   0.385,  0.493, -0.120,  1.000,  0.607, 0.557,  0.742,
#'   0.499,  0.602, -0.092,  0.607,  1.000, 0.672,  0.883,
#'   0.267,  0.396,  0.070,  0.557,  0.672, 1.000,  0.685,
#'   0.578,  0.600, -0.070,  0.742,  0.883, 0.685,  1.000),7,7)
#' vinestr095 = GaussVineMST(rmat, n, CFIbd=0.95, iprint=FALSE)
#' print(vinestr095$VineA)
#'
#' @import igraph
#' @export
GaussVineMST = function(rmat, n, CFIbd=1.01, iprint=FALSE, iprint2=FALSE)
{
  rmat = as.matrix(rmat)
  d = dim(rmat)[1]
  colnames(rmat) = rownames(rmat) = paste("V",1:d,sep="")
  treeweight = rep(NA,d-1)

  # set up objects for return
  dd=(d*(d-1))/2
  node1=rep(0,dd); node2=rep(0,dd); lev=rep(0,dd);
  condmat=matrix(0,dd,d-2)
  icomv=rep(0,d); idv=rep(0,d);
  vord=rep(0,d-1)
  treeweight=rep(0,d-1)
  pcv=rep(NA,dd);

  CFIv=rep(0,d-1)

  # simple change for singular case
  detR = det(rmat)
  if(detR<=0)
  {
    warning("Singular correlation matrix; CFI and D0, Dt not useful")
    logdet = -1.e10
  }
  else { logdet = log(detR) }

  nu0=dd
  D0=-n*logdet
  sumtreewt=0

  # tree 1
  # initialize the first graph
  # mode=upper to match C code
  g = graph.adjacency(rmat, mode="upper",weighted=TRUE,diag=FALSE)
  E(g)$name = paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep = ",")
  for (ii in 1:ecount(g))
  {
    nodes=get.edges(g,ii)
    node1[ii]=nodes[1]; node2[ii]=nodes[2]
  }
  if(iprint2) message("Initial graph:\n"); if(iprint2) print(g)

  # R indexes from 1 ; C indexes from 0
  vord[1]=1
  for(i in 2:(d-1)) { vord[i]=vord[i-1]+(d-i+1); }

  mst = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))
  treeweight[1] = sum(log(1-E(mst)$weight^2))
  sumtreewt = sumtreewt+treeweight[1]
  Dt=n*(sumtreewt-logdet); nut=((d-1)*(d-2))/2
  numer=max(0,Dt-nut); denom=max(numer+1e-5,D0-nu0)
  itree=1
  CFIv[itree]=1-numer/denom
  if(iprint2) message("\n*** tree1 ***"); if(iprint2) print(mst)

  # add info to objects: node1 node2 lev pcv
  for(i in 1:(d-1))
  {
    nodes=get.edges(mst,i)
    v1=nodes[1]; v2=nodes[2]
    eid=vord[v1]+(v2-v1-1)
    node1[eid]=v1; node2[eid]=v2;
    lev[eid]=1; pcv[eid]=rmat[v1,v2]
  }
  if(iprint2) { tmat=cbind(node1,node2,lev,pcv,condmat); message("Tree1 edges:"); print(tmat) }

  # loop through remaining trees
  if(CFIv[1]<CFIbd)
  {
    for(itree in 2:(d-1))
    {
      nres = ecount(mst)
      if(iprint2) {
        message("\n*** setting up tree ", itree, " ***")
        message("*** ecount=", nres, "***")
      }
      g = graph.full(nres) # all possible edges
      V(g)$name = E(mst)$name
      if(itree>2) V(g)$iname = E(mst)$iname
      ii=0
      for (i in 1:ecount(g))
      {
        # g is full graph so it goes thru all pairs and check proximity
        if(itree==2)
        {
          id = ends(g,i,names=F) # internal indices for edge
          temp = get.edges(mst,id) # edges of prev mst tree
          avec=temp[1,]; bvec=temp[2,]  # should be increasing
        }
        else
        {
          jj = ends(g,i,names=F) # internal indices for edge, 2-vector
          # different from itree=2 in next 4 lines
          id1=idv[jj[1]]; id2=idv[jj[2]]
          avec=c(node1[id1],node2[id1],condmat[id1,1:(itree-2)])
          bvec=c(node1[id2],node2[id2],condmat[id2,1:(itree-2)])
          avec=sort(avec); bvec=sort(bvec)
        }
        chk=proxCheck(avec,bvec)
        if(iprint2) message(chk$iok,", ", paste(avec,collapse=","),":",paste(bvec,collapse=","))
        ok = chk$iok
        if (ok)
        {
          icomv=chk$icomv[1:(itree-1)]
          icond1=chk$icond1; icond2=chk$icond2;
          if(itree==2)
          {
            E(g)[i]$name = paste(paste(c(icond1,icond2),collapse = ","),icomv,sep="|")
            tem=(rmat[icond1,icond2]-rmat[icomv,icond1]*rmat[icomv,icond2])/
              sqrt((1.-rmat[icomv,icond1]*rmat[icomv,icond1])*
                     (1.-rmat[icomv,icond2]*rmat[icomv,icond2]))
          }
          else
          {
            E(g)[i]$name = paste(paste(c(icond1,icond2),collapse=","),
                                 paste(icomv,collapse = ","),
                                 sep="|")
            tem = partCor(rmat,icomv,icond1,icond2)
          }
          E(g)[i]$weight = tem
          ii=ii+1
          eid=vord[icond1]+(icond2-icond1-1)
          pcv[eid]=tem;
          lev[eid]=-itree; # later change to positive if in mst
          condmat[eid,1:(itree-1)]=icomv;
          if(iprint2) message(ii," ", tem," ",eid)
          E(g)[i]$iname=eid
        }
        E(g)[i]$todel = !ok
      }
      g = delete.edges(g, E(g)[E(g)$todel])
      if(iprint2) { message("Graph after proximity check:"); print(g); }
      mst = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))
      treeweight[itree] = sum(log(1-E(mst)$weight^2))
      sumtreewt = sumtreewt+treeweight[itree]
      Dt=n*(sumtreewt-logdet); nut=((d-itree)*(d-itree-1))/2
      numer=max(0,Dt-nut); denom=max(numer,D0-nu0)
      CFIv[itree]=1-numer/denom
      if(iprint) message("CFI ", itree,": ", Dt," ",nut," ",numer," ",denom," ", 1-numer/denom)
      # set selected edges to have lev[eid]=itree
      #for(i in 1:(d-itree))
      # fix: do not change logic, just ensure no out-of-bounds:
      for(i in seq_len( min(d-itree, ecount(mst)) ))
      {
        eid=E(mst)[i]$iname
        v1=node1[eid]; v2=node2[eid]
        if(iprint2) message(i," ", eid," ",v1," ",v2)
        lev[eid]=itree;
        idv[i]=eid
      }
      if(iprint2) { tmat=cbind(node1,node2,lev,pcv,condmat); message("Tmat after tree ",itree,":"); print(tmat) }
      if(CFIv[itree]>=CFIbd) break;
    }
  }

  # get truncated vine array and corresponding matrix of partial correlations
  edgemat=cbind(node1,node2,lev,condmat)
  if(itree>=(d-1)) out=edges2array(d,edgemat,pcv)
  else out=edges2arrayTrunc(d,itree,edgemat,pcv,iprint=iprint)
  fitval=cumsum(treeweight)
  fitval=fitval/logdet
  if(iprint) message(CFIv)

  list(edgemat=edgemat,pcvec=pcv,VineA=out$VineA,pcmat=out$pcmat,
       treeweight=treeweight,fitval=fitval, CFIv=CFIv,ntrunc=itree)
}
