
# This code was written by Jun Dong, modified by Peter Langfelder
# datExpr: expression profiles with rows=samples and cols=genes/probesets
# power: for contruction of the weighted network
# trait: the quantitative external trait
#
#' @name networkConcepts
#' @rdname networkConcepts
#' @title Calculations of network concepts
#' @description
#' This functions calculates various network concepts (topological properties,
#' network indices) of a network calculated from expression data. See details
#' for a detailed description.
#' @inheritParams adjacency
#' @param trait optional specification of a sample trait. A vector of length
#' equal the number of samples in datExpr.
#' @details
#' This function computes various network concepts (also known as network
#' statistics, topological properties, or network indices) for a weighted
#' correlation network. The nodes of the weighted correlation network will be
#' constructed between the columns (interpreted as nodes) of the input datExpr.
#' If the option networkType="unsigned" then the adjacency between nodes i and j
#'  is defined as A[i,j]=abs(cor(datExpr[,i],datExpr[,j]))^power. In the
#'  following, we use the term gene and node interchangeably since these methods
#'  were originally developed for gene networks. The function computes the
#'  following 4 types of network concepts (introduced in Horvath and Dong 2008):
#'
#' \bold{Type I}: fundamental network concepts are defined as a function of the
#' off-diagonal elements of an adjacency matrix A and/or a node significance
#' measure GS. These network concepts can be defined for any network (not just
#' correlation networks). The adjacency matrix of an unsigned weighted
#' correlation network is given by A=abs(cor(datExpr,use="p"))^power and the
#' trait based gene significance measure is given by GS= abs(cor(datExpr,trait,
#' use="p"))^power where datExpr, trait, power are input parameters.
#'
#' \bold{Type II}: conformity-based network concepts are functions of the
#' off-diagonal elements of the conformity based adjacency matrix A.CF=CF*t(CF)
#' and/or the node significance measure. These network concepts are defined for
#' any network for which a conformity vector can be defined. Details: For any
#' adjacency matrix A, the conformity vector CF is calculated by requiring that
#' A[i,j] is approximately equal to CF[i]*CF[j]. Using the conformity one can
#' define the matrix A.CF=CF*t(CF) which is the outer product of the conformity
#' vector with itself. In general, A.CF is not an adjacency matrix since its
#' diagonal elements are different from 1. If the off-diagonal elements of A.CF
#' are similar to those of A according to the Frobenius matrix norm, then A is
#' approximately factorizable. To measure the factorizability of a network, one
#' can calculate the Factorizability, which is a number between 0 and 1 (Dong
#' and Horvath 2007). T he conformity is defined using a monotonic, iterative
#' algorithm that maximizes the factorizability measure.
#'
#' \bold{Type III}: approximate conformity based network concepts are functions
#' of all elements of the conformity based adjacency matrix A.CF (including the
#' diagonal) and/or the node significance measure GS. These network concepts are
#'  very useful for deriving relationships between network concepts in networks
#'  that are approximately factorizable.
#'
#' \bold{Type IV}: eigengene-based (also known as eigennode-based) network
#' concepts are functions of the eigengene-based adjacency matrix
#' A.E=ConformityE*t(ConformityE) (diagonal included) and/or the corresponding
#' eigengene-based gene significance measure GSE. These network concepts can
#' only be defined for correlation networks. Details: The columns (nodes) of
#' datExpr can be summarized with the first principal component, which is
#' referred to as Eigengene in coexpression network analysis. In general
#' correlation networks, it is called eigennode. The eigengene-based conformity
#' ConformityE[i] is defined as abs(cor(datE[,i], Eigengene))^power where the
#' power corresponds to the power used for defining the weighted adjacency
#' matrix A. The eigengene-based conformity can also be used to define an
#' eigengene-based adjacency matrix A.E=ConformityE*t(ConformityE). The
#' eigengene based factorizability EF(datE) is a number between 0 and 1 that
#' measures how well A.E approximates A when the power parameter equals 1.
#' EF(datE) is defined with respect to the singular values of datExpr. For a
#' trait based node significance measure GS=abs(cor(datE,trait))^power, one can
#' also define an eigengene-based node significance measure
#' GSE[i]=ConformityE[i]*EigengeneSignificance where the eigengene significance
#' abs(cor(Eigengene,trait))^power is defined as power of the absolute value of
#' the correlation between eigengene and trait. Eigengene-based network concepts
#' are very useful for providing a geometric interpretation of network concepts
#' and for deriving relationships between network concepts. For example, the
#' hub gene significance measure and its eigengene-based analog have been used
#' to characterize networks where highly connected hub genes are important with
#' regard to a trait based gene significance measure (Horvath and Dong 2008).
#' @return
#' A list with the following components:
#' @param Summary a data frame whose rows report network concepts that only
#' depend on the adjacency matrix. Density (mean adjacency), Centralization,
#' Heterogeneity (coefficient of variation of the connectivity), Mean
#' ClusterCoef, Mean Connectivity. The columns of the data frame report the 4
#' types of network concepts mentioned in the description: Fundamental concepts,
#' eigengene-based concepts, conformity-based concepts, and approximate
#' conformity-based concepts.
#' @param Size reports the network size, i.e. the number of nodes, which equals
#' the number of columns of the input data frame datExpr.
#' @param Factorizability a number between 0 and 1. The closer it is to 1, the
#' better the off-diagonal elements of the conformity based network A.CF
#' approximate those of A (according to the Frobenius norm).
#' @param Eigengene the first principal component of the standardized columns of
#' datExpr. The number of components of this vector equals the number of rows of datExpr.
#' @param VarExplained the proportion of variance explained by the first
#' principal component (the Eigengene). It is numerically different from the
#' eigengene based factorizability. While VarExplained is based on the squares
#' of the singular values of datExpr, the eigengene-based factorizability is
#' based on fourth powers of the singular values.
#' @param Conformity numerical vector giving the conformity. The number of
#' components of the conformity vector equals the number of columns in datExpr.
#' The conformity is often highly correlated with the vector of node
#' connectivities. The conformity is computed using an iterative algorithm for
#' maximizing the factorizability measure. The algorithm and related network
#' concepts are described in Dong and Horvath 2007.
#' @param ClusterCoef a numerical vector that reports the cluster coefficient
#' for each node. This fundamental network concept measures the cliquishness of
#' each node.
#' @param Connectivity a numerical vector that reports the connectivity (also
#' known as degree) of each node. This fundamental network concept is also known
#' as whole network connectivity. One can also define the scaled connectivity
#' K=Connectivity/max(Connectivity) which is used for computing the hub gene
#' significance.
#' @param MAR a numerical vector that reports the maximum adjacency ratio for
#' each node. MAR[i] equals 1 if all non-zero adjacencies between node i and the
#' remaining network nodes equal 1. This fundamental network concept is always 1
#' for nodes of an unweighted network. This is a useful measure for weighted
#' networks since it allows one to determine whether a node has high connectivity
#' because of many weak connections (small MAR) or because of strong (but few)
#' connections (high MAR), see Horvath and Dong 2008.
#' @param ConformityE a numerical vector that reports the eigengene based (aka
#' eigenenode based) conformity for the correlation network. The number of
#' components equals the number of columns of datExpr.
#' @param GS a numerical vector that encodes the node (gene) significance. The
#' i-th component equals the node significance of the i-th column of datExpr if
#' a sample trait was supplied to the function (input trait).
#' \code{GS[i]=abs(cor(datE[,i], trait, use="p"))^power}
#' @param GSE a numerical vector that reports the eigengene based gene
#' significance measure. Its i-th component is given by
#' \code{GSE[i]=ConformityE[i]*EigengeneSignificance} where the eigengene
#' \code{significance abs(cor(Eigengene,trait))^power} is defined as power of
#' the absolute value of the correlation between eigengene and trait.
#' @param Significance a data frame whose rows report network concepts that also
#' depend on the trait based node significance measure. The rows correspond to
#' network concepts and the columns correspond to the type of network concept
#' (fundamental versus eigengene based). The first row of the data frame reports
#' the network significance. The fundamental version of this network concepts is
#' the average gene significance=mean(GS). The eigengene based analog of this
#' concept is defined as mean(GSE). The second row reports the hub gene
#' significance which is defined as slope of the intercept only regression model
#' that regresses the gene significance on the scaled network connectivity K.
#' The third row reports the eigengene significance
#' \code{abs(cor(Eigengene,trait))^power}. More details can be found in Horvath
#' and Dong (2008).
#' @author
#' Jun Dong, Steve Horvath, Peter Langfelder
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC
#' Systems Biology 2007, 1:24
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @seealso
#' \code{\link{conformityBasedNetworkConcepts}} for approximate conformity-based
#' network concepts.
#' \code{\link{fundamentalNetworkConcepts}} for calculation of fundamental
#' network concepts only.
#' @export
networkConcepts = function(datExpr, power=1, trait=NULL, networkType = "unsigned")
{

  networkTypeC = charmatch(networkType, .networkTypes)
  if (is.na(networkTypeC))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  if(networkTypeC==1)
  {
        adj <- abs(cor(datExpr,use="p"))^power
  } else if (networkTypeC==2)
  {
  	adj <- abs((cor(datExpr,use="p")+1)/2)^power
  } else {
        cor = cor(datExpr,use="p")
        cor[cor < 0] = 0
  	adj <- cor^power
  }
  diag(adj)=0 # Therefore adj=A-I.

  ### Fundamental Network Concepts
  Size=dim(adj)[1]
  Connectivity=apply(adj, 2, sum) # Within Module Connectivities
  Density=sum(Connectivity)/(Size*(Size-1))
  Centralization=Size*(max(Connectivity)-mean(Connectivity))/((Size-1)*(Size-2))
  Heterogeneity=sqrt(Size*sum(Connectivity^2)/sum(Connectivity)^2-1)
  ClusterCoef=.ClusterCoef.fun(adj)

  fMAR=function(v) sum(v^2)/sum(v)
  MAR=apply(adj, 1, fMAR)
  #CONNECTIVITY=Connectivity/max(Connectivity)

  ### Conformity-Based Network Concepts
  ### Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24
  Conformity=.NPC.iterate(adj)$v1
  Factorizability=1- sum( (adj-outer(Conformity,Conformity)+ diag(Conformity^2))^2 )/sum(adj^2)
  Connectivity.CF=sum(Conformity)*Conformity-Conformity^2
  Density.CF=sum(Connectivity.CF)/(Size*(Size-1))
  Centralization.CF=Size*(max(Connectivity.CF)-mean(Connectivity.CF))/((Size-1)*(Size-2))
  Heterogeneity.CF=sqrt(Size*sum(Connectivity.CF^2)/sum(Connectivity.CF)^2-1)
  #ClusterCoef.CF=.ClusterCoef.fun(outer(Conformity,Conformity)-diag(Conformity^2) )
  ClusterCoef.CF=c(NA, Size)
  for(i in 1:Size )
     ClusterCoef.CF[i]=( sum(Conformity[-i]^2)^2 - sum(Conformity[-i]^4) )/
       ( sum(Conformity[-i])^2 - sum(Conformity[-i]^2) )

  ### Approximate Conformity-Based Network Concepts
  Connectivity.CF.App=sum(Conformity)*Conformity
  Density.CF.App=sum(Connectivity.CF.App)/(Size*(Size-1))
  Centralization.CF.App=Size*(max(Connectivity.CF.App)-mean(Connectivity.CF.App))/((Size-1)*(Size-2))
  Heterogeneity.CF.App=sqrt(Size*sum(Connectivity.CF.App^2)/sum(Connectivity.CF.App)^2-1)
  ClusterCoef.CF.App=(sum(Conformity^2)/sum(Conformity))^2

  ### Eigengene-based Network Concepts
  m1=moduleEigengenes(datExpr, colors = rep(1, Size))
  # Weighted Expression Conformity
  ConformityE=cor(datExpr,m1[[1]][,1],use="pairwise.complete.obs"); ConformityE=abs(ConformityE)^power;
  ConnectivityE=sum(ConformityE)*ConformityE; #Expression Connectivity
  DensityE=sum(ConnectivityE)/(Size*(Size-1)); #Expression Density
  CentralizationE=Size*(max(ConnectivityE)-mean(ConnectivityE))/((Size-1)*(Size-2)); #Expression Centralization
  HeterogeneityE=sqrt(Size*sum(ConnectivityE^2)/sum(ConnectivityE)^2-1); #Expression Heterogeneity
  ClusterCoefE=(sum(ConformityE^2)/sum(ConformityE))^2; ##Expression ClusterCoef
  MARE=ConformityE* sum(ConformityE^2)/sum(ConformityE)

  ### Significance measure only when trait is available.
  if(!is.null(trait)){
    EigengeneSignificance = abs(cor(trait, m1[[1]], use="pairwise.complete.obs") )^power
    EigengeneSignificance = EigengeneSignificance[1,1]
    GS= abs(cor(datExpr, trait, use="pairwise.complete.obs") )^power; GS=GS[,1]
    GSE=ConformityE * EigengeneSignificance; GSE=GSE[,1]
    ModuleSignificance=mean(GS)
    ModuleSignificanceE=mean(GSE)
    K=Connectivity/max(Connectivity)
    HubGeneSignificance=sum(GS*K)/sum(K^2)
    KE=ConnectivityE/max(ConnectivityE)
    HubGeneSignificanceE= sum(GSE*KE)/sum(KE^2)
  }

  Summary=cbind(
    c(Density, Centralization, Heterogeneity, mean(ClusterCoef), mean(Connectivity)),
    c(DensityE, CentralizationE, HeterogeneityE, mean(ClusterCoefE), mean(ConnectivityE)),
    c(Density.CF, Centralization.CF, Heterogeneity.CF, mean(ClusterCoef.CF), mean(Connectivity.CF)),
    c(Density.CF.App, Centralization.CF.App, Heterogeneity.CF.App, mean(ClusterCoef.CF.App),
mean(Connectivity.CF.App) ) )
  colnames(Summary)=c("Fundamental", "Eigengene-based", "Conformity-Based", "Approximate Conformity-based")
  rownames(Summary)=c("Density", "Centralization", "Heterogeneity", "Mean ClusterCoef", "Mean Connectivity")

  output=list(Summary=Summary, Size=Size, Factorizability=Factorizability, Eigengene=m1[[1]],
VarExplained=m1[[2]][,1], Conformity=Conformity, ClusterCoef=ClusterCoef, Connectivity=Connectivity,
MAR=MAR, ConformityE=ConformityE)
  if(!is.null(trait)){
    output$GS=GS; output$GSE=GSE
    Significance=cbind(c(ModuleSignificance, HubGeneSignificance, EigengeneSignificance),
    c(ModuleSignificanceE, HubGeneSignificanceE, NA))
    colnames(Significance)=c("Fundamental", "Eigengene-based")
    rownames(Significance)=c("ModuleSignificance", "HubGeneSignificance", "EigengeneSignificance")
    output$Significance=Significance
  }
	output
}




#
# Network functions for network concepts
#

# Function definitions



#   Cohesiveness/Conformity/Factorizability etc
# Check if adj is a valid adjacency matrix:  square matrix, non-negative entries, symmetric and no missing entries.
# Parameters:
#   adj - the input adjacency matrix
#   tol - the tolerence level to measure the difference from 0 (symmetric matrix: upper diagonal minus lower diagonal)
# Remarks:
#   1. This function is not supposed to be used directly. Instead, it should appear in function definitions.
#   2. We release the requirement that the diagonal elements be 1 or 0. Users should assign appropriate values
#      at the beginning of their function definitions.
# Usage:
# if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")

.is.adjmat = function(adj, tol=10^(-15)){
	n=dim(adj)
	is.adj=1
	if (n[1] != n[2]){ message("The adjacency matrix is not a square matrix!"); is.adj=0;}
	if ( sum(is.na(adj))>0 ){ message("There are missing values in the adjacency matrix!"); is.adj=0;}
	if ( sum(adj<0)>0 ){ message("There are negative entries in the adjacency matrix!"); is.adj=0;}
	if ( max(abs(adj-t(adj))) > tol){ message("The adjacency matrix is not symmetric!"); is.adj=0;}
			#if ( max(abs(diag(adj)-1)) > tol){ message("The diagonal elements are not all one!"); is.adj=0;}
			#The last criteria is removed because of different definitions on diagonals with other papers.
			#Always let "diagonal=1" INSIDE the function calls when using functions for Factorizability paper.
	is.adj
}



# .NPC.direct=function(adj)
# Calculates the square root of Normalized Product Connectivity (.NPC), by way of definition of .NPC. ( \sqrt{t})
# Parameters:
#   adj - the input adjacency matrix
#   tol - the tolerence level to measure the difference from 0 (zero off-diagonal elements)
# Output:
#   v1 - vector, the square root of .NPC
# Remarks:
#   1. The function requires that the off-diagonal elements of the adjacency matrix are all non-zero.
#   2. If any of the off-diagonal elements is zero, use the function .NPC.iterate().
#   3. If the adjacency matrix is 2 by 2, then a warning message is issued and vector of sqrt(adj[1,2]) is returned.
#   4. If the adjacency matrix is a ZERO matrix, then a warning message is issued and vector of 0 is returned.

.NPC.direct=function(adj){
	if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
	n=dim(adj)[1]
	if(n==2) {
		warning("The adjacecny matrix is only 2 by 2. .NPC may not be unique!")
		return(rep(sqrt(adj[1,2]),2))
	}
	diag(adj)=0
	if(!sum(adj>0)){
	  warning("The adjacency matrix is a ZERO matrix!")
	  return(rep(0,n))
	}
	diag(adj)=1
	if(sum(adj==0)) stop("There is zero off--diagonal element! Please use the function .NPC.iterate().")
	log10.prod.vec=function(vec){
		prod=0
		for(i in 1:length(vec) )
			prod=prod+log10(vec[i])
		prod
	}
	off.diag=as.vector(as.dist(adj))
	prod1=log10.prod.vec(off.diag)
	v1=rep(-666, n)
	for(i in 1:n){
		prod2=prod1-log10.prod.vec(adj[i,])
		v1[i]=10^(prod1/(n-1)-prod2/(n-2))
	}
	v1
}


# .NPC.iterate=function(adj, loop=10^(10), tol=10^(-10))
# Calculates the square root of Normalized Product Connectivity, by way of iteration algorithm. ( \sqrt{t})
# Parameters:
#   adj - the input adjacency matrix
#   loop - the maximum number of iterations before stopping the algorithm
#   tol - the tolerence level to measure the difference from 0 (zero off-diagonal elements)
# Output:
#   v1 - vector, the square root of .NPC
#   loop - integer, the number of iterations taken before convergence criterion is met
#   diff - scaler, the maximum difference between the estimates of 'v1' in the last two iterations
# Remarks:
#   1. Whenever possible, use .NPC.direct().
#   2. If the adjacency matrix is 2 by 2, then a warning message is issued.
#   3. If the adjacency matrix is a ZERO matrix, then a warning message is issued and vector of 0 is returned.

if( exists(".NPC.iterate") ) rm(.NPC.iterate)
.NPC.iterate=function(adj, loop=10^(10), tol=10^(-10)){
	if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
	n=dim(adj)[1]
	if(n==2) warning("The adjacecny matrix is only 2 by 2. .NPC may not be unique!")
	diag(adj)=0
	if(max(abs(adj))<tol){
	  warning("The adjacency matrix is a ZERO matrix!")
	  return(rep(0,n))
	}
	diff=1
	k=apply(adj, 2, sum)
	v1=k/sqrt(sum(k)) # first-step estimator for v-vector (initial value)
	i=0
	while( loop>i && diff>tol ){
	i=i+1
	diag(adj)=v1^2
	svd1=svd(adj) # Spectral Decomposition
	v2=sqrt(svd1$d[1])*abs(svd1$u[,1])
	diff=max(abs(v1-v2))
	v1=v2
	}
	list(v1=v1,loop=i,diff=diff)
}


# The function .ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix
.ClusterCoef.fun=function(adjmat1)
{
  # diag(adjmat1)=0
  no.nodes=dim(adjmat1)[[1]]
  computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
  computeSqDiagSum = function(x, vec) { sum(x^2 * vec) }
  nolinksNeighbors <- c(rep(-666,no.nodes))
  total.edge <- c(rep(-666,no.nodes))
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
  if (maxh1>1 | minh1 < 0 )
  {
     stop(paste("ERROR: the adjacency matrix contains entries that are larger",
                 "than 1 or smaller than 0: max=",maxh1,", min=",minh1))
  } else {
    nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
    subTerm = apply(adjmat1, 1, computeSqDiagSum, vec = diag(adjmat1))
    plainsum  <- apply(adjmat1, 1, sum)
    squaresum <- apply(adjmat1^2, 1, sum)
    total.edge = plainsum^2 - squaresum
    CChelp=rep(-666, no.nodes)
    CChelp=ifelse(total.edge==0,0, (nolinksNeighbors-subTerm)/total.edge)
    CChelp
  }
} # end of function



# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)
.err.bp<-function(daten,error,two.side=F)
{
   if(!is.numeric(daten)) {
        stop("All arguments must be numeric")}
   if(is.vector(daten)){
      xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1))))
   }else{
      if (is.matrix(daten)){
        xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
                            dim=c(1,length(daten))))+0:(length(daten)-1)+.5
      }else{
        stop("First argument must either be a vector or a matrix") }
   }
   MW<-0.25*(max(xval)/length(xval))
   ERR1<-daten+error
   ERR2<-daten-error
   for(i in 1:length(daten)){
      segments(xval[i],daten[i],xval[i],ERR1[i])
      segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
      if(two.side){
        segments(xval[i],daten[i],xval[i],ERR2[i])
        segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
      }
   }
}




#' @name conformityBasedNetworkConcepts
#' @rdname conformityBasedNetworkConcepts
#' @title calculation of conformity-based network concepts
#' @description
#' This function computes 3 types of network concepts (also known as network
#' indices or statistics) based on an adjacency matrix and optionally a node
#' significance measure.
#' @param adj adjacency matrix. A symmetric matrix with components between 0 and
#' 1.
#' @param GS optional node significance measure. A vector with length equal the
#' dimension of adj.
#' @details
#' This function computes 3 types of network concepts (also known as network
#' indices or statistics) based on an adjacency matrix and optionally a node
#' significance measure. Specifically, it computes I) fundamental network
#' concepts, II) conformity based network concepts, and III) approximate
#' conformity based network concepts. These network concepts are defined for any
#' symmetric adjacency matrix (weighted and unweighted). The network concepts
#' are described in Dong and Horvath (2007) and Horvath and Dong (2008). In the
#' following, we use the term gene and node interchangeably since these methods
#' were originally developed for gene networks. In the following, we briefly
#' describe the 3 types of network concepts:
#'
#' \bold{Type I}: fundamental network concepts are defined as a function of the
#' off-diagonal elements of an adjacency matrix A and/or a node significance
#' measure GS.
#'
#' \bold{Type II}: conformity-based network concepts are functions of the
#' off-diagonal elements of the conformity based adjacency matrix A.CF=CF*t(CF)
#' and/or the node significance measure. These network concepts are defined for
#' any network for which a conformity vector can be defined. Details: For any
#' adjacency matrix A, the conformity vector CF is calculated by requiring that
#' A[i,j] is approximately equal to CF[i]*CF[j]. Using the conformity one can
#' define the matrix A.CF=CF*t(CF) which is the outer product of the conformity
#' vector with itself. In general, A.CF is not an adjacency matrix since its
#' diagonal elements are different from 1. If the off-diagonal elements of A.CF
#' are similar to those of A according to the Frobenius matrix norm, then A is
#' approximately factorizable. To measure the factorizability of a network, one
#' can calculate the Factorizability, which is a number between 0 and 1 (Dong
#' and Horvath 2007). The conformity is defined using a monotonic, iterative
#' algorithm that maximizes the factorizability measure.
#'
#' \bold{Type III}: approximate conformity based network concepts are functions
#' of all elements of the conformity based adjacency matrix A.CF (including the
#' diagonal) and/or the node significance measure GS. These network concepts are
#' very useful for deriving relationships between network concepts in networks
#' that are approximately factorizable.
#' @return
#' A list with the following components:
#' @param Factorizability number between 0 and 1 giving the factorizability of
#' the matrix. The closer to 1 the higher the evidence of factorizability, that
#' is, A-I is close to outer(CF,CF)-diag(CF^2).
#' @param fundamentalNCs fundamental network concepts, that is network concepts
#' calculated directly from the given adjacency matrix adj. A list with
#' components ScaledConnectivity (giving the scaled connectivity of each node),
#' Connectivity (connectivity of each node), ClusterCoef (the clustering
#' coefficient of each node), MAR (maximum adjacency ratio of each node),
#' Density (the mean density of the network), Centralization (the centralization
#' of the network), Heterogeneity (the heterogeneity of the network). If the
#' input node significance GS is specified, the following additional components
#' are included: NetworkSignificance (network significance, the mean node
#' significance), and HubNodeSignificance (hub node significance given by the
#' linear regression of node significance on connectivity).
#' @param conformityBasedNCs network concepts based on an approximate adjacency
#' matrix given by the outer product of the conformity vector but with unit
#' diagonal. A list with components Conformity (the conformity vector) and
#' Connectivity.CF, ClusterCoef.CF, MAR.CF, Density.CF, Centralization.CF,
#' Heterogeneity.CF giving the conformity-based analogs of the above network
#' concepts.
#' @param approximateConformityBasedNCs network concepts based on an approximate
#' adjacency matrix given by the outer product of the conformity vector. A list
#' with components Conformity (the conformity vector) and Connectivity.CF.App,
#' ClusterCoef.CF.App, MAR.CF.App, Density.CF.App, Centralization.CF.App,
#' Heterogeneity.CF.App giving the conformity-based analogs of the above network
#' concepts.
#' @author
#' Steve Horvath
#' @references
#' Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC
#' Systems Biology 2007, 1:24 Horvath S, Dong J (2008) Geometric Interpretation
#' of Gene Coexpression Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @seealso
#' \code{\link{networkConcepts}}for calculation of eigennode based network
#' concepts for a correlation network.
#' \code{\link{fundamentalNetworkConcepts}} for calculcation of fundamental
#' network concepts only-
conformityBasedNetworkConcepts = function(adj, GS=NULL)
{
  if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
  diag(adj)=0 # Therefore adj=A-I.
  if (dim(adj)[[1]]<3)
    stop("The adjacency matrix has fewer than 3 rows. This network is trivial and will not be evaluated.")
  if (!is.null(GS))
  {
    if( length(GS) !=dim(adj)[[1]])
    {
       stop(paste("The length of the node significnce GS does not equal the number",
                  "of rows of the adjcency matrix. length(GS) != dim(adj)[[1]]. \n",
                  "Something is wrong with your input"))
    }
  }

	### Fundamental Network Concepts
	Size=dim(adj)[1]
	Connectivity=apply(adj, 2, sum)
	Density=sum(Connectivity)/(Size*(Size-1))
	Centralization=Size*(max(Connectivity)-mean(Connectivity))/((Size-1)*(Size-2))
	Heterogeneity=sqrt(Size*sum(Connectivity^2)/sum(Connectivity)^2-1)
	ClusterCoef=.ClusterCoef.fun(adj)
	fMAR=function(v) sum(v^2)/sum(v)
	MAR=apply(adj, 1, fMAR)
	### Conformity-Based Network Concepts
	Conformity=.NPC.iterate(adj)$v1
	Factorizability=1- sum( (adj-outer(Conformity,Conformity)+ diag(Conformity^2))^2 )/sum(adj^2)
	Connectivity.CF=sum(Conformity)*Conformity-Conformity^2
	Density.CF=sum(Connectivity.CF)/(Size*(Size-1))
	Centralization.CF=Size*(max(Connectivity.CF)-mean(Connectivity.CF))/((Size-1)*(Size-2))
	Heterogeneity.CF=sqrt(Size*sum(Connectivity.CF^2)/sum(Connectivity.CF)^2-1)
	#ClusterCoef.CF=.ClusterCoef.fun(outer(Conformity,Conformity)-diag(Conformity^2) )
	ClusterCoef.CF=c(NA, Size)
  for(i in 1:Size )
    ClusterCoef.CF[i]=( sum(Conformity[-i]^2)^2 - sum(Conformity[-i]^4) )/
                          ( sum(Conformity[-i])^2 - sum(Conformity[-i]^2) )
  MAR.CF=ifelse(sum(Conformity,na.rm=T)-Conformity==0, NA,
                Conformity*(sum(Conformity^2,na.rm=T)-Conformity^2)/(sum(Conformity,na.rm=T)-Conformity))

	### Approximate Conformity-Based Network Concepts
	Connectivity.CF.App=sum(Conformity)*Conformity
	Density.CF.App=sum(Connectivity.CF.App)/(Size*(Size-1))
	Centralization.CF.App=Size*(max(Connectivity.CF.App)-mean(Connectivity.CF.App))/((Size-1)*(Size-2))
	Heterogeneity.CF.App=sqrt(Size*sum(Connectivity.CF.App^2)/sum(Connectivity.CF.App)^2-1)

  if(sum(Conformity,na.rm=T)==0)
  {
      warning(paste("The sum of conformities equals zero.\n",
                    "Maybe you used an input adjacency matrix with lots of zeroes?\n",
                    "Specifically, sum(Conformity,na.rm=T)==0."))
      MAR.CF.App= rep(NA,Size)
      ClusterCoef.CF.App= rep(NA,Size)
  } #end of if
  if(sum(Conformity,na.rm=T) !=0)
  {
    MAR.CF.App=Conformity*sum(Conformity^2,na.rm=T) /sum(Conformity,na.rm=T)
    ClusterCoef.CF.App=rep((sum(Conformity^2)/sum(Conformity))^2,Size)
  }# end of if
  output=list(
              Factorizability =Factorizability,
              fundamentalNCs=list(
                  ScaledConnectivity=Connectivity/max(Connectivity,na.rm=T),
                  Connectivity=Connectivity,
                  ClusterCoef=ClusterCoef,
                  MAR=MAR,
                  Density=Density,
                  Centralization =Centralization,
                  Heterogeneity= Heterogeneity),
              conformityBasedNCs=list(
                  Conformity=Conformity,
                  Connectivity.CF=Connectivity.CF,
                  ClusterCoef.CF=ClusterCoef.CF,
                  MAR.CF=MAR.CF,
                  Density.CF=Density.CF,
                  Centralization.CF =Centralization.CF,
                  Heterogeneity.CF= Heterogeneity.CF),
              approximateConformityBasedNCs=list(
                  Conformity=Conformity,
                  Connectivity.CF.App= Connectivity.CF.App,
                  ClusterCoef.CF.App=ClusterCoef.CF.App,
                  MAR.CF.App=MAR.CF.App,
                  Density.CF.App= Density.CF.App,
                  Centralization.CF.App =Centralization.CF.App,
                  Heterogeneity.CF.App= Heterogeneity.CF.App))
  if ( !is.null(GS) )
  {
    output$FundamentalNC$NetworkSignificance = mean(GS,na.rm=T)
    K = Connectivity/max(Connectivity)
    output$FundamentalNC$HubNodeSignificance = sum(GS * K,na.rm=T)/sum(K^2,na.rm=T)
  }
  output
} # end of function




#' @name fundamentalNetworkConcepts
#' @rdname fundamentalNetworkConcepts
#' @title Calculation of fundamental network concepts from an adjacency matrix
#' @description
#' This function computes fundamental network concepts (also known as network
#' indices or statistics) based on an adjacency matrix and optionally a node
#' significance measure. These network concepts are defined for any symmetric
#' adjacency matrix (weighted and unweighted). The network concepts are
#' described in Dong and Horvath (2007) and Horvath and Dong (2008).
#' Fundamental network concepts are defined as a function of the off-diagonal
#' elements of an adjacency matrix adj and/or a node significance measure GS.
#' @inheritParams conformityBasedNetworkConcepts
#' @return
#' A list of elements:
#' @param Connectivity a  numerical vector that reports the connectivity (also
#' known as degree) of each node. This fundamental network concept is also known
#' as whole network connectivity. One can also define the scaled connectivity
#' K=Connectivity/max(Connectivity) which is used for computing the hub gene
#' significance.
#' @param ScaledConncectivity the Connectivity vector scaled by the highest
#' connectivity in the network, i.e., Connectivity/max(connectivity).
#' @param ClusterCoef a numerical vector that reports the cluster coefficient
#' for each node. This fundamental network concept measures the cliquishness of
#' each node
#' @param MAR a numerical vector that reports the maximum adjacency ratio for
#' each node. MAR[i] equals 1 if all non-zero adjacencies between node i and the
#' remaining network nodes equal 1. This fundamental network concept is always 1
#' for nodes of an unweighted network. This is a useful measure for weighted
#' networks since it allows one to determine whether a node has high connectivity
#' because of many weak connections (small MAR) or because of strong (but few)
#' connections (high MAR), see Horvath and Dong 2008.
#' @param Density the density of the network
#' @param Centralization the centralization of the network
#' @param Heterogeneity the heterogeneity of the network
#' @author
#' Steve Horvath
#' @references
#' Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC
#' Systems Biology 2007, 1:24
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @seealso
#' \code{\link{conformityBasedNetworkConcepts}} for calculation of conformity
#' based network concepts for a network adjacency matrix;
#' \code{\link{networkConcepts}}, for calculation of conformity based and
#' eigennode based network concepts for a correlation network.
fundamentalNetworkConcepts=function(adj,GS=NULL)
{
   if(!.is.adjmat(adj)) stop("The input matrix is not a valid adjacency matrix!")
   diag(adj)=0 # Therefore adj=A-I.

   if (dim(adj)[[1]]<3) stop("The adjacency matrix has fewer than 3 rows. This network is trivial and will not be evaluated.")

if (!is.null(GS)) { if( length(GS) !=dim(adj)[[1]]){ stop("The length of the node significnce GS does not
equal the number of rows of the adjcency matrix. length(GS) unequal dim(adj)[[1]]. GS should be a vector whosecomponents correspond to the nodes.")}}

	Size=dim(adj)[1]

### Fundamental Network Concepts
	Connectivity=apply(adj, 2, sum) # Within Module Connectivities
	Density=sum(Connectivity)/(Size*(Size-1))
	Centralization=Size*(max(Connectivity)-mean(Connectivity))/((Size-1)*(Size-2))
	Heterogeneity=sqrt(Size*sum(Connectivity^2)/sum(Connectivity)^2-1)
	ClusterCoef=.ClusterCoef.fun(adj)
	fMAR=function(v) sum(v^2)/sum(v)
	MAR=apply(adj, 1, fMAR)
	ScaledConnectivity=Connectivity/max(Connectivity,na.rm=T)

	output=list(
                  Connectivity=Connectivity,
                  ScaledConnectivity=ScaledConnectivity,
                  ClusterCoef=ClusterCoef, MAR=MAR,
                  Density=Density, Centralization =Centralization,
                  Heterogeneity= Heterogeneity)

   if ( !is.null(GS) ) {
        output$NetworkSignificance = mean(GS,na.rm=T)
        output$HubNodeSignificance = sum(GS * ScaledConnectivity,na.rm=T)/sum(ScaledConnectivity^2,na.rm=T)
   }
   output
} # end of function
