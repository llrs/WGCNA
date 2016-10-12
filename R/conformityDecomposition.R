#


#' Conformity and module based decomposition of a network adjacency matrix.
#'
#' The function calculates the conformity based approximation \code{A.CF} of an
#' adjacency matrix and a factorizability measure codeFactorizability. If a
#' module assignment \code{Cl} is provided, it also estimates a corresponding
#' intermodular adjacency matrix. In this case, function automatically carries
#' out the module- and conformity based decomposition of the adjacency matrix
#' described in chapter 2 of (Horvath 2011).
#'
#' We distinguish two situation depending on whether or not \code{Cl} equals
#' \code{NULL}.  1) Let us start out assuming that \code{Cl = NULL}. In this
#' case, the function calculates the conformity vector for a general, possibly
#' non-factorizable network \code{adj} by minimizing a quadratic (sums of
#' squares) loss function. The conformity and factorizability for an adjacency
#' matrix is defined in (Dong and Horvath 2007, Horvath and Dong 2008) but we
#' briefly describe it in the following. A network is called exactly
#' factorizable if the pairwise connection strength (adjacency) between 2
#' network nodes can be factored into node specific contributions, named node
#' 'conformity', i.e. if \code{adj[i,j]=Conformity[i]*Conformity[j]}. The
#' conformity turns out to be highly related to the network connectivity (aka
#' degree). If \code{adj} is not exactly factorizable, then the function
#' \code{conformityDecomposition} calculates a conformity vector of the exactly
#' factorizable network that best approximates \code{adj}. The factorizability
#' measure \code{Factorizability} is a number between 0 and 1. The higher
#' \code{Factorizability}, the more factorizable is \code{adj}. Warning: the
#' algorithm may only converge to a local optimum and it may not converge at
#' all. Also see the notes below.
#'
#' 2) Let us now assume that \code{Cl} is not NULL, i.e. it specifies the
#' module assignment of each node. Then the function calculates a module- and
#' CF-based approximation of \code{adj} (explained in chapter 2 in Horvath
#' 2011). In this case, the function calculates a conformity vector
#' \code{Conformity} and a matrix \code{IntermodularAdjacency} such that
#' \code{adj[i,j]} is approximately equal to
#' \code{Conformity[i]*Conformity[j]*IntermodularAdjacency[module.index[i],module.index[j]]}
#' where \code{module.index[i]} is the row of the matrix
#' \code{IntermodularAdjacency} that corresponds to the module assigned to node
#' i. To estimate \code{Conformity} and a matrix \code{IntermodularAdjacency},
#' the function attempts to minimize a quadratic loss function (sums of
#' squares). Currently, the function only implements a heuristic algorithm for
#' optimizing the objective function (chapter 2 of Horvath 2011). Another, more
#' accurate Majorization Minorization (MM) algorithm for the decomposition is
#' implemented in the function \code{propensityDecomposition} by Ranola et al
#' (2011).
#'
#' @param adj a symmetric numeric matrix (or data frame) whose entries lie
#' between 0 and 1.
#' @param Cl a vector (or factor variable) of length equal to the number of
#' rows of \code{adj}. The variable assigns each network node (row of
#' \code{adj}) to a module. The entries of \code{Cl} could be integers or
#' character strings.
#' @return \item{A.CF}{a symmetric matrix that approximates the input matrix
#' \code{adj}. Roughly speaking, the i,j-the element of the matrix equals
#' \code{Conformity[i]*Conformity[j]*IntermodularAdjacency[module.index[i],module.index[j]]}
#' where \code{module.index[i]} is the row of the matrix
#' \code{IntermodularAdjacency} that corresponds to the module assigned to node
#' i. } \item{Conformity}{a numeric vector whose entries correspond to the rows
#' of codeadj. If \code{Cl=NULL} then \code{Conformity[i]} is the conformity.
#' If \code{Cl} is not NULL then \code{Conformity[i]} is the intramodular
#' conformity with respect to the module that node i belongs to. }
#' \item{IntermodularAdjacency}{ a symmetric matrix (data frame) whose rows and
#' columns correspond to the number of modules specified in \code{Cl}.
#' Interpretation: it measures the similarity (adjacency) between the modules.
#' In this case, the rows (and columns) of \code{IntermodularAdjacency}
#' correspond to the entries of \code{Cl.level}. } \item{Factorizability}{ is a
#' number between 0 and 1. If \code{Cl=NULL} then it equals 1, if (and only if)
#' \code{adj} is exactly factorizable. If \code{Cl} is a vector, then it
#' measures how well the module- and CF based decomposition approximates
#' \code{adj}.  } \item{Cl.level}{ is a vector of character strings which
#' correspond to the factor levels of the module assignment \code{Cl}.
#' Incidentally, the function automatically turns \code{Cl} into a factor
#' variable. The components of Conformity and
#' \code{IntramodularFactorizability} correspond to the entries of
#' \code{Cl.level}. } \item{IntramodularFactorizability}{ is a numeric vector
#' of length equal to the number of modules specified by \code{Cl}. Its entries
#' report the factorizability measure for each module. The components
#' correspond to the entries of \code{Cl.level}.} \item{listConformity}{}
#' @note Regarding the situation when \code{Cl=NULL}. One can easily show that
#' the conformity vector is not unique if \code{adj} contains only 2 nodes.
#' However, for more than 2 nodes the conformity is uniquely defined when
#' dealing with an exactly factorizable weighted network whose entries
#' \code{adj[i,j]} are larger than 0. In this case, one can get explicit
#' formulas for the conformity (Dong and Horvath 2007).
#' @author Steve Horvath
#' @seealso
#' \code{\link{conformityBasedNetworkConcepts}}
#' \code{\link{propensityDecomposition}}
#' @references Dong J, Horvath S (2007) Understanding Network Concepts in
#' Modules. BMC Systems Biology 2007, June 1:24 Horvath S, Dong J (2008)
#' Geometric Interpretation of Gene Co-Expression Network Analysis. PloS
#' Computational Biology. 4(8): e1000117. PMID: 18704157 Horvath S (2011)
#' Weighted Network Analysis. Applications in Genomics and Systems Biology.
#' Springer Book. ISBN: 978-1-4419-8818-8 Ranola JMO, Langfelder P, Song L,
#' Horvath S, Lange K (2011) An MM algorithm for the module- and propensity
#' based decomposition of a network. Currently a draft.
#' @keywords misc
#' @examples
#'
#'
#' # assume the number of nodes can be divided by 2 and by 3
#' n=6
#' # here is a perfectly factorizable matrix
#' A=matrix(1,nrow=n,ncol=n)
#' # this provides the conformity vector and factorizability measure
#' conformityDecomposition(adj=A)
#' # now assume we have a class assignment
#' Cl=rep(c(1,2),c(n/2,n/2))
#' conformityDecomposition(adj=A,Cl=Cl)
#' # here is a block diagonal matrix
#' blockdiag.A=A
#' blockdiag.A[1:(n/3),(n/3+1):n]=0
#' blockdiag.A[(n/3+1):n , 1:(n/3)]=0
#' block.Cl=rep(c(1,2),c(n/3,2*n/3))
#' conformityDecomposition(adj= blockdiag.A,Cl=block.Cl)
#'
#' # another block diagonal matrix
#' blockdiag.A=A
#' blockdiag.A[1:(n/3),(n/3+1):n]=0.3
#' blockdiag.A[(n/3+1):n , 1:(n/3)]=0.3
#' block.Cl=rep(c(1,2),c(n/3,2*n/3))
#' conformityDecomposition(adj= blockdiag.A,Cl=block.Cl)
#'
#'
#' @export conformityDecomposition
conformityDecomposition = function (adj, Cl = NULL)
{
    if (  is.null(dim(adj) )) stop("Input adj is not a matrix or data frame. ")
    if ( dim(adj)[[1]] < 3)   stop("The adjacency matrix has fewer than 3 rows. This network is trivial and will not be evaluated.")
    if (!.is.adjmat(adj))
        stop("The input matrix is not a valid adjacency matrix!")
    diag(adj) = 0
    if (!is.null(Cl)) {
        if (length(Cl) != dim(adj)[[1]]) {
            stop(paste("The length of the class assignment Cl does not equal the number",
                       "of rows of the adjcency matrix. length(Cl) != dim(adj)[[1]]. \n",
                       "Something is wrong with your input"))
        }
        if (sum(is.na(Cl))>0 ) stop("Cl must not contain missing values (NA)." )
    }

    A.CF=matrix(0, nrow=dim(adj)[[1]], ncol= dim(adj)[[2]] )
    diag(A.CF)=1
    if ( is.null(Cl) )  {
        Conformity = .NPC.iterate(adj)$v1
        if (sum(adj^2,na.rm=T)==0) {Factorizability=NA} else {
            A.CF=outer(Conformity, Conformity) - diag(Conformity^2)
            Factorizability = 1 - sum((adj - A.CF)^2)/sum(adj^2)}
        diag(A.CF)=1
        output = list(A.CF=A.CF, Conformity=data.frame( Conformity),  IntermodularAdjacency=1, Factorizability = Factorizability)
    }

    if ( !is.null(Cl) )  {
        Cl=factor(Cl)
        Cl.level=levels( Cl )
        if ( length(Cl.level)>100 ) warning(paste("Your class assignment variable Cl contains",  length(Cl.level), "different classes. I assume this is a proper class assignment variable. But if not, stop the calculation, e.g. by using the Esc key on your keybord."))

        Conformity=rep(NA, length(Cl) )
        listConformity=list()
        IntramodularFactorizability=rep(NA, length(Cl.level) )
        IntermodularAdjacency= matrix(0,nrow=length(Cl.level),ncol=length(Cl.level) )
        diag(IntermodularAdjacency)=1
        IntermodularAdjacency=data.frame(IntermodularAdjacency)
        dimnames(IntermodularAdjacency)[[1]]=as.character(Cl.level)
        dimnames(IntermodularAdjacency)[[2]]=as.character(Cl.level)

        numeratorFactorizability=0
        for (i in 1:length(Cl.level) ) {
            restclass= Cl== Cl.level[i]
            if (sum(restclass)==1) { A.help=0; CF.help =0;   Conformity[restclass]=CF.help;
            A.CF[restclass,restclass]=CF.help*CF.help - CF.help^2
            }
            if (sum(restclass)==2) {
                A.help=adj[restclass,restclass];diag(A.help)=0
                CFvalue=sqrt(adj[restclass,restclass][1,2]);
                CF.help= c(CFvalue , CFvalue )
                Conformity[restclass]=CF.help
                A.CF[restclass,restclass]=outer(CF.help, CF.help) - diag(CF.help^2)
            }

            if (sum(restclass)>2) {
                A.help=adj[restclass,restclass];diag(A.help)=0 ;
                CF.help = .NPC.iterate(A.help )$v1
                Conformity[restclass]=CF.help
                A.CF[restclass,restclass]=outer(CF.help, CF.help) - diag(CF.help^2)
            }
            if (length(CF.help)>1) {numeratorFactorizability= numeratorFactorizability+sum(   (A.help-outer(CF.help,CF.help)+diag(CF.help^2) )^2 )}

            listConformity[[i]] =CF.help
            if (sum(A.help^2,na.rm=T)==0 |  length(CF.help)==1 ) {IntramodularFactorizability[i]=NA} else {
                IntramodularFactorizability[i] = 1 - sum((A.help - outer(CF.help, CF.help) +
                                                              diag(CF.help^2))^2)/sum(A.help^2,na.rm=T)  }
        } # end of for loop over i

        if ( length(Cl.level)==1) {IntermodularAdjacency[1,1]=1} else {
            for (i in 1:(length(Cl.level)-1)  ) {
                for (j in (i+1):length(Cl.level) ) {
                    restclass1= Cl== Cl.level[i]
                    restclass2= Cl== Cl.level[j]
                    A.inter=adj[restclass1,restclass2]
                    mean.CF1=mean(listConformity[[i]], na.rm=T)
                    mean.CF2=mean(listConformity[[j]], na.rm=T)
                    if (   mean.CF1* mean.CF2    != 0 ) {
                        IntermodularAdjacency[i,j]= mean(A.inter,na.rm=T)/(mean.CF1* mean.CF2)
                        IntermodularAdjacency[j,i]= IntermodularAdjacency[i,j]
                    }

                    if (  length(listConformity[[i]])==1 |  length(listConformity[[j]])==1    )   {
                        numeratorFactorizability=  numeratorFactorizability+ 2*sum( (A.inter- IntermodularAdjacency[i,j]* listConformity[[i]] * listConformity[[j]])^2  )
                        A.CF[restclass1,restclass2]= IntermodularAdjacency[i,j]* listConformity[[i]] *listConformity[[j]]
                        A.CF[restclass2,restclass1]= IntermodularAdjacency[j,i]* listConformity[[j]] *listConformity[[i]]
                    } else {
                        numeratorFactorizability=  numeratorFactorizability+ 2*sum( (A.inter- IntermodularAdjacency[i,j]*outer( listConformity[[i]] , listConformity[[j]]) )^2  )
                        A.CF[restclass1,restclass2]= IntermodularAdjacency[i,j]* outer(listConformity[[i]], listConformity[[j]]  )
                        A.CF[restclass2,restclass1]= IntermodularAdjacency[j,i]* outer(listConformity[[j]], listConformity[[i]]  )

                    } # end of else
                } # end of for (i in
            } # end of for (j in
            diag(adj)=NA
            Factorizability=  1-numeratorFactorizability/ sum(adj^2,na.rm=T)

        } # end of if else statement

        diag(A.CF)=1
        output = list(A.CF=A.CF,  Conformity=Conformity, IntermodularAdjacency= IntermodularAdjacency, Factorizability = Factorizability,  Cl.level =Cl.level, IntramodularFactorizability=  IntramodularFactorizability, listConformity=listConformity )
    } # end of   if ( !is.null(Cl) )
    output
}
