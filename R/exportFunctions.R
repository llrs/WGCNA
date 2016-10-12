# Functions for exporting networks to various network visualization software


#' Export network data in format readable by VisANT
#' 
#' Exports network data in a format readable and displayable by the VisANT
#' software.
#' 
#' The adjacency matrix is checked for validity. The entries can be negative,
#' however. The adjacency matrix is expected to also have valid \code{names} or
#' \code{dimnames[[2]]} that represent the probe names of the corresponding
#' edges.
#' 
#' Whether the output is a weighted network or not, only edges whose (absolute
#' value of) adjacency are above \code{threshold} will be included in the
#' output. If \code{maxNConnections} is given, at most \code{maxNConnections}
#' will be included in the output.
#' 
#' If \code{probeToGene} is given, it is expected to have two columns, the
#' first one corresponding to the probe names, the second to their
#' corresponding gene names that will be used in the output.
#' 
#' @param adjMat adjacency matrix of the network to be exported.
#' @param file character string specifying the file name of the file in which
#' the data should be written. If not given, no file will be created. The file
#' is in a plain text format.
#' @param weighted logical: should the exported network by weighted?
#' @param threshold adjacency threshold for including edges in the output.
#' @param maxNConnections maximum number of exported adjacency edges. This can
#' be used as another filter on the exported edges.
#' @param probeToGene optional specification of a conversion between probe
#' names (that label columns and rows of \code{adjacency}) and gene names (that
#' should label nodes in the output).
#' @return A data frame containing the network information suitable as input to
#' VisANT. The same data frame is also written into a file specified by
#' \code{file}, if given.
#' @author Peter Langfelder
#' @references VisANT software is available from http://visant.bu.edu/.
#' @keywords misc
#' @export exportNetworkToVisANT
exportNetworkToVisANT = function(
  adjMat,
  file = NULL,
  weighted = TRUE,
  threshold = 0.5,
  maxNConnections = NULL,
  probeToGene = NULL)
{
  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0
  nRow = nrow(adjMat)
  checkAdjMat(adjMat, min = -1, max = 1)
  probes = dimnames(adjMat)[[1]]
  if (!is.null(probeToGene))
  {
    probes2genes = match(probes, probeToGene[,1])
    if (sum(is.na(probes2genes)) > 0)
      stop("Error translating probe names to gene names: some probe names could not be translated.")
    probes = probeToGene[probes2genes, 2]
  }

  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
  colMat = matrix(c(1:nRow), nRow, nRow)

  adjDst = as.dist(adjMat)
  dstRows = as.dist(rowMat)
  dstCols = as.dist(colMat)

  if (is.null(maxNConnections)) maxNConnections = length(adjDst)

  ranks = rank(-abs(adjDst), na.last = TRUE, ties.method = "first")
  edges = abs(adjDst) > threshold & ranks <= maxNConnections
  nEdges = sum(edges)

  visAntData = data.frame (
     from = probes[dstRows[edges]],
     to = probes[dstCols[edges]],
     direction = rep(0, nEdges),
     method = rep("M0039", nEdges),
     weight = if (weighted) adjDst[edges] else rep(1, nEdges)
     )

  if (!is.null(file))
    write.table(visAntData, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  invisible(visAntData)
}



#' Export network to Cytoscape
#' 
#' This function exports a network in edge and node list files in a format
#' suitable for importing to Cytoscape.
#' 
#' If the corresponding file names are supplied, the edge and node data is
#' written to the appropriate files. The edge and node data is also returned as
#' return value (see below).
#' 
#' @param adjMat adjacency matrix giving connection strengths among the nodes
#' in the network.
#' @param edgeFile file name of the file to contain the edge information.
#' @param nodeFile file name of the file to contain the node information.
#' @param weighted logical: should the exported network be weighted?
#' @param threshold adjacency threshold for including edges in the output.
#' @param nodeNames names of the nodes. If not given, \code{dimnames} of
#' \code{adjMat} will be used.
#' @param altNodeNames optional alternate names for the nodes, for example gene
#' names if nodes are labeled by probe IDs.
#' @param nodeAttr optional node attribute, for example module color. Can be a
#' vector or a data frame.
#' @param includeColNames logical: should column names be included in the
#' output files? Note that Cytoscape can read files both with and without
#' column names.
#' @return A list with the following componens: \item{egdeData}{a data frame
#' containing the edge data, with one row per edge} \item{nodeData}{a data
#' frame containing the node data, with one row per node}
#' @author Peter Langfelder
#' @seealso \code{\link{exportNetworkToVisANT}}
#' @keywords misc
#' @export exportNetworkToCytoscape
exportNetworkToCytoscape = function(
  adjMat,
  edgeFile = NULL,
  nodeFile = NULL,
  weighted = TRUE,
  threshold = 0.5,
  nodeNames = NULL,
  altNodeNames = NULL,
  nodeAttr = NULL,
  includeColNames = TRUE)
{
  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0
  nRow = nrow(adjMat)
  checkAdjMat(adjMat, min = -1, max = 1)
  if (is.null(nodeNames)) nodeNames = dimnames(adjMat)[[1]]
  if (is.null(nodeNames))
    stop("Cannot determine node names: nodeNames is NULL and adjMat has no dimnames.")
  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
  colMat = matrix(c(1:nRow), nRow, nRow)

  if (!is.null(nodeAttr))
  {
    if (is.null(dim(nodeAttr))) nodeAttr = data.frame(nodeAttribute = nodeAttr)
    nodeAttr = as.data.frame(nodeAttr)
  } else nodeAttr = data.frame(nodeAttribute = rep(NA, ncol(adjMat)))


  adjDst = as.dist(adjMat)
  dstRows = as.dist(rowMat)
  dstCols = as.dist(colMat)

  edges = abs(adjDst) > threshold
  nEdges = sum(edges)
  edgeData = data.frame (
     fromNode = nodeNames[dstRows[edges]],
     toNode = nodeNames[dstCols[edges]],
     weight = if (weighted) adjDst[edges] else rep(1, nEdges),
     direction = rep("undirected", nEdges),
     fromAltName = if (is.null(altNodeNames)) rep("NA", nEdges) else altNodeNames[dstRows[edges]],
     toAltName = if (is.null(altNodeNames)) rep("NA", nEdges) else altNodeNames[dstCols[edges]]
     )

  nodesPresent = rep(FALSE, ncol(adjMat))
  nodesPresent[dstRows[edges]] = TRUE
  nodesPresent[dstCols[edges]] = TRUE
  nNodes = sum(nodesPresent)
  nodeData = data.frame (
     nodeName = nodeNames[nodesPresent],
     altName = if (is.null(altNodeNames)) rep("NA", nNodes) else altNodeNames[nodesPresent],
     nodeAttr[nodesPresent, ]
     )

  if (!is.null(edgeFile))
    write.table(edgeData, file = edgeFile, quote = FALSE, row.names = FALSE, col.names = includeColNames,
                sep = "\t")

  if (!is.null(nodeFile))
    write.table(nodeData, file = nodeFile, quote = FALSE, row.names = FALSE, col.names = includeColNames,
                sep = "\t")

  list(edgeData = edgeData, nodeData = nodeData)
}
