#' roxygen2::roxygenise()

#' the generateCellPopulation function
#'
#' @title generateCellPopulation
#' @param numCells is a vector of length N_cellpop containing the number of cells each should contain
#' @param numProtsPerCell is a N_cellpop x N_prot matrix containing the number of each type of protein appearing in each cellpop (rows are cell populations and columns are proteins).
#' @param plot logical if you should plot the cell information
#' @return \code{igraph} output of this function is the igraph network object for protein network, without any edges

#' @examples
#' numCells = c(25,25,25)
#' numProtsPerCell = rbind(c(5,5,5),c(5,5,5), c(5,5,5))
#' proteinGraph = generateCellPopulation(numCells,numProtsPerCell)
#' plot(proteinGraph)
#'
#' @export

generateCellPopulation = function(numCells, numProtsPerCell, plot = TRUE) {
  # numCells is a vector of length N_cellpop containing the number of cells each should contain
  # numProtsPerCell is a N_cellpop x N_prot matrix containing the number of each type of protein appearing in each cellpop (rows are cell populations and columns are proteins).
  # output of this function is the igraph network object for protein network, without any edges

  # temproary
  # numCells = c(20,20)
  # numProtsPerCell = rbind(c(5,2,4),c(5,5,1))

  require(reshape)
  require(igraph)

  N_cellpop = length(numCells)
  N_prot = ncol(numProtsPerCell)

  if (is.null(names(numCells))) {
    names(numCells) <- paste0("cellPopulation_",1:length(numCells))
  }

  if (is.null(rownames(numProtsPerCell))) {
    rownames(numProtsPerCell) <- names(numCells)
  }

  if (is.null(colnames(numProtsPerCell))) {
    colnames(numProtsPerCell) <- paste0("protein_",1:ncol(numProtsPerCell))
  }

  # build protein nodes information
  m = as.matrix(melt(numProtsPerCell))
  m = m[m[,3] != "0",]
  ms = split(m,1:nrow(m))

  mseach = do.call(rbind,lapply(ms,function(x){
    t(replicate(as.numeric(x[3]), x[1:2], simplify = TRUE))
  }))

  mseachs = split(mseach,1:nrow(mseach))

  proteinNodes = do.call(rbind,lapply(mseachs,function(x){
    n = numCells[x[1]]
    cellnames = paste0(x[1],"_cell_",1:n)
    cbind(cellnames,rep(x[1],n),rep(x[2],n))
  }))
  colnames(proteinNodes) <- c("cellnames", "cellPopulation", "protein")
  proteinNodes <- proteinNodes[order(proteinNodes[,1], proteinNodes[,3]),]
  # rownames(proteinNodes) <- apply(proteinNodes,1,function(x)paste0(x[1],"_",x[3]))

  if (plot) {
    require(ggplot2)

    proteinNodesDF = as.data.frame(proteinNodes)
    g = ggplot(proteinNodesDF,aes(x = protein, fill = cellPopulation)) + geom_bar(position="dodge") +
      theme_minimal() +
      ggtitle("Total number of proteins in simulation")
    print(g)
  }

  # each row in this dataframe should be a node in the network, with the corresponding
  # metadata
  # proteinNodes

  proteinGraph = make_empty_graph(n = nrow(proteinNodes), directed = FALSE)
  V(proteinGraph)$cellnames <- proteinNodes[,1]
  V(proteinGraph)$cellPopulation <- proteinNodes[,2]
  V(proteinGraph)$protein <- proteinNodes[,3]

  return(proteinGraph)
}


#' the proteinGraphToCellGraph function
#'
#' @title proteinGraphToCellGraph
#' @param proteinGraph is an igraph object specifying the protein-protein bindings
#' @param forPlotting optimise the igraph object for plotting (default TRUE, can be set to FALSE for speed)
#' @return \code{igraph} cellGraph object network of cell-cell bindings

#' @examples
#'
#'
#' @export

proteinGraphToCellGraph = function(proteinGraph, forPlotting = TRUE) {
  # this function also includes stylistic things, so that you
  # can just run plot(cellGraph)

  require(igraph)

  numMap = as.numeric(factor(unique(V(proteinGraph)$cellnames,
                                    levels = unique(V(proteinGraph)$cellnames))))
  names(numMap) <- unique(V(proteinGraph)$cellnames)
  numericMap = cbind(V(proteinGraph)$cellnames,
                     numMap[V(proteinGraph)$cellnames]
  )

  cellGraph = simplify(contract(proteinGraph, mapping = numericMap[,2],
                                vertex.attr.comb=function(x)x[1]))

  if (forPlotting) {
    V(cellGraph)$name <- names(numMap)
    V(cellGraph)$label = ""
    # at least 5 populations are hardcoded
    if (length(unique(V(cellGraph)$cellPopulation)) > 5) warning("not enough colours specified")
    pal <- c("#21D921", "#D92121", "#2121D9", "#FFFF4D", "#FF9326")[1:length(unique(V(cellGraph)$cellPopulation))]
    names(pal) <- unique(V(cellGraph)$cellPopulation)
    V(cellGraph)$color = pal[as.numeric(factor(V(cellGraph)$cellPopulation))]
    V(cellGraph)$size = 7

    cellPopulationPairs = t(apply(get.edgelist(cellGraph),1,function(x) sort(V(cellGraph)$cellPopulation[V(cellGraph)$cellnames %in% x])))

    interpolatedColor = apply(cellPopulationPairs,1,function(x){
      colorRampPalette(c(pal[x[1]],pal[x[2]]))(3)[2]
    })

    E(cellGraph)$color = interpolatedColor
    E(cellGraph)$width = 2

    if (FALSE) {
      plot(cellGraph, layout=layout_with_graphopt)
      legend("topright", col = pal[unique(as.numeric(factor(V(cellGraph)$cellPopulation)))], pch = 16,
             legend = unique(V(cellGraph)$cellPopulation), bty = "n")
    }
  }

  return(cellGraph)
}

#' the getproteinStatus function
#'
#' @title getproteinStatus
#' @param proteinGraph is an igraph object specifying the protein-protein bindings
#' @return \code{vector} character vector of "Bound" or "Unbound"
#' @examples
#'
#'
#' @export

getproteinStatus = function(proteinGraph) {
  proteinStatus = ifelse(degree(proteinGraph)>0,"Bound","Unbound")
  # names(proteinStatus = V(proteinGraph)$names)
  return(proteinStatus)
}


#' the add function
#'
#' @title add
#' @param x objects to be added
#' @return added objects
#' @examples
#'
#'
#' @export

add <- function(x) Reduce("+", x)

#' the speedDateCells function
#'
#' @title speedDateCells
#' @param cellGraph is an igraph object specifying the cell-cell bindings
#' @return \code{matrix} a set of cell pairs which should speed date

#' @examples
#'
#'
#' @export

speedDateCells = function(cellGraph) {

  require(igraph)

  # input is cellGraph taken from proteinGraphToCellGraph()
  # output is a pairing of cells
  # g = graph.adjacency(cellNetwork,mode="undirected")
  d_list = replicate(5,{
    positions = layout_with_graphopt(cellGraph)
    rownames(positions) <-  V(cellGraph)$name
    d = as.matrix(dist(positions))
    diag(d) <- NA
    return(d)
  }, simplify = FALSE)
  d = add(d_list)
  # heatmap(d)
  speedDate = matrix("",nrow=floor(vcount(cellGraph)/2),ncol=2)
  for (i in 1:(nrow(speedDate)-1)) {

    # randomly select from rownames(d)
    cell1 = sample(rownames(d),1)

    # rank of distance, rank 1 is closest
    r = rank(d[cell1,!colnames(d) %in% c(cell1)])

    # do 1/distance based probability
    cell2 = sample(names(r),1,prob=(length(r)-r))

    speedDate[i,] <- c(cell1,cell2)

    # remove those cells from the d matrix
    d = d[!rownames(d) %in% c(cell1,cell2),!colnames(d) %in% c(cell1,cell2)]
  }
  # speedDate[nrow(speedDate),] <- rownames(d)
  return(speedDate)
}



#' the bindProteins function
#'
#' @title bindProteins
#' @param proteinGraph is an igraph object specifying the protein-protein bindings
#' @param speedDate is the matrix of cells that should speed date
#' @param bindingAffinity is the matrix of binding affinities for the different proteins
#' @param propSpeedDating proportion of speed dating pairs that should actually bind
#' @param bindingLength number of timesteps in which the highest affinity proteins should bind
#' @return \code{igraph} a proteinGraph object with new protein-protein bindings specified

#' @examples
#'
#'
#' @export

bindProteins = function(proteinGraph,
                        speedDate,
                        bindingAffinity,
                        propSpeedDating = 0.75,
                        bindingLength = 5) {
  # perform the protein pairing with equal prob to the available proteins
  # if proteins are found to be available then they bind for a
  # set length of time depending on the protein and its affinity value
  # do not allow all speedDating cells to bind at once
  # can lead to weird oscillatory behaviour

  require(igraph)

  proteinStatus = getproteinStatus(proteinGraph)

  for (i in 1:ceiling(propSpeedDating*nrow(speedDate))) {
    cells = speedDate[i,]

    # available proteins in first of the cell pair, i.e. proteins in that cell with degree zero
    availableProteins1 = which(V(proteinGraph)$cellnames %in% cells[1] &
                                 degree(proteinGraph) == 0)

    # available proteins in first of the cell pair, i.e. proteins in that cell with degree zero
    availableProteins2 = which(V(proteinGraph)$cellnames %in% cells[2] &
                                 degree(proteinGraph) == 0)

    if (length(availableProteins1) == 0|length(availableProteins2) == 0) next

    protein1 = sample(availableProteins1,1)
    protein2 = sample(availableProteins2,1)

    proteintype1 = V(proteinGraph)$protein[protein1]
    proteintype2 = V(proteinGraph)$protein[protein2]

    pairAffinity = bindingAffinity[proteintype1,proteintype2]

    # number of timesteps where this pair is bound, an integer between 1 and bindingLength
    time = max(1,round(pairAffinity*bindingLength))

    proteinGraph = add.edges(proteinGraph, c(protein1,protein2), attr = list(time = time))

  }

  return(proteinGraph)

}


#' the bindProteins function
#'
#' @title bindProteins
#' @param proteinGraph is an igraph object specifying the protein-protein bindings
#' @return \code{igraph} a proteinGraph object with a countdown on the number of timesteps remaining

#' @examples
#'
#'
#' @export

countDown = function(proteinGraph) {

  require(igraph)

  E(proteinGraph)$time <- E(proteinGraph)$time - 1
  proteinGraph = delete.edges(proteinGraph, which(E(proteinGraph)$time == 0))
  return(proteinGraph)
}


#' the mixingIndex function
#'
#' @title mixingIndex
#' @param proteinGraph is an igraph object specifying the protein-protein bindings (only one of proteinGraph or cellGraph should be specified)
#' @param cellGraph is an igraph object specifying the cell-cell bindings (only one of proteinGraph or cellGraph should be specified)
#' @return \code{list} object containing the t-index and tabulation of number of edges between each pair of proteins

#' @examples
#'
#'
#' @export

mixingIndex = function(proteinGraph = NULL, cellGraph = NULL) {

  require(igraph)

  if (is.null(cellGraph)) cellGraph = proteinGraphToCellGraph(proteinGraph)

  cellPopulationPairs = t(apply(get.edgelist(cellGraph),1,function(x) sort(V(cellGraph)$cellPopulation[V(cellGraph)$cellnames %in% x])))

  # mixing t_index is the proportion of edges that are different between cell populations vs same pairs in cell populations
  t_index = mean(cellPopulationPairs[,1] != cellPopulationPairs[,2])

  levels = unique(c(cellPopulationPairs))

  # tabulation is a cellpopxcellpop table of the number of edges occurring for each of the scenarios
  tabulation = table(factor(cellPopulationPairs[,1],levels = levels),
                     factor(cellPopulationPairs[,2],levels = levels))

  return(list(t_index = t_index,
              tabulation = tabulation))
}


#' the cellAggregator function
#'
#' @title cellAggregator
#' @param numCells vector of the number of cells per cell population
#' @param numProtsPerCell matrix (number of cell populations x number of protein types) specifying number of proteins per cell
#' @param bindingAffinity matrix of binding affinities between proteins (should be symmetric matrix)
#' @param timesteps (default 100) number of timesteps
#' @param burnIn (default 0.75) consider the last burnIn proportion of timesteps for index calculation
#' @param verbose (default TRUE) print messages
#' @param includeListOfGraphs (default TRUE) include list of graphs in results object
#' @param plot (default FALSE) if TRUE prints a breakdown of the input parameters
#' @param ... arguments go to bindProteins()
#' @param cellGraph is an igraph object specifying the cell-cell bindings (only one of proteinGraph or cellGraph should be specified)
#' @return \code{list} cellAggregationResult object containing the input parameters: numCells, numProtsPerCell, bindingAffinity, timesteps, burnIn, verbose, includeListOfGraphs; and output objects: listofCellGraphs (list of cell graphs per timestep), listofProteinGraphs (list of protein graphs per timestep), listoftIndex (list of t index per timestep), listofTabulation (list of tabulation of proteins per timestep), and t_index (final t index summary value across entire simulation).

#' @examples
#'
#' numCells = c(25,25)
#' numProtsPerCell = rbind(c(5,0,0),c(5,5,5))
#' bindingAffinity = cbind(c(1,0,0), c(0,1,0), c(0,0,1))
#'
#' cellAggregationResult = cellAggregator(
#'  numCells,
#'  numProtsPerCell,
#'  bindingAffinity,
#'  timesteps = 100,
#'  burnIn = 0.75,
#'  verbose = TRUE,
#'  includeListOfGraphs = TRUE,
#'  plot = TRUE
#' )
#'
#' cellAggregationResult$t_index
#' unlist(cellAggregationResult$listoftIndex)
#' cellAggregationBarplot(cellAggregationResult)
#'
#' @export


cellAggregator <- function(
  numCells,
  numProtsPerCell,
  bindingAffinity,
  timesteps = 100,
  burnIn = 0.75,
  verbose = TRUE,
  includeListOfGraphs = TRUE,
  plot = FALSE, ...
) {

  require(igraph)

  # ... arguments go to bindProteins

  # burnIn is actually the proportion of the last t_index values that we want to average over
  # so if we choose burnIn = 1 then actually it's averaging over all the timesteps

  listofCellGraphs = list()
  listofProteinGraphs = list()
  listoftIndex = list()
  listofTabulation = list()

  # initialise
  proteinGraph = generateCellPopulation(numCells,numProtsPerCell)

  # ensure bindingAffinity and protein names match
  rownames(bindingAffinity) = colnames(bindingAffinity) = unique(V(proteinGraph)$protein)

  for (timestep in 1:timesteps) {
    if (verbose) print(timestep)

    # speed dating
    cellGraph = proteinGraphToCellGraph(proteinGraph)
    speedDate = speedDateCells(cellGraph)

    # proteins bind
    proteinGraph = bindProteins(proteinGraph, speedDate, bindingAffinity, ...)

    # gather statistics and store objects
    m = mixingIndex(proteinGraph)
    listoftIndex[[timestep]] <- m[["t_index"]]
    listofTabulation[[timestep]] <- m[["tabulation"]]

    if (includeListOfGraphs) {
      listofProteinGraphs[[timestep]] <- proteinGraph
      cellGraph = proteinGraphToCellGraph(proteinGraph)
      listofCellGraphs[[timestep]] <- cellGraph
    }

    if (plot & (timestep %% 10 == 0)) {
      if (!includeListOfGraphs) {
        cellGraph = proteinGraphToCellGraph(proteinGraph)
      }
      plot(cellGraph)
      title(paste0("Timestep = ", timestep))
    }

    if (timestep == timesteps) break
    # countdown ahead of next timestep
    # (don't want to coutndown on the last timestep though)
    proteinGraph <- countDown(proteinGraph)
  }

  # summarise the t_index after a burn-in period
  t_index = mean(unlist(listoftIndex)[ceiling((1-burnIn)*timesteps):timesteps])

  cellAggregationResult = list(
    numCells = numCells,
    numProtsPerCell = numProtsPerCell,
    bindingAffinity = bindingAffinity,
    timesteps = timesteps,
    burnIn = burnIn,
    verbose = verbose,
    includeListOfGraphs = includeListOfGraphs,
    listofCellGraphs = listofCellGraphs,
    listofProteinGraphs = listofProteinGraphs,
    listoftIndex = listoftIndex,
    listofTabulation = listofTabulation,
    t_index = t_index)

  return(cellAggregationResult)
}




#' the cellAggregationBarplot function
#'
#' @title cellAggregationBarplot
#' @param cellAggregationResult result object obtained from running cellAggregator()
#' @return \code{ggplot} ggplot object of barplot of number of edges of each type per timestep

#' @examples
#'
#' numCells = c(25,25)
#' numProtsPerCell = rbind(c(5,0,0),c(5,5,5))
#' bindingAffinity = cbind(c(1,0,0), c(0,1,0), c(0,0,1))
#'
#' cellAggregationResult = cellAggregator(
#'  numCells,
#'  numProtsPerCell,
#'  bindingAffinity,
#'  timesteps = 100,
#'  burnIn = 0.75,
#'  verbose = TRUE,
#'  includeListOfGraphs = TRUE,
#'  plot = TRUE
#' )
#'
#' cellAggregationResult$t_index
#' unlist(cellAggregationResult$listoftIndex)
#' cellAggregationBarplot(cellAggregationResult)
#'
#' @export

cellAggregationBarplot = function(cellAggregationResult) {

  require(igraph)

  cellGraph = cellAggregationResult$listofCellGraphs[[length(cellAggregationResult$listofCellGraphs)]]

  pal <- c("#21D921", "#D92121", "#2121D9", "#FFFF4D", "#FF9326")
  names(pal)[1:length(unique(V(cellGraph)$cellPopulation))] <- unique(V(cellGraph)$cellPopulation)
  # interpolatedColor = apply(cellPopulationPairs,1,function(x){
  #   colorRampPalette(c(pal[x[1]],pal[x[2]]))(3)[2]
  # })

  n = nrow(cellAggregationResult$listofTabulation[[1]])
  tableOfEdges = do.call(rbind,lapply(cellAggregationResult$listofTabulation,melt))
  tableOfEdges = cbind(timestep = rep(1:(nrow(tableOfEdges)/(n*n)), each = n*n),
                       apply(tableOfEdges[,1:2],1,paste0,collapse = "_"),
                       tableOfEdges[,3],
                       tableOfEdges[,1:2])
  colnames(tableOfEdges)[1:3] <- c("timestep", "cellPops", "NumberEdges")
  # tableOfEdges = as.data.frame(tableOfEdges)
  # tableOfEdges$NumberEdges = as.numeric(as.character(tableOfEdges[,3]))
  tableOfEdges$col = apply(tableOfEdges[4:5],1,function(x){
    colorRampPalette(c(pal[x[1]],pal[x[2]]))(3)[2]
  })
  # tableOfEdges$cellPops = factor(tableOfEdges$cellPops, levels = unique(tableOfEdges$cellPops))

  colpal = unique(tableOfEdges[,c(2,6)])
  colpallette = colpal[,2]
  names(colpallette) <- colpal[,1]

  g = ggplot(tableOfEdges, aes(x = timestep, y = NumberEdges, fill = cellPops)) +
    geom_col(position = "stack") + theme_minimal() +
    scale_fill_manual(values = colpallette) +
    theme(legend.position = "bottom")
  g
}
