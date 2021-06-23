#' Simulate dose-response data
#' 
#' Simulate dose-response single-cell or single-nuclei RNA-seq count data
#' for a dose-response series of samples.
#' 
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{splatParams}} from the Splatter package for details.
#' @param method which simulation method to use. Not currently in use
#' @param sparsify logical. Whether to convert assay data to sparse matrices
#'        in order to reduce the final object size/
#' @param verbose logical. Whether to print progress messages.
#' @param dose.names vector of the doses to simulate in ascending numerical order.
#' @param dose.prob vector of cell proportions for each dose. Must add up to 1 and 
#'        be of the same length as dose.names.
#' 
#' @details 
#' 
#' \enumerate{
#'     \item step1
#'     \item step2
#'     \item step3
#' }
#' 
#' The final output...
#' 
#' @return SingleCellExperiment object containing simulated counts and log-normalized
#'         values.
#' 
#' @references
#' Nault R., ....
#' 
#' @examples 
#' 
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @import splatter
#' @export
splatSimulateDR = function(params = newSplatParams(), 
                           method = c('doseresp'), 
                           sparsify = TRUE, 
                           verbose = FALSE,
                           dose.names = c(0), 
                           dose.prob = c(1),
                           lib.scale = 1.3, ...){
  checkmate::assertClass(params, "SplatParams")
  seed = getParam(params, "seed")
  set.seed(seed)

  # Get the parameters we are going to use
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nBatches <- getParam(params, "nBatches")
  batch.cells <- getParam(params, "batchCells")
  
  # Set up name vectors
  if (verbose) {message("Creating simulation object...")}
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  batch.names <- paste0("Batch", seq_len(nBatches))
  
  # Create SingleCellExperiment to store simulation
  cells =  data.frame(Cell = cell.names)
  rownames(cells) = cell.names
  features = data.frame(Gene = gene.names)
  rownames(features) = gene.names
  sim = SingleCellExperiment(rowData = features, colData = cells,
                              metadata = list(Params = params))
  
  # Create the dose groups
  doses = sample(dose.names, nCells, prob = dose.prob, replace = TRUE)
  colData(sim)$Dose = factor(doses, levels = dose.names)
  
  # Create the batches which will serve as replicates
  if (nBatches > 1) {
    replicates.groups = split(c(1:nBatches), c(1:length(dose.names)))
    replicates.values = sapply(as.numeric(colData(sim)$Dose), 
                               function(x) sample(replicates.groups[[x]], 1))
    colData(sim)$Batch <- as.factor(paste0('Batch',replicates.values))
  }
  
  if (verbose) {message("Simulating library sizes...")}
  sim = splatter:::splatSimLibSizes(sim, params)
  sim$ExpLibSize = sim$ExpLibSize * lib.scale

  if (verbose) {message("Simulating Gene means...")}
  sim = splatter:::splatSimGeneMeans(sim, params)
  
  if (nBatches > 1) {
    if (verbose) {message("Simulating batch effects...")}
    sim <- splatter:::splatSimBatchEffects(sim, params)
  }
  
  if (verbose) {message("Simulating cell means...")}
  sim = splatter:::splatSimBatchCellMeans(sim, params)
  

  if (verbose) {message("Simulating dose-response models...")}
  sim = splatSimDoseResponse(sim, params)
  
  if (verbose) {message("Simulating dose-response models, part 2...")}
  sim = splatSimDoseResponseModel(sim, params, verbose = verbose)

  if (verbose) {message("Simulating BCV...")}
  sim <- splatter:::splatSimBCVMeans(sim, params)
  
  if (verbose) {message("Simulating counts...")}
  sim <- splatter:::splatSimTrueCounts(sim, params)
  
  if (verbose) {message("Simulating dropout (if needed)...")}
  sim <- splatter:::splatSimDropout(sim, params)
  assays(sim)$counts[is.na(assays(sim)$counts)] = 0
  assays(sim)$logcounts = log1p(t(t(assays(sim)$counts)/colSums(assays(sim)$counts))*10000)
  
  return(sim)
}

#' Simulate maximal |fold-change|
#'
#' Simulate the gene-wise scaling factors for differential expression. Scaling 
#' factor is determined following a log-normal distribution using 
#' \code{\link{getLNormFactors}}. These values are use to determine peak 
#' absolute fold change.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with scaling factors for differential expression.
#'
#' @name splatSimDE
#' 
#' @importFrom SummarizedExperiment rowData
#' @export
splatSimDoseResponse <- function(sim, params) {
  
  nGenes <- getParam(params, "nGenes")
  nGroups <- getParam(params, "nGroups")
  de.prob <- getParam(params, "de.prob")
  de.downProb <- getParam(params, "de.downProb")
  de.facLoc <- getParam(params, "de.facLoc")
  de.facScale <- getParam(params, "de.facScale")
  means.gene <- rowData(sim)$GeneMean
  
  
  de.facs <- splatter:::getLNormFactors(nGenes, de.prob, de.downProb, de.facLoc, de.facScale)
  rowData(sim)$DE_idx = de.facs
  
  return(sim)
}

#EDIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' Simulate group differential expression
#'
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taken to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated differential expression.
#'
#' @name splatSimDE
#' 
#' @importFrom SummarizedExperiment rowData
#' @export
splatSimDoseResponseModel = function(sim, params, models.prob = rep(1/6, 6), verbose = FALSE){
  # TODO: Document this block and the following 2
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nDoses = length(unique(colData(sim)$Dose))
  doses <- colData(sim)$Dose
  dose.names <- c(0,0.01,0.03,0.1,0.3,1,3,10,30)
  exp.lib.sizes <- colData(sim)$ExpLibSize
  batch.means.cell <- assays(sim)$BatchCellMeans
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  
  de.ans = rowData(sim)$DE_idx
  de.idx = which(de.ans != 1)
  models.num = floor(models.prob*rep(length(de.idx)))
  while (sum(models.num) < length(de.idx)){
    rand.idx = sample(1:length(models.num), 1)
    models.num[rand.idx] = models.num[rand.idx] + 1
  }
  models = rep(c('Hill','Power','Linear', 'Exp', 'Exp2', 'ExpB'), models.num)
  
  rowData(sim)$Model = 'Unchanged'
  rowData(sim)[sample(de.idx, sum(models.num), replace = FALSE),'Model'] = models
  
  #Create mat
  dose.means = data.frame(matrix(ncol = length(dose.names), nrow = nGenes))
  colnames(dose.means) = paste0("DEFac",dose.names)
  
  m.list = list()
  #Change to apply
  for (ix in 1:nrow(rowData(sim))){
    if (ix%%100 == 0){
      it.percent = (ix/nrow(rowData(sim)))*100
      message(paste(it.percent, '%----'), appendLF = TRUE)
    }
    ix.name = as.character(rowData(sim)[ix,'Gene'])
    if (rowData(sim)[ix, 'Model'] == 'Hill'){
      o = splatsimHill(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], verbose = verbose)
      dose.means[ix,] = o$resp
      o$resp = split(o$resp, dose.names)
      o$gene = ix.name
      m.list[[ix.name]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'Exp'){
      o = splatsimExp(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], verbose = verbose)
      dose.means[ix,] = o$resp
      o$resp = split(o$resp, dose.names)
      o$gene = ix.name
      m.list[[ix.name]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'Exp2'){
      o = splatsimExp(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], power = TRUE, verbose = verbose)
      dose.means[ix,] = o$resp
      o$resp = split(o$resp, dose.names)
      o$gene = ix.name
      m.list[[ix.name]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'ExpB'){
      o = splatsimExpB(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], verbose = verbose)
      dose.means[ix,] = o$resp
      o$resp = split(o$resp, dose.names)
      o$gene = ix.name
      m.list[[ix.name]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'Power'){
      o = splatsimPower(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], verbose = verbose)
      if (is.null(o$resp)){
        print(paste('power_',ix))
        dose.means[ix,] = rep(rowData(sim)[ix,'GeneMean'], 9)
      } else {
        dose.means[ix,] = o$resp
        o$resp = split(o$resp, dose.names)
        o$gene = ix.name
        m.list[[ix.name]] = o
      }
    } else if (rowData(sim)[ix, 'Model'] == 'Linear'){
      o = splatsimLinear(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], verbose = verbose)
      dose.means[ix,] = o$resp
      o$resp = split(o$resp, dose.names)
      o$gene = ix.name
      m.list[[ix.name]] = o
    } else {
      dose.means[ix,] = rep(rowData(sim)[ix,'GeneMean'], 9)
    }
  }
  model.params = do.call(dplyr::bind_rows, lapply(m.list, as.data.frame))
  rownames(model.params) = model.params$gene
  metadata(sim)$modelFits = model.params
  
  # TODO: Document this and the next blocks
  rowData(sim) = cbind(rowData(sim), dose.means)
  dose.facs.gene = dose.means/rowData(sim)$GeneMean
  dose.facs.gene[which(is.na(dose.facs.gene)),'DEFac0'] = 1.0000
  cell.facs.gene = as.matrix(dose.facs.gene[, paste0('DEFac', doses)])
  cell.means.gene <- batch.means.cell * cell.facs.gene
  cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
  
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  assays(sim)$BaseCellMeans <- base.means.cell
  
  colnames(cell.means.gene) <- cell.names
  rownames(cell.means.gene) <- gene.names
  assays(sim)$ScaledCellMeans <- cell.means.gene
  
  return(sim)
}

#EDIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' Simulate group differential expression
#'
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taken to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated differential expression.
#'
#' @name splatSimDE
#' 
#' @importFrom SummarizedExperiment rowData
#' @export
calcFC = function(sim){
  dose_vec = sort(unique(colData(sim)$Dose))
  m0 = rowMeans(as.matrix(logcounts(sim)[,which(colData(sim)$Dose == 0)]))
  fc = list()
  for (dose in dose_vec[-1]){
    temp.means = rowMeans(as.matrix(logcounts(sim)[,which(colData(sim)$Dose == dose)]))
    fc[[dose]] = temp.means-m0
  }
  if (length(dose_vec) > 2){
  fc.out = do.call(cbind, lapply(fc, as.data.frame))
  colnames(fc.out) = paste0('calculatedFC',names(fc))
  } else {
    fc.out <- data.frame(fc[[dose]])
  }
  return(fc.out)
}

#EDIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' Simulate group differential expression
#'
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taken to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated differential expression.
#'
#' @name splatSimDE
#' 
#' @importFrom SummarizedExperiment rowData
#' @export
calcZeroP = function(sim){
  dose_vec = sort(unique(colData(sim)$Dose))
  pz = list()
  for (dose in dose_vec){
    temp.df = as.matrix(logcounts(sim)[,which(colData(sim)$Dose == dose)])
    pz[[dose]] = apply(temp.df, 1, function(x) sum(x == 0)/length(x))
  }
  pz.out = do.call(cbind, lapply(pz, as.data.frame))
  colnames(pz.out) = paste0('percent.zero',names(pz))
  return(pz.out)
}


#' Simulate batch effects
#'
#' Simulate batch effects. Batch effect factors for each batch are produced
#' using \code{\link{getLNormFactors}} and these are added along with updated
#' means for each batch.
#'
#' @param sim SingleCellExperiment to add batch effects to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated batch effects.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
splattSimReplicates <- function(sim, params) {
  
  nGenes <- getParam(params, "nGenes")
  nBatches <- getParam(params, "nBatches")
  batch.facLoc <- getParam(params, "batch.facLoc")
  batch.facScale <- getParam(params, "batch.facScale")
  means.gene <- rowData(sim)$GeneMean
  
  for (idx in seq_len(nBatches)) {
    batch.facs <- splatter:::getLNormFactors(nGenes, 1, 0.5, batch.facLoc[idx],
                                  batch.facScale[idx])
    
    rowData(sim)[[paste0("BatchFacBatch", idx)]] <- batch.facs
  }
  
  return(sim)
}
