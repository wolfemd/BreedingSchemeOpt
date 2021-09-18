#' Run burn-in breeding scheme simulations
#'
#' Allows users to run simulation and then continue again later. Output is direct input for \code{runSchemesPostBurnIn}.
#' Runs potentially multiple replications and optionally in parallel.
#'
#' @param nReplications Integer number of replications of the specific breeding scheme to run
#' @param nSimCores Integer, number of cores to optionally execute replicate simulations in parallel
#' @param bsp  A list of breeding scheme parameters.
#' @param nBurnInCycles Integer number of cycles to as 'burn-in' using the \code{selCritPop} and \code{selCritPipe} settings.
#' @param iniFunc string, Function to initialize the breeding program.
#' @param productFunc string, Function to advance the product pipeline by one generation
#' @param popImprovFunc string, Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS. Will speed up, for example, genomic predictions when using selCritGRM. Careful to balance with other forms of parallel processing.
#' @param nThreadsMacs2 uses the nThreads argument in \code{runMacs2}, parallelizes founder sim by chrom.
#' @param selCritPop string, overrides the selCrit in \code{bsp} for the burn-in stage.
#' @param selCritPipe string, overrides the selCrit in \code{bsp} for the burn-in stage.
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#'
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#'
#' @export
runBurnInSchemes<-function(bsp,
                           nBurnInCycles,
                           selCritPop="selCritIID",
                           selCritPipe="selCritIID",
                           iniFunc="initializeScheme",
                           productFunc="productPipeline",
                           popImprovFunc="popImprov1Cyc",
                           nReplications=1,nSimCores=1,
                           nBLASthreads=NULL,nThreadsMacs2=NULL){

  require(furrr); plan(multisession, workers = nSimCores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

  simulations<-tibble(SimRep=1:nReplications) %>%
    mutate(burnInSim=future_map(SimRep,function(SimRep,...){
      if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
      cat("******", SimRep, "\n")

      # This initiates the founding population
      bsp[["initializeFunc"]] <- get(iniFunc)
      bsp[["productPipeline"]] <- get(productFunc)
      bsp[["populationImprovement"]] <- get(popImprovFunc)

      initList <- bsp$initializeFunc(bsp,nThreadsForMacs=nThreadsMacs2)
      SP <- initList$SP
      bsp <- initList$bsp
      records <- initList$records

      ## set the selection criteria for burn-in
      bsp[["selCritPipeAdv"]] <- get(selCritPipe)
      bsp[["selCritPopImprov"]] <- get(selCritPop)

      # Burn-in cycles
      cat("\n"); cat("Burn-in cycles"); cat("\n")
      for (cycle in 1:nBurnInCycles){
        cat(cycle, " ")
        records <- bsp$productPipeline(records, bsp, SP)
        records <- bsp$populationImprovement(records, bsp, SP)
      }
      return(list(records=records,
                  bsp=bsp,
                  SP=SP))
    },
    bsp=bsp,
    nBurnInCycles=nBurnInCycles,
    selCritPop=selCritPop,
    selCritPipe=selCritPipe,
    iniFunc=iniFunc,
    productFunc=productFunc,
    popImprovFunc=popImprovFunc,
    nBLASthreads=nBLASthreads,
    nThreadsMacs2=nThreadsMacs2))
  plan(sequential)
  return(simulations)
}
