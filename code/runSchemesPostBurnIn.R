#' Continue running breeding schemes post burn-in
#'
#' Allows users to continue a simulation, potentially initiated or 'burned-in' with \code{runBurnInSchemes}.
#' Allows user to optionally change the \code{bsp} and selection criterion.
#' Input \code{simulations} are the output of e.g. \code{runBurnInSchemes}.
#' Runs potentially multiple replications and optionally in parallel.
#'
#' @param simulations tibble, each row is a simulation, 2 columns: SimRep and burnInSim.
#' SimRep is an identifier. burnInSim is a list with 3 named elements:
#' "records", "bsp" and "SP"
#' @param nSimCores Integer, number of cores to optionally execute replicate simulations in parallel
#' @param newBSP optional, so you can specify a different bsp for post-burn in sims.
#' @param nPostBurnInCycles Integer number of cycles to run the \code{selCritPop} and \code{selCritPipe} settings.
#' @param productFunc string, Function to advance the product pipeline by one generation
#' @param popImprovFunc string, Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS. Will speed up, for example, genomic predictions when using selCritGRM. Careful to balance with other forms of parallel processing.
#' @param selCritPop string, overrides the selCrit in \code{bsp} for the post burn-in stage.
#' @param selCritPipe string, overrides the selCrit in \code{bsp} for the post burn-in stage.
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#'
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#'
#' @export
runSchemesPostBurnIn<-function(simulations,
                               newBSP=NULL, # so you can change the scheme entirely after burn-in
                               nPostBurnInCycles,
                               selCritPop="parentSelCritBLUP",
                               selCritPipe="productSelCritBLUP",
                               productFunc="productPipelinePostBurnIn",
                               popImprovFunc="popImprovByParentSel",
                               nSimCores=1,
                               nBLASthreads=NULL){

  require(furrr); plan(multisession, workers = nSimCores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

  simulations<-simulations %>%
    dplyr::mutate(SimOutput=future_map2(SimRep,burnInSim,function(SimRep,burnInSim,...){
      # debug
      # burnInSim<-simulations$burnInSim[[1]]
      if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
      cat("******", SimRep, "\n")

      # This CONTINUES where previous sims left off
      ## no initialize step
      ## Keep burn-in stage sim params "SP"
      SP<-burnInSim$SP
      ## specify a potentially new bsp object
      ## (keep checks stored in burn-in stage's bsp)
      if(!is.null(newBSP)){
        bsp<-newBSP; bsp$checks<-burnInSim$bsp$checks
        bsp[["burnInBSP"]]<-burnInSim$bsp
        # years are indexed starting at year==0,
        ## so for 10 burn-in cycles, max value should be 9, store for later
        bsp[["maxYearBurnInStage"]]<-max(burnInSim$records$stageOutputs$year)
      } else { bsp<-burnInSim$bsp }
      ## 'historical' records from burn-in
      records<-burnInSim$records
      ## override burn-in specified product and population improvement funcs
      bsp[["productPipeline"]] <- get(productFunc)
      bsp[["populationImprovement"]] <- get(popImprovFunc)
      bsp[["selCritPipeAdv"]] <- get(selCritPipe)
      bsp[["selCritPopImprov"]] <- get(selCritPop)

      # Post burn-in cycles
      cat("\n"); cat("Post burn-in cycles"); cat("\n")
      for (cycle in 1:nPostBurnInCycles){
        cat(cycle, " ")
        records <- bsp$productPipeline(records, bsp, SP)
        records <- bsp$populationImprovement(records, bsp, SP)
      }

      # Finalize the stageOutputs
      records <- AlphaSimHlpR:::lastCycStgOut(records, bsp, SP)

      return(list(records=records,
                  bsp=bsp,
                  SP=SP))
    },
    nPostBurnInCycles=nPostBurnInCycles,
    selCritPop=selCritPop,
    selCritPipe=selCritPipe,
    productFunc=productFunc,
    popImprovFunc=popImprovFunc,
    nBLASthreads=nBLASthreads,
    newBSP=newBSP))
  plan(sequential)
  return(simulations)
}
