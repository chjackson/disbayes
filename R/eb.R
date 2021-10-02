## Empirical Bayes method
## Call disbayes to estimate posterior modes of hyperparameters
## Return the modes so that they can subsequently be fixed to fit the desired model 

eb_find_modes <- function(dbcall, dbfn=disbayes, hp_fixed_arg){ 
    dbcall$method <- "opt"
    dbcall$hessian <- FALSE
    dbcall$draws <- 0
    dbcall$hp_fixed <- hp_fixed_arg
    dbcall[[1]] <- NULL
    ## Don't allow extra args to be passed through (e.g. to customise the sampler)
    fargs <- names(dbcall)[names(dbcall) %in% names(formals(dbfn))]
    dbcall <- dbcall[fargs]
    ebres <- do.call(dbfn, dbcall)
    ebres$fit$par
}

expand_hpfixed <- function(hp, hp_fixed){
  hpfnew <- vector(nrow(hp), mode="list")
  names(hpfnew) <- hp$pars
  for (i in seq_along(hp$pars)) {
    hpfnew[[i]] <- if (is.null(hp_fixed[[hp$pars[i]]])) FALSE else hp_fixed[[hp$pars[i]]]
  }
  hpfnew
}

.disbayes_hier_hp <- data.frame(
  pars      = c("scf","sinc","scfmale","sd_int","sd_slope"),
  row.names = c("scf","sinc","scfmale","sd_int","sd_slope"),
  stannames = c("lambda_cf[1]","lambda_inc[1]","lambda_cf_male[1]", "sd_inter[1]","sd_slope[1]"),
  stringsAsFactors=FALSE
)

.disbayes_hp <- data.frame(
  pars      = c("scf","sinc"),
  row.names = c("scf","sinc"),
  stannames = c("lambda_cf[1]","lambda_inc[1]"),
  stringsAsFactors=FALSE
)

eb_disbayes <- function(hplist, hp_fixed, dbcall, dbfn, method, dotargs){
  hp <- hplist
  hp$include <- TRUE
  hp_fixed <- expand_hpfixed(hp, hp_fixed)
  hp$vals <- 1
  for (i in seq_along(hp_fixed)) if (is.numeric(hp_fixed[[i]])) hp$vals[i] <- hp_fixed[[i]]
  normalapprox_wanted <- ((method=="opt") && isTRUE(dotargs$hessian) && 
                            !is.null(dotargs$draws) && dotargs$draws> 1)
  unc_wanted <- (method %in% c("mcmc","vb") || normalapprox_wanted)
  hp$eb <- sapply(hp_fixed, isTRUE) # parameters to do empirical Bayes on
  if (unc_wanted && any(hp$eb)) {
    hp_fixed_arg <- hp_fixed
    hp_fixed_arg[hp$eb] <- FALSE
    modes <- eb_find_modes(as.list(dbcall), dbfn=dbfn, hp_fixed_arg)
    hp$vals[hp$eb] <- modes[hp$stannames[hp$eb]]
  } else modes <- NULL 
  tmp <- sapply(hp_fixed, function(x){is.numeric(x) || isTRUE(x)})
  hp$isfixed <- tmp
  hp
}

