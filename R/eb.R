## Empirical Bayes method
## Call disbayes to estimate posterior modes of hyperparameters
## Return the modes so that they can subsequently be fixed to fit the desired model 
## args: names of the disbayes arguments indicating hyperparameters

## TODO can we abstract this, it's used for the RE SD in hier.  then can do it in hier for smoothness too
## arguments: call list, db_argnames, stan_parnames, stan_fixedparnames
## TODO problems if outer call uses mcmc, and we try to recurse disbayes(method="opt", chains=...)
## do we need to pass opt args to mode-finding function?  don't think so. if so could just
## remove all args that are not formals of disbayes. 

## If this recursion causes problems, could just write out all the arguments for safety

eb_find_modes <- function(dbcall, dbfn=disbayes, args){ 
    dbcall$method <- "opt"
    dbcall$hessian <- FALSE
    dbcall$draws <- 0
    for (arg_ind in args)
        dbcall[[arg_ind]] <- NULL 
    dbcall[[1]] <- NULL
    ebres <- do.call(dbfn, dbcall)
    ebres$fit$par
}
