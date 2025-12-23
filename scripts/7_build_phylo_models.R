# Build phylogenetic logistic regression models

# required library
library(dplyr)
library(phytools)
library(phylolm)
library(parallel) # if using parallelization

# ---------------------------------------------
# Functions
#----------------------------------------------
if (TRUE) {
  # Generate dataset for model fitting
  subset_data <- function(data, nest) {
    list_var <- c("Migration", "Egg", "Ground", "Cooperative", "Clutch", "pred", "PC1", "PC2")
    if (nest == "enclosed") {
      rs <- data %>% mutate(Y = Enclosed>0, N = Open>0 & Enclosed==0) %>% filter(Y+N > 0) %>%
        mutate_at(list_var, ~(scale(.) %>% as.vector))
    } else if (nest == "dome") {
      rs <- data %>% mutate(Y = Dome>0, N = Open>0 & Enclosed==0) %>% filter(Y+N > 0) %>%
        mutate_at(list_var, ~(scale(.) %>% as.vector))
    } else if (nest == "cavity") {
      rs <- data %>% mutate(Y = Cavity>0, N = Open>0 & Enclosed==0) %>% filter(Y+N > 0) %>%
        mutate_at(list_var, ~(scale(.) %>% as.vector))
    }
    return(rs)
  }
  
  # Select initial beta values
  find_best_ini_beta <- function(data, tree, n.random.beta, logfile) {
    # Generate n sets of random initial beta values
    set.seed(123)
    if (model_version == "full") {
      beta <- matrix(runif(n.random.beta * 16, -.2, .2), ncol = 16) # 16 betas per set
    } else {
      beta <- matrix(runif(n.random.beta * 9, -.2, .2), ncol = 9)
    }
    
    # Build models using the initial beta values
    mlist <- list()
    for (i in 1:n.random.beta) {
      cat(Sys.time(), paste0("Fitting model ", i, "\n"), file = logfile, append = TRUE) # print progress in log file
      tryCatch({
        withCallingHandlers(
          {
            mlist[[i]] <- phyloglm(formula, data=data, phy=tree, method="logistic_MPLE", full.matrix=T, btol=50, boot=0, start.beta=beta[i,])
          }, 
          warning = function(w) {
            msg <- paste(Sys.time(), "Warning in model", i, ":", conditionMessage(w), "\n")
            cat(msg, file = logfile, append = TRUE)
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
      }, 
      error = function(e) {
        msg <- paste(Sys.time(), "Error in model", i, ":", conditionMessage(e), "\n")
        cat(msg, file = logfile, append = TRUE)
      })
    }
    
    # Select the set of beta values with the lowest AIC
    aic <- NULL
    for (i in 1:length(mlist)) {
      if (is.null(mlist[[i]])) {
        aic <- c(aic, NA)
      } else {
        aic <- c(aic, mlist[[i]]$aic)
      }
    }
    idx.best <- which(aic == min(aic, na.rm = T))
    beta.best <- beta[idx.best, ]
    cat(Sys.time(), paste0("Best beta is set ", idx.best, "\n"), file = logfile, append = TRUE) # print progress in log file
    return(list(beta.best, idx.best, aic, mlist))
  }
  
  # Fit model using given initial beta values
  fit_model <- function(subset, nest, contree, n.random.beta, nboot, path_output, btol = 20) {
    # print progres to log file
    dir.create(file.path(path_output), showWarnings = FALSE)
    logfile <- paste0(path_output, "log_", model_version, "_", sub("dt_", "", deparse(substitute(subset))), "_", nest, ".txt")
    msg <- paste(Sys.time(), "Starting model.\n")
    cat(msg, file = logfile, append = TRUE)

    # get data subset
    data <- subset_data(subset, nest)
    tree <- contree %>% drop.tip(.$tip.label[!.$tip.label %in% rownames(data)])
    
    # find best fit initial beta values
    beta.ini <- find_best_ini_beta(data, tree, n.random.beta, logfile)
    saveRDS(beta.ini, paste0(path_output, "beta_model_", sub("dt_", "", deparse(substitute(subset))), "_", nest, ".rds")) # optional
    cat(Sys.time(), "Initial beta computed. Running final model.\n", file = logfile, append = TRUE) # print progress in log file

    # fit model using the best initial betas and save results
    set.seed(123)
    m <- phyloglm(formula, data=data, phy=tree, method="logistic_MPLE", start.beta=beta.ini[[1]], btol=btol, boot=nboot, full.matrix=T)
    saveRDS(m, paste0(path_output, "model_", sub("dt_", "", deparse(substitute(subset))), "_", nest, ".rds"))

    cat(Sys.time(), "Finished and saved model.\n", file = logfile, append = TRUE) # print progress in log file
    return(m)
  }
  
  # Fit model using multiple phylogenetic tree (supplementary analysis) with given initial beta values
  fit_model_multi_trees <- function(subset, path_output, nest, trees, beta.ini, cores=NULL) {
    # print progres to log file
    dir.create(file.path(path_output), showWarnings = FALSE)
    logfile <- paste0(path_output, "log_multi_", model_version, "_", deparse(substitute(subset)), "_", nest, ".txt")
    cat(Sys.time(), "Starting model.\n", file = logfile, append = TRUE)

    # btol value
    if (data_version == "" & model_version == "clim") {
      if (deparse(substitute(subset)) == "dt_all") {
        if (nest == "enclosed") btol <- 80
        if (nest == "dome") btol <- 70
      }
      if (deparse(substitute(subset)) == "dt_np") {
        if (nest == "enclosed") btol <- 30
        if (nest == "dome") btol <- 50
        if (nest == "cavity") btol <- 30
      }
    } else {
      btol <- 20
    }

    phyloreg <- function(i) {
      warnings <- character()

      result <- tryCatch({
        withCallingHandlers(
          {
            cat(Sys.time(), paste0("Fitting model using tree ", i, "\n"), file = logfile, append = TRUE) # print progress in log file
            tree <- trees[[i]] %>% drop.tip(.$tip.label[!.$tip.label %in% rownames(data)])
            m <- phyloglm(formula, data=data, phy=tree, method="logistic_MPLE", start.beta=beta.ini, btol=btol, boot=0)
            m
          }, 
          warning = function(w) {
            msg <- paste(Sys.time(), "Warning in tree", i, ":", conditionMessage(w), "\n")
            cat(msg, file = logfile, append = TRUE)
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
      }, 
      error=function(e) {
        msg <- paste(Sys.time(), "Error in tree", i, ":", conditionMessage(e), "\n")
        cat(msg, file = logfile, append = TRUE)
        return(NULL)
      })
      return(list(i = i, model = result, warnings = warnings))
    }
    
    # get data subset
    data <- subset_data(subset, nest)

    # fit models
    if (is.null(cores)) { # no parallelization
      m <- lapply(1:length(trees), phyloreg)
    } else {
      m <- mclapply(1:length(trees), phyloreg, mc.cores = cores) # (parallelized, available on Unix OS)
    }
    saveRDS(m, paste0(path_output, "multitrees_model_", sub("dt_", "", deparse(substitute(subset))), "_", nest, ".rds")) # optional
    cat(Sys.time(), "Finished and saved models.\n", file = logfile, append = TRUE) # print progress in log file
  }
}

# ---------------------------------------------
# Set parameters
#----------------------------------------------
data_version <- ""; model_version <- "pred" # global
data_version <- ""; model_version <- "clim" # global
data_version <- "_north"; model_version <- "pred" # Northern hemisphere
data_version <- "_north"; model_version <- "clim" # Northern hemisphere

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

if (TRUE) {
  # output folder
  path_output <- paste0("output/model", data_version, "_", model_version, "/")

  # Concesnsus tree model parameter
  N_BETA_RAND <- 300 # number of sets of random initial beta values to select from
  N_BOOT <- 1000 # number of bootstrap in the final model fitting

  # Multiple tree model parameter
  ntree <- 1000 # number of trees

  # Create three subsets
  dt_all <- read.csv(paste0("data/processed/trait", data_version, "/sp_traits_pca.csv"), row.names = 1) # all apseices
  dt_np <- dt_all %>% filter(Order != "Passeriformes") # non-passerines
  dt_psr <- dt_all %>% filter(Order == "Passeriformes") # passerines

  # Model formula
  if (model_version == "full") {
    formula <- formula("Y ~ pred + PC1 + PC2 + I(PC2^2) +
                    Ground + Cooperative + Migration + Egg + Clutch +
                    pred*Ground + pred*Cooperative + pred*Clutch +
                    PC1*Egg + PC2*Egg + I(PC2^2)*Egg")
  } else if (model_version == "pred") {
    formula <- formula("Y ~ pred + Ground + Cooperative + Migration + Clutch +
                    pred*Ground + pred*Cooperative + pred*Clutch")
  } else if (model_version == "clim") {
    formula <- formula("Y ~ PC1 + PC2 + I(PC2^2) + Migration + Egg +
                    PC1*Egg + PC2*Egg + I(PC2^2)*Egg")
  }

  ## import phylogenetic trees
  # For main analysis (import consensus tree)
  contree <- readRDS("data/input/tree/consensus_tree_ultrametric.rds")

  # For 1000 tree analysis (import and select candidate phylogenetic trees)
  tree1k <- readRDS("data/input/tree/tree1k.rds") # use Sample1kTree.R to generate file
  set.seed(123)
  trees <- tree1k[sample(1000, ntree)]
  ncores <- 5 # number of cores to run (For no parallelization, set ncores to NULL)
}

#---------------------------------------
# Fit models
#---------------------------------------
#### Fit models and save results for each subset and model type (Main analysis)
# Fit models (parallelized version, available on Unix OS)
# Warning: long runtime (proportional to n.random.beta + nboot)

# tasks <- list(
#   function() fit_model(dt_np,  "enclosed", contree, N_BETA_RAND, N_BOOT, path_output), # 20 min
#   function() fit_model(dt_np,  "dome",     contree, N_BETA_RAND, N_BOOT, path_output), # 20 min
#   function() fit_model(dt_np,  "cavity",   contree, N_BETA_RAND, N_BOOT, path_output), # 25 min
#   function() fit_model(dt_psr, "enclosed", contree, N_BETA_RAND, N_BOOT, path_output), # 25 min
#   function() fit_model(dt_psr, "dome",     contree, N_BETA_RAND, N_BOOT, path_output), # 25 min
#   function() fit_model(dt_psr, "cavity",   contree, N_BETA_RAND, N_BOOT, path_output), # 25 min
#   function() fit_model(dt_all, "enclosed", contree, N_BETA_RAND, N_BOOT, path_output), # 30 min
#   function() fit_model(dt_all, "dome",     contree, N_BETA_RAND, N_BOOT, path_output), # 35 min
#   function() fit_model(dt_all, "cavity",   contree, N_BETA_RAND, N_BOOT, path_output)  # 40 min
# )
# tic <- proc.time()
# results <- mclapply(tasks, function(f) f(), mc.cores = 5) # specify number of cores used (13 hr with 5 cores on iMac) # 4hr
# proc.time() - tic
# results # Check results

tic <- proc.time() # 4.5 hrs (1 core)
fit_model(dt_all, "enclosed", contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_all, "dome",     contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_all, "cavity",   contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_np,  "enclosed", contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_np,  "dome",     contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_np,  "cavity",   contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_psr, "enclosed", contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_psr, "dome",     contree, N_BETA_RAND, N_BOOT, path_output)
fit_model(dt_psr, "cavity",   contree, N_BETA_RAND, N_BOOT, path_output)
proc.time() - tic

#### Fit models using 1,000 candidate phylogenetic trees (Supplementary analysis)
# Fit model and save results
# Warning: long runtime (proportional to ntree)
tic <- proc.time() # 3 hrs (5 cores)
fit_model_multi_trees(dt_all, path_output, "enclosed", trees, beta.ini = readRDS(paste0(path_output, "beta_model_all_enclosed.rds"))[[1]], ncores)
fit_model_multi_trees(dt_all, path_output, "dome",     trees, beta.ini = readRDS(paste0(path_output, "beta_model_all_dome.rds"))[[1]], ncores)
fit_model_multi_trees(dt_all, path_output, "cavity",   trees, beta.ini = readRDS(paste0(path_output, "beta_model_all_cavity.rds"))[[1]], ncores)
fit_model_multi_trees(dt_np,  path_output, "enclosed", trees, beta.ini = readRDS(paste0(path_output, "beta_model_np_enclosed.rds"))[[1]], ncores)
fit_model_multi_trees(dt_np,  path_output, "dome",     trees, beta.ini = readRDS(paste0(path_output, "beta_model_np_dome.rds"))[[1]], ncores)
fit_model_multi_trees(dt_np,  path_output, "cavity",   trees, beta.ini = readRDS(paste0(path_output, "beta_model_np_cavity.rds"))[[1]], ncores)
fit_model_multi_trees(dt_psr, path_output, "enclosed", trees, beta.ini = readRDS(paste0(path_output, "beta_model_psr_enclosed.rds"))[[1]], ncores)
fit_model_multi_trees(dt_psr, path_output, "dome",     trees, beta.ini = readRDS(paste0(path_output, "beta_model_psr_dome.rds"))[[1]], ncores)
fit_model_multi_trees(dt_psr, path_output, "cavity",   trees, beta.ini = readRDS(paste0(path_output, "beta_model_psr_cavity.rds"))[[1]], ncores)
proc.time() - tic
