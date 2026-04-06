
# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(viridis)

# Simulate Data -----------------------------------------------------------

simulate_data <- function(n, betas, predictors, error_sd = 1) {
  # Step 1: Initialize dataframe
  df <- data.frame(row_id = seq_len(n))
  
  # Step 2: Generate predictors
  #Loops through each entry in predictors and is generated based on the spec type listed for that variable
  for (pred_name in names(predictors)) {
    spec <- predictors[[pred_name]]
    
    #Normal: draws random values from a normal distribution
    if (spec$type == "normal") {
      df[[pred_name]] <- rnorm(n, mean = spec$mean, sd = spec$sd)
      #Binary: draws Bernoulli outcomes
    } else if (spec$type == "binary") {
      df[[pred_name]] <- rbinom(n, 1, spec$prob)
      #Randomly samples factor levels using specified probabilities (named or unnamed)
    } else if (spec$type == "categorical") {
      #Determine the levels
      levels_vec <- if (!is.null(names(spec$probs))) names(spec$probs) else seq_along(spec$probs)
      df[[pred_name]] <- factor(
        sample(levels_vec, n, replace = TRUE, prob = spec$probs),
        levels = levels_vec
      )
    }
  }
  
  # Step 3: Constructs the inear predictor and checks coefficients
  
  X <- model.matrix(~ . - row_id, data = df)
  
  if (length(betas) != ncol(X)) {
    stop("Length of betas (", length(betas), ") does not match ncol(X) (", ncol(X), ").\n",
         "Columns: ", paste(colnames(X), collapse = ", "))
  }
  
  #Step 4: Generate linear predictors and outcome
  eta <- as.numeric(X %*% betas) #True model predicted mean without random error
  y <- eta + rnorm(n, 0, error_sd) #Adds random Gaussion noise with SD = error_sd
  df$Y <- y
  
  # Step 5: Return Structured output
  return(list(
    data = df, #full simulated dataset
    model_matrix = X, #The design matrix used
    betas = betas, #The input coefficient vector
    eta = eta, #The true underlying means
    error_sd = error_sd #The residual SD used
  ))
}


# Simulate Clustered Data -------------------------------------------------

simulate_clustered_data <- function(n_clusters,
                                    obs_per_cluster = 1,
                                    predictors,
                                    betas,
                                    slope_var,
                                    error_sd = 1,
                                    icc = 0.2,
                                    tau1 = 0.5,
                                    cov01 = 0) {
  if (!slope_var %in% names(predictors)) {
    stop("slope_var must be one of the predictors.")
  }
  
  # Random effect covariance matrix
  tau0 <- sqrt((icc * error_sd^2) / (1 - icc))
  D <- matrix(c(tau0^2, cov01, cov01, tau1^2), nrow = 2)
  
  # Generate random intercepts/slopes per cluster
  u <- MASS::mvrnorm(n_clusters, mu = c(0, 0), Sigma = D)
  colnames(u) <- c("u0", "u1")
  
  # Initialize list to collect cluster data
  data_list <- vector("list", n_clusters)
  
  for (i in 1:n_clusters) {
    # Generate predictors for obs_per_cluster observations
    cluster_df <- data.frame(
      participant = i,
      lapply(predictors, function(spec) {
        if (spec$type == "normal") {
          rnorm(obs_per_cluster, mean = spec$mean, sd = spec$sd)
        } else if (spec$type == "binary") {
          rbinom(obs_per_cluster, 1, spec$prob)
        } else if (spec$type == "categorical") {
          sample(names(spec$probs), obs_per_cluster, replace = TRUE, prob = spec$probs)
        }
      })
    )
    
    # Construct model matrix for fixed effects
    X <- model.matrix(~ . - participant, data = cluster_df)
    
    if (length(betas) != ncol(X)) {
      stop("Length of betas does not match number of predictors (including intercept).")
    }
    
    # Fixed part
    eta <- as.numeric(X %*% betas)
    
    # Add random effects and residual
    cluster_df$Y <- eta + u[i, "u0"] + u[i, "u1"] * cluster_df[[slope_var]] +
      rnorm(obs_per_cluster, 0, error_sd)
    
    data_list[[i]] <- cluster_df
  }
  
  # Combine all clusters
  full_df <- do.call(rbind, data_list)
  
  # Return structured list
  return(list(
    data = full_df,
    betas = betas,
    u0 = u[, "u0"],
    u1 = u[, "u1"],
    D = D,
    error_sd = error_sd,
    icc = icc,
    n_clusters = n_clusters,
    obs_per_cluster = obs_per_cluster
  ))
}


# Set Y Distribution ------------------------------------------------------

set_outcome_distribution <- function(data, outcome = "Y", dist = "normal", ...) {
  n <- nrow(data)
  args <- list(...)
  
  data[[outcome]] <- switch(
    dist,
    "normal"   = rnorm(n, mean = args$mean %||% 0, sd = args$sd %||% 1),
    "binary"   = rbinom(n, 1, prob = args$prob %||% 0.5),
    "poisson"  = rpois(n, lambda = args$lambda %||% 1),
    "gamma"    = rgamma(n, shape = args$shape %||% 2, scale = args$scale %||% 1),
    "uniform"  = runif(n, min = args$min %||% 0, max = args$max %||% 1),
    stop("Unsupported distribution type.")
  )
  
  message("Replaced outcome with ", dist, " distribution.")
  return(data)
}


# Add Level 2 Covariate ---------------------------------------------------

add_level2_covariate <- function(data, id_var, cov_name = "Z",
                                 dist = "normal", mean = 0, sd = 1,
                                 beta = 0, outcome_var = NULL) {
  
  # 1. Generate one Z per cluster
  ids <- unique(data[[id_var]])
  if (dist == "normal") {
    Z <- rnorm(length(ids), mean = mean, sd = sd)
  } else if (dist == "uniform") {
    cov_values <- runif(n_clusters, min = mean, max = sd)
  } else if (dist == "binary") {
    Z <- rbinom(length(ids), 1, mean)
  }
  
  cov_df <- data.frame(id = ids, Z = Z)
  names(cov_df)[2] <- cov_name
  
  # 2. Merge back into full data
  data <- merge(data, cov_df, by.x = id_var, by.y = "id")
  
  # 3. Optionally inject it into Y (if outcome_var + beta provided)
  if (!is.null(outcome_var) && beta != 0) {
    data[[outcome_var]] <- data[[outcome_var]] + beta * data[[cov_name]]
  }
  
  return(data)
}


# Add Bounded Variable ----------------------------------------------------

add_bounded_variable <- function(data,
                                 id_var,
                                 var_name = "time",
                                 lower = 0,
                                 upper = 10,
                                 n_points = NULL,
                                 random = FALSE,
                                 center = FALSE,
                                 scale = FALSE,
                                 outcome_var = NULL,
                                 beta = 0) {
  # Work on a copy of the data
  df <- data
  
  df <- do.call(rbind, lapply(split(df, df[[id_var]]), function(cluster_data) {
    n_obs <- nrow(cluster_data)
    n <- if (!is.null(n_points)) n_points else n_obs
    
    # Generate time values
    if (!random) {
      t_vals <- seq(lower, upper, length.out = n)
    } else {
      t_vals <- sort(runif(n, min = lower, max = upper))
    }
    
    # Center and/or scale if requested
    if (center) t_vals <- t_vals - mean(t_vals)
    if (scale) t_vals <- as.numeric(scale(t_vals, center = FALSE))
    
    # Add to dataframe
    cluster_data[[var_name]] <- t_vals
    
    # Optionally modify Y based on time effect
    if (!is.null(outcome_var) && beta != 0) {
      cluster_data[[outcome_var]] <- cluster_data[[outcome_var]] + beta * cluster_data[[time_name]]
    }
    
    return(cluster_data)
  }))
  
  return(df)
}


# Adding Correlation to Data ----------------------------------------------

add_correlation <- function(df, predictors, cor_matrix, betas = NULL, group_var = "Group") {
  # df: dataframe from simulate_data
  # predictors: vector of predictor names to correlate (e.g. c('X1', 'X2'))
  # cor_matrix: desired correlation matrix (square, symmetric, valid)
  # betas: optional named vector of coefficients (e.g. c("(Intercept)"=1, "X1"=0.5, "X2"=0.3, "Group2"=0.2, "Group3"=-0.2))
  # group_var: name of categorical variable (default "Group")
  
  if (length(predictors) < 2) {
    warning("add_correlation() needs at least two predictors.")
    return(df)
  }
  
  # --- Step 1: Extract predictor data ---
  X <- df[, predictors, drop = FALSE]
  
  # --- Step 2: Compute means & SDs ---
  means <- colMeans(X)
  sds <- apply(X, 2, sd)
  
  # --- Step 3: Build covariance matrix ---
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
  
  # --- Step 4: Generate correlated data ---
  X_corr <- MASS::mvrnorm(n = nrow(df), mu = means, Sigma = cov_matrix)
  df[, predictors] <- X_corr
  
  # --- Step 5: Recompute Y if betas are supplied ---
  if (!is.null(betas)) {
    # build linear predictor η based on formula-like structure
    intercept <- ifelse("(Intercept)" %in% names(betas), betas["(Intercept)"], 0)
    eta <- rep(intercept, nrow(df))
    
    for (pred in predictors) {
      if (pred %in% names(betas)) {
        eta <- eta + betas[pred] * df[[pred]]
      }
    }
    
    # handle categorical group effects if present
    if (group_var %in% names(df)) {
      levels_group <- unique(df[[group_var]])
      for (lvl in levels_group[-1]) {  # skip reference
        term <- paste0(group_var, lvl)
        if (term %in% names(betas)) {
          eta <- eta + betas[term] * ifelse(df[[group_var]] == lvl, 1, 0)
        }
      }
    }
    
    # finally regenerate Y with fresh residuals
    df$Y <- eta + rnorm(nrow(df), 0, 1)
  }
  
  return(df)
}



# Adding Heteroskedasticity to Data ---------------------------------------

add_heteroskedasticity <- function(data, eta, predictor, sigma0 = 1, factor = 0.5) {
  # Note: This modifier is designed only for continuous outcomes with approximately normal errors.
  # For other models with non-normal outcomes, this model will not work.
  
  
  # df: dataset from simulate_data
  # eta: linear predictor X %*% beta
  # predictor: predictor name that drives the variance
  # sigma0: baseline sd
  # factor: how strongly variance increases with predictor (multiplicative coef)
  sigma_i <- sigma0 * (1 + factor * data[[predictor]])
  y <- eta + rnorm(nrow(data), mean = 0, sd = sigma_i)
  data$Y <- y
  
  data
}

# Adding Clustering and Random Slopes to Data -----------------------------

library(MASS)

add_cluster_random_slope <- function(sim,
                                     obs_per_participant = 1,
                                     icc,
                                     tau1,
                                     cov01 = 0,
                                     slope_var = NULL) {
  # Extract data
  df <- sim$data
  n_clusters <- nrow(df)
  
  # Check slope_var
  if (!slope_var %in% names(df)) {
    stop(paste0("slope_var '", slope_var, "' not found in data."))
  }
  
  # Random effect covariance matrix
  tau0 <- sqrt((icc * sim$error_sd^2) / (1 - icc))  # derive tau0 from ICC
  D <- matrix(c(tau0^2, cov01, cov01, tau1^2), nrow = 2)
  
  # Generate random intercepts and slopes per cluster
  u <- MASS::mvrnorm(n_clusters, mu = c(0, 0), Sigma = D)
  colnames(u) <- c("u0", "u1")
  
  # Expand data by obs_per_participant
  df <- df[rep(1:n_clusters, each = obs_per_participant), , drop = FALSE]
  
  # Assign cluster/participant IDs
  df$participant <- rep(1:n_clusters, each = obs_per_participant)
  
  # Update outcome Y with random effects
  df$Y <- df$Y + u[df$participant, "u0"] + u[df$participant, "u1"] * df[[slope_var]]
  
  # Save updates back into sim
  sim$data <- df
  sim$u0 <- u[, "u0"]
  sim$u1 <- u[, "u1"]
  sim$D <- D
  sim$icc <- icc
  sim$n_clusters <- n_clusters
  sim$obs_per_participant <- obs_per_participant
  
  return(sim)
}


# Scaling Predictors in the Data ------------------------------------------

scale_predictors <- function(data, outcome = "Y") {
  # Automatically detect numeric predictors, excluding outcome
  numeric_cols <- setdiff(names(data)[sapply(data, is.numeric)], outcome)
  
  # Scale each numeric predictor
  data[numeric_cols] <- lapply(data[numeric_cols], scale)
  
  return(data)
}


# Centering and Scaling Predictors ----------------------------------------

center_scale_vars <- function(data, id_var, group_mean = NULL, grand_mean = NULL, scale = TRUE, suffix = "_cs") {
  
  df <- as.data.frame(data)  # ensure base data.frame compatibility
  
  # --- Group-mean centering ---
  if (!is.null(group_mean)) {
    for (var in group_mean) {
      gm_name <- paste0(var, suffix)
      # Ensure both columns are numeric vectors
      x <- as.numeric(df[[var]])
      g <- df[[id_var]]
      
      df[[gm_name]] <- ave(x, g, FUN = function(xi) xi - mean(xi, na.rm = TRUE))
      
      if (scale) {
        df[[gm_name]] <- as.numeric(scale(df[[gm_name]], center = FALSE, scale = TRUE))
      }
    }
  }
  
  # --- Grand-mean centering ---
  if (!is.null(grand_mean)) {
    for (var in grand_mean) {
      gm_name <- paste0(var, suffix)
      df[[gm_name]] <- as.numeric(df[[var]] - mean(df[[var]], na.rm = TRUE))
      
      if (scale) {
        df[[gm_name]] <- as.numeric(scale(df[[gm_name]], center = FALSE, scale = TRUE))
      }
    }
  }
  
  return(df)
}


# Introduce Missing Data --------------------------------------------------

# Main function to add missing data
add_missing_data <- function(data,
                             vars = NULL,
                             cat_vars = NULL,
                             prop_missing,
                             mechanism = c("MCAR", "MAR", "MNAR"),
                             ref_var = NULL) {
  mechanism <- match.arg(mechanism)
  
  # Combine numeric + categorical for validation
  total_vars <- c(vars, cat_vars)
  
  if (length(prop_missing) != length(total_vars)) {
    stop("Length of prop_missing must match number of vars + cat_vars.")
  }
  
  if (mechanism == "MAR" && is.null(ref_var)) {
    stop("MAR mechanism requires specifying a reference variable (ref_var).")
  }
  
  df <- data
  
  # Apply missingness to numeric variables
  if (!is.null(vars)) {
    for (i in seq_along(vars)) {
      var <- vars[i]
      prop <- prop_missing[i]
      
      mask <- switch(
        mechanism,
        "MCAR" = mcar_mask(nrow(df), prop),
        "MAR"  = mar_mask_numeric(df[[ref_var]], prop),
        "MNAR" = mar_mask_numeric(df[[var]], prop)
      )
      
      df[[var]][mask] <- NA
    }
  }
  
  # Apply missingness to categorical variables
  if (!is.null(cat_vars)) {
    for (j in seq_along(cat_vars)) {
      var <- cat_vars[j]
      prop <- prop_missing[length(vars) + j]
      
      mask <- switch(
        mechanism,
        "MCAR" = mcar_mask(nrow(df), prop),
        "MAR"  = mar_mask_categorical(df[[ref_var]], prop),
        "MNAR" = mar_mask_categorical(df[[var]], prop)
      )
      
      tmp <- as.character(df[[var]])
      orig_levels <- levels(df[[var]])
      tmp[mask] <- NA
      df[[var]] <- factor(tmp, levels = orig_levels)
    }
  }
  
  #message("Missing data added via ", mechanism, " mechanism.")
  return(df)
}

#Helper function MCAR
mcar_mask <- function(n, prop) {
  sample(c(TRUE, FALSE), size = n, replace = TRUE, prob = c(prop, 1 - prop))
}

#Helper function MAR/MNAR for numeric variables
mar_mask_numeric <- function(x, prop) {
  p <- plogis(scale(x))  # Logistic transform
  thresh <- quantile(p, probs = 1 - prop, na.rm = TRUE)
  mask <- p > thresh
  return(mask)
}

#Helper function MAR/MNAR for categorical variables
mar_mask_categorical <- function(x, prop) {
  x <- factor(x)               # ensure factor
  levs <- levels(x)
  probs <- runif(length(levs))
  probs <- probs / sum(probs) * prop * length(levs)  # scale so avg = prop
  
  # Create mask
  mask <- vapply(seq_along(x), function(i) {
    val <- x[i]
    if (is.na(val)) return(FALSE)        # keep existing NAs untouched
    idx <- match(val, levs)
    runif(1) < probs[idx]
  }, logical(1))
  
  return(mask)
}

generate_data <- function(subj_n = 30, trial_n = 50, ICC = 0.2) {
  # --- 1. SUBJECT-LEVEL GENERATION (Between-Subjects) ---
  # We define traits that are constant for each person (Anxiety & Experience)
  subject_preds <- list(
    STAIT = list(type = "normal", mean = 50, sd = 15),
    PCL  = list(type = "normal", mean = 30, sd = 10),
    group = list(type = "binary", prob = 0.5) # 0 = Non-Exp, 1 = Exp
  )
  
  # Generate 30 participants
  psych_data <- simulate_data(n = 30, betas = rep(0, 4), predictors = subject_preds)$data %>%
    rename(ID = row_id)
  
  # Define Correlation Matrix for Subject Traits
  # [1] STAIT, [2] PCL, [3] Base_Speed_Factor
  # Note: We set STAIT/PCL to 0.8 to intentionally create Multicollinearity
  cor_mat_subj <- matrix(c(1.0, 0.95, 0.3,
                           0.95, 1.0, 0.2,
                           0.3, 0.2, 1.0), nrow = 3)
  
  # Create a 'Base_Speed_Factor' (Subject-level random intercept/baseline)
  psych_data$Base_Speed_Factor <- rnorm(30, 1, 0.1) 
  
  # Inject the correlations into the subject data
  psych_data <- add_correlation(psych_data, 
                                predictors = c("STAIT", "PCL", "Base_Speed_Factor"), 
                                cor_matrix = cor_mat_subj)
  
  # --- 2. TRIAL-LEVEL GENERATION (Within-Subjects) ---
  # These variables vary across the 50 trials for every participant
  trial_preds <- list(
    distance   = list(type = "normal", mean = 250, sd = sqrt(5000)),
    Weapon_vel = list(type = "normal", mean = 5,   sd = 1),
    Head_vel   = list(type = "normal", mean = 6,   sd = sqrt(1.4)),
    Occlusion  = list(type = "normal", mean = 20,  sd = 10)
  )
  
  # Use clustered simulation to create 1500 rows (30 subjects x 50 trials)
  # 'icc' ensures trials within an ID are related (Violation of Independence)
  trial_sim <- simulate_clustered_data(
    n_clusters = subj_n, 
    obs_per_cluster = trial_n, 
    predictors = trial_preds,
    betas = rep(0.01, 5), # Placeholder coefficients
    slope_var = "distance",
    icc = ICC,           
    error_sd = 0.45 
  )
  
  behavioral_data <- trial_sim$data %>% 
    rename(ID = participant) %>%
    group_by(ID) %>%
    mutate(Trial = row_number()) %>% 
    ungroup()
  
  # --- 3. FINAL ASSEMBLY & DISTRIBUTION MAPPING ---
  # Merge Subject and Trial data, then apply relationships to the outcome (RT)
  data <- behavioral_data %>%
    left_join(psych_data, by = "ID") %>%
    mutate(
      # 1. Map ALL variables to a latent linear predictor (Log Scale)
      # This ensures every variable has a mathematical relationship with RT
      RT_latent = exp(
        -0.7 +                             # Intercept
          (STAIT * 0.004) +                  # Higher anxiety = Slower RT
          (PCL * 0.002) +                   # Trauma symptoms = Slower RT
          (distance * 0.0005) +              # Further distance = Slower RT
          (Weapon_vel * 0.01) +              # Faster weapon = Slower RT (complexity)
          (Head_vel * 0.005) +               # Head movement = Slower RT
          (Occlusion * 0.002) +              # Less visibility = Slower RT
          (ifelse(group == 1, -0.08, 0))     # Experience = Faster RT
      ) * Base_Speed_Factor,               # Individual 'speediness' baseline
      
      # 2. Assign proper distributions to variables
      # Reaction Time: Gamma (Skewed)
      RT = rgamma(n(), shape = 20, scale = RT_latent / 20),
      
      # Shots: Poisson (Counts)
      Shots = rpois(n(), lambda = 2.5 + (group * 1.5)), 
      
      # Target Type: Binary/Factor
      Type = factor(rbinom(n(), 1, 0.5), labels = c("Known", "Unknown")),
      
      # Cleanup Group labels
      group = factor(group, labels = c("Non-Experienced", "Experienced"))
    ) %>%
    # Remove intermediate calculation columns for a clean final object
    dplyr::select(-RT_latent, -Base_Speed_Factor, -Y.x, -Y.y)
  
  return(data)
}