# Dependency Check --------------------------------------------------------

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required for multivariate normal simulation. Please install it.")
}

# Note: 'stats' is a base package, so requireNamespace isn't strictly needed,
# but we use functions like rnorm, rbinom, and model.matrix from it.

#' Simulate Basic Linear Data
#'
#' @description Generates a foundation dataset with various predictor types 
#' (Normal, Binary, Categorical) and a continuous outcome Y based on a 
#' linear combination of those predictors.
#'
#' @param n Integer. Number of observations.
#' @param betas Numeric vector. The "true" coefficients, including the Intercept.
#' @param predictors Named list. Each element defines a predictor (type, and parameters).
#' @param error_sd Numeric. Standard deviation of the Gaussian noise (residual error).
#'
#' @return A list containing the generated 'data', the 'model_matrix' (X), 
#' underlying 'eta' (linear predictor), and the 'betas' used.
#' @export
simulate_data <- function(n, betas, predictors, error_sd = 1) {
  # Step 1: Initialize dataframe
  df <- data.frame(row_id = seq_len(n))
  
  # Step 2: Generate predictors
  for (pred_name in names(predictors)) {
    spec <- predictors[[pred_name]]
    
    if (spec$type == "normal") {
      df[[pred_name]] <- rnorm(n, mean = spec$mean, sd = spec$sd)
    } else if (spec$type == "binary") {
      df[[pred_name]] <- rbinom(n, 1, spec$prob)
    } else if (spec$type == "categorical") {
      # Determine levels from names or indices
      levels_vec <- if (!is.null(names(spec$probs))) names(spec$probs) else seq_along(spec$probs)
      df[[pred_name]] <- factor(
        sample(levels_vec, n, replace = TRUE, prob = spec$probs),
        levels = levels_vec
      )
    }
  }
  
  # Step 3: Construct model matrix (X)
  # Exclude row_id from the design matrix
  X <- model.matrix(~ . - row_id, data = df)
  
  # Validation: Check if beta length matches the design matrix columns
  if (length(betas) != ncol(X)) {
    stop("Length of betas (", length(betas), ") does not match ncol(X) (", ncol(X), ").\n",
         "Columns: ", paste(colnames(X), collapse = ", "))
  }
  
  # Step 4: Generate linear predictor (eta) and outcome (Y)
  eta <- as.numeric(X %*% betas) 
  y <- eta + rnorm(n, 0, error_sd) 
  df$Y <- y
  
  return(list(
    data = df,
    model_matrix = X,
    betas = betas,
    eta = eta,
    error_sd = error_sd  # Crucial for the ICC derivation in the next step
  ))
}


#' Simulate Clustered Data from Scratch
#'
#' @description Generates a multilevel dataset with random intercepts and 
#' random slopes for a specific predictor across clusters.
#'
#' @param n_clusters Integer. Number of clusters/participants.
#' @param obs_per_cluster Integer. Number of observations per cluster.
#' @param predictors Named list. Specification for predictors.
#' @param betas Numeric vector. Fixed effects coefficients.
#' @param slope_var Character. The name of the predictor that varies randomly.
#' @param error_sd Numeric. Level-1 residual standard deviation.
#' @param icc Numeric. Intraclass Correlation Coefficient (used to derive intercept variance).
#' @param tau1 Numeric. Standard deviation of the random slope.
#' @param cov01 Numeric. Covariance between random intercept and random slope.
#'
#' @export
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
  
  # Step 1: Random effect covariance matrix (D)
  # Derive tau0 (intercept SD) from ICC
  tau0 <- sqrt((icc * error_sd^2) / (1 - icc))
  D <- matrix(c(tau0^2, cov01, cov01, tau1^2), nrow = 2)
  
  # Step 2: Generate random effects per cluster
  u <- MASS::mvrnorm(n_clusters, mu = c(0, 0), Sigma = D)
  colnames(u) <- c("u0", "u1")
  
  data_list <- vector("list", n_clusters)
  
  # Step 3: Loop through clusters
  for (i in 1:n_clusters) {
    cluster_df <- data.frame(
      participant = i,
      lapply(predictors, function(spec) {
        if (spec$type == "normal") {
          rnorm(obs_per_cluster, mean = spec$mean, sd = spec$sd)
        } else if (spec$type == "binary") {
          rbinom(obs_per_cluster, 1, spec$prob)
        } else if (spec$type == "categorical") {
          # Validation for categorical probabilities
          if (is.null(spec$probs)) stop("Categorical predictor must include 'probs'.")
          if (any(spec$probs < 0)) stop("Probabilities must be non-negative.")
          if (abs(sum(spec$probs) - 1) > 1e-6) stop("Probabilities must sum to 1.")
          
          levels_vec <- if (!is.null(names(spec$probs))) {
            names(spec$probs)
          } else {
            paste0("L", seq_along(spec$probs))
          }
          sample(levels_vec, obs_per_cluster, replace = TRUE, prob = spec$probs)
        }
      })
    )
    
    # Construct fixed effects model matrix
    X <- model.matrix(~ . - participant, data = cluster_df)
    
    if (length(betas) != ncol(X)) {
      stop("Length of betas does not match number of predictors (including intercept).")
    }
    
    # Calculate Y with fixed effects + random intercept + random slope
    eta <- as.numeric(X %*% betas)
    
    # Crucial: cast slope_var to numeric to avoid factor math errors
    slope_val <- as.numeric(cluster_df[[slope_var]])
    
    cluster_df$Y <- eta + u[i, "u0"] + u[i, "u1"] * slope_val +
      rnorm(obs_per_cluster, 0, error_sd)
    
    data_list[[i]] <- cluster_df
  }
  
  full_df <- do.call(rbind, data_list)
  
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


#' Add Clustering and Random Slopes to an Existing Simulation
#'
#' @description Takes a simulation object (from simulate_data) and expands it 
#' by adding clusters (participants), random intercepts, and random slopes.
#'
#' @param sim List. The output object from 'simulate_data()'.
#' @param obs_per_participant Integer. Observations per cluster.
#' @param icc Numeric. Intraclass Correlation Coefficient.
#' @param tau1 Numeric. Standard deviation of the random slope.
#' @param cov01 Numeric. Covariance between random intercept and slope.
#' @param slope_var Character. Predictor name that carries the random slope.
#'
#' @export
add_cluster_random_slope <- function(sim,
                                     obs_per_participant = 1,
                                     icc,
                                     tau1,
                                     cov01 = 0,
                                     slope_var = NULL) {
  df <- sim$data
  n_clusters <- nrow(df)
  
  if (!slope_var %in% names(df)) {
    stop(paste0("slope_var '", slope_var, "' not found in data."))
  }
  
  # Derive tau0 from ICC
  tau0 <- sqrt((icc * sim$error_sd^2) / (1 - icc))
  D <- matrix(c(tau0^2, cov01, cov01, tau1^2), nrow = 2)
  
  # Generate random effects
  u <- MASS::mvrnorm(n_clusters, mu = c(0, 0), Sigma = D)
  colnames(u) <- c("u0", "u1")
  
  # Expand data rows for clustering
  df <- df[rep(1:n_clusters, each = obs_per_participant), , drop = FALSE]
  df$participant <- rep(1:n_clusters, each = obs_per_participant)
  
  # Update outcome Y: Original Y + u0 + (u1 * slope_var)
  # Ensure slope_var is numeric
  slope_val <- as.numeric(df[[slope_var]])
  df$Y <- df$Y + u[df$participant, "u0"] + u[df$participant, "u1"] * slope_val
  
  # Update the simulation object
  sim$data <- df
  sim$u0 <- u[, "u0"]
  sim$u1 <- u[, "u1"]
  sim$D <- D
  sim$icc <- icc
  sim$n_clusters <- n_clusters
  sim$obs_per_participant <- obs_per_participant
  
  return(sim)
}



#' Set Outcome Distribution
#'
#' @description Replaces the existing outcome variable with values drawn from 
#' a specified probability distribution.
#'
#' @param data Dataframe. The simulated dataset.
#' @param outcome Character. The name of the column to replace.
#' @param dist Character. "normal", "binary", "poisson", "gamma", or "uniform".
#' @param ... Additional arguments passed to the distribution functions (e.g., lambda, shape).
#'
#' @return The dataframe with the modified outcome variable.
#' @export
set_outcome_distribution <- function(data, outcome = "Y", dist = "normal", ...) {
  n <- nrow(data)
  args <- list(...)
  
  # Helper to handle default values without needing rlang/%||%
  get_arg <- function(name, default) {
    if (name %in% names(args)) return(args[[name]]) else return(default)
  }
  
  data[[outcome]] <- switch(
    dist,
    "normal"  = rnorm(n, mean = get_arg("mean", 0), sd = get_arg("sd", 1)),
    "binary"  = rbinom(n, 1, prob = get_arg("prob", 0.5)),
    "poisson" = rpois(n, lambda = get_arg("lambda", 1)),
    "gamma"   = rgamma(n, shape = get_arg("shape", 2), scale = get_arg("scale", 1)),
    "uniform" = runif(n, min = get_arg("min", 0), max = get_arg("max", 1)),
    stop("Unsupported distribution type.")
  )
  
  message("Replaced outcome with ", dist, " distribution.")
  return(data)
}


#' Add Level-2 Covariate
#'
#' @description Generates a covariate that varies only at the cluster level and 
#' optionally injects its effect into the outcome Y.
#'
#' @param data Dataframe.
#' @param id_var Character. The column identifying clusters (e.g., "participant").
#' @param cov_name Character. Name for the new covariate.
#' @param dist Character. "normal", "uniform", or "binary".
#' @param mean Numeric. Mean for normal/binary, or min for uniform.
#' @param sd Numeric. SD for normal, or max for uniform.
#' @param beta Numeric. The effect size to add to the outcome Y.
#' @param outcome Character. The outcome column to modify.
#'
#' @export
add_level2_covariate <- function(data, id_var, cov_name = "Z",
                                 dist = "normal", mean = 0, sd = 1,
                                 beta = 0, outcome = NULL) {
  
  ids <- unique(data[[id_var]])
  n_ids <- length(ids)
  
  Z <- if (dist == "normal") {
    rnorm(n_ids, mean = mean, sd = sd)
  } else if (dist == "uniform") {
    runif(n_ids, min = mean, max = sd)
  } else if (dist == "binary") {
    rbinom(n_ids, 1, mean)
  } else {
    stop("Unsupported distribution.")
  }
  
  cov_df <- data.frame(id = ids, temp_Z = Z)
  names(cov_df) <- c(id_var, cov_name)
  
  data <- merge(data, cov_df, by = id_var)
  
  if (!is.null(outcome) && beta != 0) {
    data[[outcome]] <- data[[outcome]] + beta * data[[cov_name]]
  }
  
  return(data)
}


#' Add Bounded Variable
#'
#' @description Adds a variable (like time) that is bounded between two values.
#'
#' @param data Dataframe.
#' @param id_var Character. Cluster ID variable.
#' @param var_name Character. Name for the new variable.
#' @param lower Numeric. Lower bound.
#' @param upper Numeric. Upper bound.
#' @param n_points Integer. If provided, overrides cluster size.
#' @param random Logical. If TRUE, samples randomly; if FALSE, uses seq().
#' @param center Logical. If TRUE, centers the variable.
#' @param scale Logical. If TRUE, scales to unit variance.
#' @param outcome Character. Outcome to modify.
#' @param beta Numeric. Effect size.
#'
#' @export
add_bounded_variable <- function(data, id_var, var_name = "time",
                                 lower = 0, upper = 10, n_points = NULL,
                                 random = FALSE, center = FALSE, scale = FALSE,
                                 outcome = NULL, beta = 0) {
  
  # Split-Apply-Combine pattern
  data_list <- lapply(split(data, data[[id_var]]), function(cluster_data) {
    n_obs <- nrow(cluster_data)
    n <- if (!is.null(n_points)) n_points else n_obs
    
    t_vals <- if (!random) {
      seq(lower, upper, length.out = n)
    } else {
      sort(runif(n, min = lower, max = upper))
    }
    
    if (center) t_vals <- t_vals - mean(t_vals)
    if (scale) t_vals <- as.numeric(scale(t_vals, center = FALSE))
    
    cluster_data[[var_name]] <- t_vals
    
    if (!is.null(outcome) && beta != 0) {
      cluster_data[[outcome]] <- cluster_data[[outcome]] + beta * cluster_data[[var_name]]
    }
    return(cluster_data)
  })
  
  return(do.call(rbind, data_list))
}


#' Add Correlation to Predictors
#'
#' @description Uses a multivariate normal approach to impose a specific 
#' correlation structure on existing numeric predictors.
#'
#' @param df Dataframe.
#' @param predictors Vector. Names of columns to correlate.
#' @param cor_matrix Matrix. A valid correlation matrix.
#' @param betas Named vector. If provided, regenerates Y to maintain effect sizes.
#' @param group_var Character. Categorical variable to handle in Y regeneration.
#'
#' @export
add_correlation <- function(df, predictors, cor_matrix, betas = NULL, group_var = "Group") {
  
  if (length(predictors) < 2) {
    warning("add_correlation() needs at least two predictors.")
    return(df)
  }
  
  # Extract numeric data
  X <- as.matrix(df[, predictors, drop = FALSE])
  means <- colMeans(X)
  sds <- apply(X, 2, sd)
  
  # Reconstruct covariance matrix
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
  
  # Generate correlated data
  X_corr <- MASS::mvrnorm(n = nrow(df), mu = means, Sigma = cov_matrix)
  df[, predictors] <- X_corr
  
  # Optional: Regenerate Y because changing X changes the relationship with Y
  if (!is.null(betas)) {
    intercept <- if ("(Intercept)" %in% names(betas)) betas["(Intercept)"] else 0
    eta <- rep(intercept, nrow(df))
    
    for (pred in predictors) {
      if (pred %in% names(betas)) eta <- eta + betas[pred] * df[[pred]]
    }
    
    # Handle categorical groups
    if (group_var %in% names(df)) {
      lvls <- unique(df[[group_var]])
      for (lvl in lvls[-1]) {
        term <- paste0(group_var, lvl)
        if (term %in% names(betas)) {
          eta <- eta + betas[term] * (df[[group_var]] == lvl)
        }
      }
    }
    df$Y <- eta + rnorm(nrow(df), 0, 1)
  }
  
  return(df)
}


#' Add Heteroskedasticity
#'
#' @description Modifies Y such that the residual variance is a function 
#' of a specific predictor (variance increases/decreases with X).
#'
#' @param data Dataframe.
#' @param eta Numeric vector. The linear predictor (X %*% beta).
#' @param predictor Character. The name of the predictor driving the variance.
#' @param sigma0 Numeric. Baseline standard deviation.
#' @param factor Numeric. Magnitude of the variance increase.
#'
#' @export
add_heteroskedasticity <- function(data, eta, predictor, sigma0 = 1, factor = 0.5) {
  # sigma_i = baseline * (1 + factor * X)
  sigma_i <- sigma0 * (1 + factor * data[[predictor]])
  data$Y <- eta + rnorm(nrow(data), mean = 0, sd = sigma_i)
  
  return(data)
}


#' Scale All Numeric Predictors
#' @export
scale_predictors <- function(data, outcome = "Y") {
  numeric_cols <- setdiff(names(data)[sapply(data, is.numeric)], outcome)
  data[numeric_cols] <- lapply(data[numeric_cols], scale)
  return(data)
}

#' Advanced Centering and Scaling
#' @description Supports group-mean and grand-mean centering.
#' @export
center_scale_vars <- function(data, id_var, group_mean = NULL, 
                              grand_mean = NULL, scale = TRUE, suffix = "_cs") {
  
  df <- as.data.frame(data)
  
  # Group-mean centering (Within-cluster)
  if (!is.null(group_mean)) {
    for (var in group_mean) {
      gm_name <- paste0(var, suffix)
      x <- as.numeric(df[[var]])
      g <- df[[id_var]]
      df[[gm_name]] <- ave(x, g, FUN = function(xi) xi - mean(xi, na.rm = TRUE))
      if (scale) df[[gm_name]] <- as.numeric(scale(df[[gm_name]], center = FALSE))
    }
  }
  
  # Grand-mean centering
  if (!is.null(grand_mean)) {
    for (var in grand_mean) {
      gm_name <- paste0(var, suffix)
      df[[gm_name]] <- as.numeric(df[[var]] - mean(df[[var]], na.rm = TRUE))
      if (scale) df[[gm_name]] <- as.numeric(scale(df[[gm_name]], center = FALSE))
    }
  }
  
  return(df)
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
