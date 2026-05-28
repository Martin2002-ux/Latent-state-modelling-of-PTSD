
run = function(phase_def,
               n = 1,
               seed = 42,
               my_world = NULL,
               
               alpha0 = "0.05",
               alpha1 = "0.05",
               alpha2 = "0.05",
               gamma  = "0.05",
               eta    = "0.2",
               delta  = "0.6",
               sigma0 = "0.5",
               tau = "10",
               
               chi    = "1",
               ITI = NULL,
               ezITI = NULL,
               contextShift = NULL,
               
               center = 0.5,
               
               initialV = NULL,
               initialex = NULL,
               lsb0 = NULL,
               
               ncop   = 10,
               par_L = 1,
               pseudo = 4,
               block = NULL,
               probabilistic = FALSE,
               features = NULL,
               
               group = "ctrl",
               prev = 0,
               
               plot=FALSE,
               plot_all=FALSE,
               plot_label=FALSE,
               plot_arm=FALSE,
               plot_mu=FALSE,
               plot_shown=FALSE,
               plot_cue=FALSE,
               plot_sd=FALSE,
               plot_raw=FALSE,
               plot_raw_all=FALSE,
               plot_q=FALSE,
               fixed_axis=TRUE
) {
  
  all_results <- vector("list", n)
  set.seed(seed)
  
  # ============================================================
  # Parameter sampling helpers
  # ============================================================
  
  sample_par <- function(x, name) {
    if (!is.character(x)) stop(
      name, " must be passed as a string. ",
      "Use e.g. ", name, " = \"0.05\" for a fixed value or ",
      name, " = \"rnorm(1, 0.05, 0.01)\" for a distribution."
    )
    val <- eval(parse(text = x))[1]
    if (!is.numeric(val) || !is.finite(val)) stop(
      name, " must evaluate to a single finite numeric value, but got: ", deparse(val), "."
    )
    val
  }
  
  for (i in 1:n) {
    
    # ============================================================
    # 1. Create participant parameters for this participant
    # ============================================================
    tau_ind    <- sample_par(tau,    "tau")
    alpha0_ind <- sample_par(alpha0, "alpha0")
    alpha1_ind <- sample_par(alpha1, "alpha1")
    alpha2_ind <- sample_par(alpha2, "alpha2")
    gamma_ind  <- sample_par(gamma,  "gamma")
    eta_ind    <- sample_par(eta,    "eta")
    delta_ind  <- sample_par(delta,  "delta")
    sigma0_ind <- sample_par(sigma0, "sigma0")
    chi_ind    <- sample_par(chi, "chi")
    
    # Parameter range checks
    if (tau_ind <= 0)              stop("tau must be positive, got: ", tau_ind, ".")
    if (alpha0_ind < 0 || alpha0_ind > 1) stop("alpha0 must be in [0, 1], got: ", alpha0_ind, ".")
    if (alpha1_ind < 0 || alpha1_ind > 1) stop("alpha1 must be in [0, 1], got: ", alpha1_ind, ".")
    if (alpha2_ind < 0 || alpha2_ind > 1) stop("alpha2 must be in [0, 1], got: ", alpha2_ind, ".")
    if (gamma_ind  < 0 || gamma_ind  > 1) stop("gamma must be in [0, 1], got: ", gamma_ind,  ".")
    if (eta_ind    < 0)            stop("eta must be non-negative, got: ", eta_ind, ".")
    if (delta_ind  < 0)            stop("delta must be non-negative, got: ", delta_ind, ".")
    if (sigma0_ind <= 0)           stop("sigma0 must be positive, got: ", sigma0_ind, ".")
    if (chi_ind    < 1)            stop("chi must be at least 1, got: ", chi_ind, ".")
    
    # ============================================================
    # 2. Initialize task
    # ============================================================
    
    out   <- TaskInit(phase_def = phase_def,
                      pseudo = pseudo,
                      par_L = par_L,
                      ncop = ncop,
                      center = center,
                      contextShift = contextShift,
                      block = block,
                      probabilistic = probabilistic,
                      features = features)
    world <- out$world
    par   <- out$par
    all_cue_types <- out$cue_names
    
    # ============================================================
    # 3. Set model parameters
    # ============================================================
    
    par$tau <- tau_ind
    
    par$ls <- list(
      alpha0 = alpha0_ind,
      alpha1 = alpha1_ind,
      alpha2 = alpha2_ind,
      gamma  = gamma_ind,
      ncop   = ncop,
      eta    = eta_ind,
      delta  = delta_ind,
      chi    = chi_ind,
      sigma0 = sigma0_ind
    )
    
    # Make sure ITI is in the correct form
    ntrials <- sum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
    

    # --- 1. Disallow using both ITI and ezITI ---
    if (!is.null(ITI) && !is.null(ezITI)) {
      stop("Cannot use both ITI and ezITI at the same time. Please specify only one.")
    }
    
    # --- 2. Validate ITI if provided directly ---
    if (!is.null(ITI)) {
      if (!is.list(ITI) || length(ITI) != n) stop(paste0(
        "ITI must be a list with exactly n (", n, ") vectors, ",
        "but has length ", length(ITI), "."
      ))
      for (j in seq_along(ITI)) {
        if (length(ITI[[j]]) != ntrials) stop(paste0(
          "ITI[[", j, "]] must have length equal to the number of trials (", ntrials, "), ",
          "but has length ", length(ITI[[j]]), "."
        ))
      }
    }
    
    # --- 3. Validate ezITI if used ---
    if (!is.null(ezITI)) {
      if (!is.list(ezITI) || length(ezITI) == 0) {
        stop("ezITI must be a non-empty list of numeric vectors c(trial, ITI).")
      }
      for (j in seq_along(ezITI)) {
        v <- ezITI[[j]]
        if (!is.numeric(v) || length(v) != 2) {
          stop(paste0(
            "Each ezITI[[", j, "]] must be a numeric vector of length 2: c(trial, ITI)."
          ))
        }
        trial_raw <- v[1]
        iti_gap   <- v[2]
        
        if (!is.numeric(trial_raw) || length(trial_raw) != 1 ||
            trial_raw != round(trial_raw)) {
          stop("In ezITI[[", j, "]], trial index ", trial_raw,
               " is not an integer. Use whole trial numbers, e.g. 56, not 56.5.")
        }
        trial_idx <- as.integer(trial_raw)
        if (trial_idx < 1 || trial_idx >= ntrials) {
          stop("In ezITI[[", j, "]], trial index ", trial_idx,
               " must be between 1 and ", ntrials - 1, ".")
        }
        if (!is.finite(iti_gap) || iti_gap <= 1) {
          stop("In ezITI[[", j, "]], the ITI gap ", iti_gap,
               " must be greater than 1 (the default gap is 1).")
        }
      }
    }
    
    # --- 4. Construct par$ITI for this participant ---
    if (!is.null(ITI)) {
      # Use direct per-participant ITI
      par$ITI <- ITI[[i]]
      
    } else if (!is.null(ezITI)) {
      # Build from ezITI (global schedule, same for all participants)
      if (ntrials == 1) {
        par$ITI <- 1
      } else {
        gaps <- rep(1, ntrials - 1)  # default gaps
        
        for (spec in ezITI) {
          trial_idx <- as.integer(spec[1])
          iti_gap   <- spec[2]
          gaps[trial_idx] <- iti_gap
        }
        
        par$ITI <- numeric(ntrials)
        par$ITI[1] <- 1
        for (t in 2:ntrials) {
          par$ITI[t] <- par$ITI[t - 1] + gaps[t - 1]
        }
      }
      
    } else {
      par$ITI <- NULL
    }
    
    par$initialex <- initialex
    par$initialV  <- initialV
    par$lsb0      <- lsb0
    
    L <- par$L * par$ls$ncop
    
    # Cannot use both ITI and ezITI
    if (!is.null(ITI) && !is.null(ezITI)) {
      stop("Cannot use both ITI and ezITI at the same time. Please specify only one.")
    }
    
    # Validate ITI (if provided directly)
    if (!is.null(ITI)) {
      if (length(ITI) != n) stop(paste0(
        "ITI must be a list with exactly n (", n, ") vectors, ",
        "but has length ", length(ITI), "."
      ))
      for (j in seq_along(ITI)) {
        if (length(ITI[[j]]) != ntrials) stop(paste0(
          "ITI[[", j, "]] must have length equal to the number of trials (", ntrials, "), ",
          "but has length ", length(ITI[[j]]), "."
        ))
      }
    }
    
    # ezITI: one global schedule, same for all participants.
    # It must be a non-empty list; each element is c(trial, ITI_gap).
    if (!is.null(ezITI)) {
      if (!is.list(ezITI) || length(ezITI) == 0) {
        stop("ezITI must be a non-empty list of numeric vectors c(trial, ITI).")
      }
      for (j in seq_along(ezITI)) {
        v <- ezITI[[j]]
        if (!is.numeric(v) || length(v) != 2) {
          stop(paste0(
            "Each ezITI[[", j, "]] must be a numeric vector of length 2: c(trial, ITI)."
          ))
        }
      }
    }
    
    #replace generated world with a custom one inserted by the user
    if (!is.null(my_world)) {
      
      # Must be a list
      if (!is.list(my_world)) stop(
        "my_world must be a list with elements $c_vec and $win."
      )
      
      # Must contain both required elements
      if (is.null(my_world$c_vec)) stop(
        "my_world must contain element $c_vec."
      )
      if (is.null(my_world$win)) stop(
        "my_world must contain element $win."
      )
      
      # c_vec must be a 5-dimensional array of the correct size
      if (!is.array(my_world$c_vec) || length(dim(my_world$c_vec)) != 5) stop(
        "my_world$c_vec must be a 5-dimensional array (D x C x L x A x ntrials)."
      )
      if (!all(dim(my_world$c_vec) == c(par$D, par$C, L, par$A, par$ntrials))) stop(paste0(
        "my_world$c_vec has wrong dimensions. ",
        "Expected (", paste(c(par$D, par$C, L, par$A, par$ntrials), collapse = " x "), "), ",
        "got (", paste(dim(my_world$c_vec), collapse = " x "), ")."
      ))
      
      # win must be a 3-dimensional array of the correct size
      if (!is.array(my_world$win) || length(dim(my_world$win)) != 3) stop(
        "my_world$win must be a 3-dimensional array (D x A x ntrials)."
      )
      if (!all(dim(my_world$win) == c(par$D, par$A, par$ntrials))) stop(paste0(
        "my_world$win has wrong dimensions. ",
        "Expected (", paste(c(par$D, par$A, par$ntrials), collapse = " x "), "), ",
        "got (", paste(dim(my_world$win), collapse = " x "), ")."
      ))
      
      world <- my_world
      # Rebuild type_seq with per-phase structure matching phase_def
      flat_seq <- unlist(my_world$type_seq)
      phase_lengths <- sapply(phase_def, function(ph) sum(unlist(ph$trials)))
      phase_ends_idx <- cumsum(phase_lengths)
      phase_starts_idx <- c(1, phase_ends_idx[-length(phase_ends_idx)] + 1)
      out$type_seq <- setNames(lapply(seq_along(phase_def), function(ph) {
        as.list(flat_seq[phase_starts_idx[ph]:phase_ends_idx[ph]])
      }), names(phase_def))
    }
    
    # Validate chi if non-null (must be an integer)
    if (!is.null(par$ls$chi)) {
      if (!is.numeric(par$ls$chi) || length(par$ls$chi) != 1) stop("chi must be a single numeric value.")
      if (par$ls$chi < 1) stop("chi must be at least 1.")
      if (par$ls$chi != round(par$ls$chi)) stop("chi must be a whole number.")
      if (is.null(par$ITI) && chi != "1") warning(
        "chi has no effect when ITI is NULL. ",
        "Rumination iterations are only applied when inter-trial intervals are provided via the ITI argument."
      )
    }
    
    # Validate pre-specified initial V
    if(!is.null(par$initialex)) {
      if (!is.array(initialex) || !identical(dim(initialex), c(par$D, par$C, L))) stop(paste0(
        "initialex must be a numeric array of dimension D x C x L (",
        par$D, " x ", par$C, " x ", L, "), ",
        "but has dimensions ", paste(dim(initialex), collapse = " x "), "."
      ))
      if (any(abs(initialex) > 0.5)) warning(
        "initialex contains values outside [-0.5, 0.5]. Values are on the centered reward scale ",
        "and should lie within this range unless rewards are non-binary."
      )
    }
    # ============================================================
    # 4. Storage objects
    # ============================================================
    
    arm_chosen   = numeric(par$ntrials)
    win_received = matrix(NA, nrow = par$ntrials, ncol = par$D)  # D outcomes per trial
    lsb          = matrix(NA, par$ntrials, L)
    Lmax         = numeric(par$ntrials)
    q            = numeric(par$ntrials)
    mu = array(NA, dim = c(par$ntrials, par$D, par$A))
    mu_raw = array(NA, dim = c(par$ntrials, par$D, L, par$A))
    B_history <- array(NA, dim = c(par$ntrials, par$C, par$C, par$D, L))
    
    # Belief-weighted V: trial x D x cue
    V_weighted_history <- array(NA, dim = c(par$ntrials, par$D, par$C))
    dimnames(V_weighted_history) <- list(
      trial = NULL,
      dim   = paste0("D", 1:par$D),
      cue   = paste0("cue", 1:par$C)
    )
    
    # Belief-weighted V for shown cues only: NA for cues not presented on that trial
    V_weighted_shown <- array(NA, dim = c(par$ntrials, par$D, par$C))
    dimnames(V_weighted_shown) <- list(
      trial = NULL,
      dim   = paste0("D", 1:par$D),
      cue   = paste0("cue", 1:par$C)
    )
    
    # Raw V for all latent states: trial x D x cue x state
    V_raw_history <- array(NA, dim = c(par$ntrials, par$D, par$C, L))
    dimnames(V_raw_history) <- list(
      trial = NULL,
      dim   = paste0("D", 1:par$D),
      cue   = paste0("cue", 1:par$C),
      state = paste0("state", 1:L)
    )
    
    # V for each cue type (e.g. A, AX)
    V_cue <- setNames(vector("list", length(all_cue_types)), all_cue_types)
    for (ct in all_cue_types) {
      V_cue[[ct]] <- matrix(NA_real_, nrow = 0, ncol = par$D)
      colnames(V_cue[[ct]]) <- paste0("D", 1:par$D)
    }
    
    #used to create V_cue
    type_seq_flat <- unlist(out$type_seq)
    
    # For each feature cue type, find all feature indices that are subsets of it,
    # so that later code can sum them to get the total value for that cue type.
    # ct = "A"  -> indices for ["A"]
    # ct = "AB" -> indices for ["A", "B", "AB"] (all features whose letters ⊆ {A,B})
    cue_type_cue_idx <- setNames(lapply(all_cue_types, function(ct) {
      ct_letters <- strsplit(ct, "")[[1]]
      which(sapply(par$base_cues, function(f) {
        f_letters <- strsplit(f, "")[[1]]
        all(f_letters %in% ct_letters)
      }))
    }), all_cue_types)
    
    # ============================================================
    # 5. Run simulation
    # ============================================================
    
    agent <- LatentState(par, NULL, NULL)
    
    for (t in 1:par$ntrials) {
      par$t <- t
      
      world_this <- list(
        win   = array(world$win[, , t],       dim = c(par$D, par$A)),
        c_vec = array(world$c_vec[, , , , t], dim = c(par$D, par$C, L, par$A))
      )
      
      agent <- LatentState(par, agent, world_this)
      
      # Store trial-level results
      arm_chosen[t]      <- agent$arm
      win_received[t, ]  <- agent$win   # D-length vector
      lsb[t, ]           <- agent$lsb
      Lmax[t]            <- agent$Lmax
      q[t]               <- agent$q_plot
      mu[t, , ] <- agent$mu
      mu_raw[t, , , ] <- agent$mu_arr
      B_history[t, , , , ] <- agent$B
      
      # Belief-weighted V: loop over D and C
      for (d in 1:par$D) {
        for (c in 1:par$C) {
          V_weighted_history[t, d, c] <- sum(agent$V_pre[d, c, ] * agent$lsb_pre)
        }
      }
      
      c_trial <- world$c_vec[1, , 1, , t]
      if (is.null(dim(c_trial))) dim(c_trial) <- c(par$C, par$A)
      cue_shown <- rowSums(c_trial) > 0
      for (d in 1:par$D) {
        V_weighted_shown[t, d, cue_shown] <- V_weighted_history[t, d, cue_shown]
      }
      
      V_shown <- setNames(lapply(1:par$C, function(c) {
        vals <- V_weighted_shown[, , c]  # ntrials x D, may drop to vector if D == 1
        dim(vals) <- c(par$ntrials, par$D)  # always enforce 2D shape
        shown_idx <- !is.na(vals[, 1])
        out <- vals[shown_idx, , drop = FALSE]  # prevent drop at subsetting
        if (length(dim(out)) < 2) dim(out) <- c(sum(shown_idx), par$D)  # fallback guard
        out
      }), par$base_cues)
      
      # Raw V: loop over D and L
      for (d in 1:par$D) {
        for (l in 1:L) {
          V_raw_history[t, d, , l] <- agent$V_pre[d, , l]  # C-length vector
        }
      }
      
      # V_cue: identify the chosen arm's cue segment, then append summed V per dimension
      chosen_arm_seg <- strsplit(type_seq_flat[t], "_")[[1]][agent$arm]
      if (chosen_arm_seg != "0") {
        cue_idx <- cue_type_cue_idx[[chosen_arm_seg]]
        row_d   <- sapply(1:par$D, function(d) sum(V_weighted_history[t, d, cue_idx]))
        V_cue[[chosen_arm_seg]] <- rbind(V_cue[[chosen_arm_seg]], row_d)
      }
    }
    
    # ============================================================
    # 5.5 Create long format output for ANOVA use and general output
    # ============================================================
    
    phase_ends_trial <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
    phase_starts_trial <- c(1, phase_ends_trial[-length(phase_ends_trial)] + 1)
    
    # Build long-format data for this participant using D = 1 only
    long <- vector("list", length(phase_def))
    for (phase_i in seq_along(phase_def)) {
      
      phase_rows <- do.call(rbind, lapply(all_cue_types, function(ct) {
        v_mat <- V_cue[[ct]]
        if (nrow(v_mat) == 0) return(NULL)
        
        # Find which trials in type_seq_flat contain this cue type as a chosen arm
        cue_trial_idx <- which(sapply(seq_along(type_seq_flat), function(t_idx) {
          chosen_seg <- strsplit(type_seq_flat[t_idx], "_")[[1]][arm_chosen[t_idx]]
          chosen_seg == ct
        }))
        
        # Presentation index of the last presentation before this phase starts
        pres_start <- if (phase_i == 1) 1 else
          sum(cue_trial_idx <= phase_starts_trial[phase_i] - 1) + 1
        
        # Presentation index of the last presentation within this phase
        pres_end <- sum(cue_trial_idx <= phase_ends_trial[phase_i])
        
        if (pres_start > pres_end) {
          return(NULL)
        }
        
        v_phase <- v_mat[pres_start:pres_end, 1] # D1 only
        
        #create the appropriate block vector based on block input
        block_vec <- if (!is.null(block)) {
          if (length(block) == 1) rep(block, length(phase_def)) else block
        } else NULL
        
        n_pres  <- length(v_phase)
        
        #(cue trial/total trials) * block size
        pres_per_block <- if (!is.null(block_vec)) {
          (phase_def[[phase_i]]$trials[[ct]] / 
             sum(unlist(phase_def[[phase_i]]$trials))) * block_vec[phase_i]
        } else NULL
        
        block_col <- if (!is.null(block_vec)) {
          factor(rep(seq_len(n_pres / pres_per_block), each = pres_per_block))
        } else NULL
        
        long_df <- data.frame(
          subject      = factor(i + prev),
          group        = factor(group),
          cue_type     = factor(ct),
          presentation = factor(seq_len(length(v_phase))),
          V            = v_phase
        )
        
        if (!is.null(block_col)) {
          long_df$block <- block_col
        }
        
        long_df
        
      }))
      
      long[[phase_i]] <- phase_rows
    }
    
    
    #Make V into a list format which is more user friendly
    V_list <- setNames(lapply(1:par$C, function(c) {
      mat <- V_weighted_history[, , c]
      dim(mat) <- c(par$ntrials, par$D)
      mat
    }), par$base_cues)
    
    #Make V_raw into a list format to be more user friendly
    V_raw_list <- setNames(lapply(1:par$C, function(c) {
      mat <- V_raw_history[, , c, ]
      dim(mat) <- c(par$ntrials, par$D, L)
      mat
    }), par$base_cues)
    
    all_results[[i]] <- list(
      arm_chosen   = arm_chosen,
      win_received = win_received,
      lsb          = lsb,
      Lmax         = Lmax,
      q            = q,
      mu = mu,
      mu_raw = mu_raw,
      V = V_weighted_history,
      V_list=V_list,
      V_raw_list=V_raw_list,
      V_shown = V_shown,
      V_raw      = V_raw_history,
      V_cue        = V_cue,
      B            = B_history,
      final_agent  = agent,
      par                = par,
      world              = world,
      type_seq           = out$type_seq,
      long = long
    )
  }
  
  # ============================================================
  # 6. Plot V_shown per cue (across participants if n > 1)
  # ============================================================
  base_cues <- par$base_cues
  if (plot == TRUE) {
    n_cues <- par$C
    
    plot_ids <- if (isTRUE(plot_all)) seq_len(n) else 1L
    
    for (i in plot_ids) {
      for (d in 1:par$D) {
        
        x_list <- lapply(1:n_cues, function(c) {
          1:all_results[[i]]$par$ntrials
        })
        
        y_list <- lapply(1:n_cues, function(c) {
          all_results[[i]]$V[, d, c]
        })
        
        y_all <- unlist(y_list)
        
        par(mar = c(5, 4, 4, 10))
        
        plot(NA,
             xlim = c(1, all_results[[i]]$par$ntrials),
             ylim = if (isTRUE(fixed_axis)) c(-0.5, 0.5) else range(y_all, na.rm = TRUE),             xlab = "Trial",
             ylab = "V",
             xaxt = "n",
             main = paste0("V — D", d, paste0(" — Participant ", i)))
        
        abline(v = 1:all_results[[i]]$par$ntrials, col = "gray90", lty = 1)
        
        for (c in 1:n_cues) {
          lines(x_list[[c]], y_list[[c]], col = c)
        }
        
        # Vertical lines at phase boundaries
        phase_ends <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
        abline(v = phase_ends[-length(phase_ends)] + 0.5, col = "gray40", lty = 2)
        
        legend(x = par("usr")[2],
               y = par("usr")[4],
               legend = par$base_cues,
               col    = 1:n_cues,
               lty    = 1,
               bty    = "n",
               xpd    = TRUE)
        
        if (isTRUE(plot_label)) {
          type_seq_flat <- unlist(all_results[[i]]$type_seq)
          win_flat <- all_results[[i]]$win_received  # ntrials x D
          arm_flat <- all_results[[i]]$arm_chosen
          
          # Color each label by reward on chosen arm
          label_cols <- sapply(1:all_results[[i]]$par$ntrials, function(t) {
            reward <- sum(win_flat[t, ])
            if (reward > 0) "darkgreen" else "red"
          })
          
          axis(1,
               at       = 1:all_results[[i]]$par$ntrials,
               labels   = type_seq_flat,
               las      = 2,
               cex.axis = 0.6,
               col.axis = "black",
               tick     = TRUE)
          
          # Recolor each label individually
          for (t in 1:all_results[[i]]$par$ntrials) {
            axis(1,
                 at       = t,
                 labels   = type_seq_flat[t],
                 las      = 2,
                 cex.axis = 0.6,
                 col.axis = label_cols[t],
                 tick     = FALSE)
          }
          
          if (isTRUE(plot_arm)) {
            mtext(arm_flat,
                  side = 1,
                  at   = 1:all_results[[i]]$par$ntrials,
                  line = 4.5,
                  las  = 2,
                  cex  = 0.5,
                  col  = label_cols)
          }
          
        } else {
          axis(1, at = seq(0, all_results[[i]]$par$ntrials, by = 10))
          axis(1, at = 1:all_results[[i]]$par$ntrials, labels = FALSE, tcl = -0.3)
        }
      }
    }
  }
  
  
  # ============================================================
  # 7. Plot V_shown: per-cue per-dimension, mean across participants
  # ============================================================
  if (isTRUE(plot_shown)) {
    
    # For each cue, find which presentation indices correspond to phase boundaries
    type_seq_flat <- unlist(all_results[[1]]$type_seq)
    phase_lengths <- sapply(phase_def, function(ph) sum(unlist(ph$trials)))
    phase_ends_trial <- cumsum(phase_lengths)  # trial indices where each phase ends
    
    for (d in 1:par$D) {
      for (c in 1:par$C) {
        
        # Collect V_shown[[c]][, d] across all participants, then average
        mat <- sapply(1:n, function(i) all_results[[i]]$V_shown[[c]][, d])
        if (is.null(dim(mat))) dim(mat) <- c(nrow(as.matrix(mat)), n)
        y_mean <- rowMeans(mat, na.rm = TRUE)
        y_sd   <- apply(mat, 1, function(x) sd(x, na.rm = TRUE))
        
        # Find which presentations of this cue fall at phase boundaries
        # by checking which trial indices (from type_seq) contain this cue
        cue_trial_idx <- which(sapply(type_seq_flat, function(nm) {
          segs <- strsplit(nm, "_")[[1]]
          cue_letters <- strsplit(base_cues[c], "")[[1]]
          any(sapply(segs, function(seg) {
            seg_letters <- strsplit(seg, "")[[1]]
            all(cue_letters %in% seg_letters)
          }))
        }))
        
        
        # Phase boundary presentation indices: last presentation before each phase end
        phase_boundary_pres <- sapply(phase_ends_trial[-length(phase_ends_trial)], function(end_t) {
          sum(cue_trial_idx <= end_t)
        })
        phase_boundary_pres <- phase_boundary_pres[phase_boundary_pres > 0 &
                                                     phase_boundary_pres < nrow(mat)]
        
        y_lo <- y_mean - y_sd
        y_hi <- y_mean + y_sd
        ylim_range <- if (isTRUE(plot_sd) && n > 1) {
          sd_range <- range(c(y_lo, y_hi), na.rm = TRUE)
          if (isTRUE(fixed_axis)) c(min(-0.5, sd_range[1]), max(0.5, sd_range[2])) else sd_range
        } else if (isTRUE(fixed_axis)) c(-0.5, 0.5) else range(mat, na.rm = TRUE)
        plot(NA,
             xlim = c(1, nrow(mat)),
             ylim = ylim_range,
             xlab = "Presentation number",
             ylab = "V (±SD)",
             xaxt = "n",
             main = paste0("V shown — D", d, " — Cue ", par$base_cues[c],
                           if (n > 1) paste0(" (mean of ", n, " participants)") else ""))
        
        # Faint vertical grid lines
        abline(v = 1:nrow(mat), col = "gray90", lty = 1)
        
        # SE ribbon
        if (isTRUE(plot_sd) && n > 1) {
          xs <- c(1:nrow(mat), nrow(mat):1)
          ys <- c(y_hi, rev(y_lo))
          polygon(xs, ys, col = adjustcolor("black", alpha.f = 0.15), border = NA)
        }
        
        # Mean line
        lines(1:nrow(mat), y_mean, col = "black", lwd = 2)
        
        # Phase boundary lines
        if (length(phase_boundary_pres) > 0) {
          abline(v = phase_boundary_pres + 0.5, col = "gray40", lty = 2)
        }
        
        axis(1, at = seq(0, nrow(mat), by = 10))
        axis(1, at = 1:nrow(mat), labels = FALSE, tcl = -0.3)
      }
    }
  }
  # ============================================================
  # 8. Plot V_raw_history: per-latent-state
  # ============================================================
  if (isTRUE(plot_raw)) {
    
    plot_ids_raw <- if (isTRUE(plot_raw_all)) seq_len(n) else 1L
    
    for (i in plot_ids_raw) {
      
      res      <- all_results[[i]]
      n_states <- max(res$Lmax)
      V_raw    <- res$V_raw      # ntrials x D x C x L
      
      for (l in 1:n_states) {
        for (d in 1:res$par$D) {
          
          # Extract V for this state and dimension: ntrials x C
          v_mat <- V_raw[, d, , l]
          if (is.null(dim(v_mat))) dim(v_mat) <- c(res$par$ntrials, res$par$C)
          
          # Belief proportion for this state across trials
          lsb_state <- res$lsb[, l]
          
          # First trial where belief > 0 for this state
          first_active <- which(lsb_state > 0)[1]
          
          y_all <- as.vector(v_mat)
          
          par(mar = if (isTRUE(plot_label)) c(7, 4, 4, 12) else c(5, 4, 4, 12))
          
          plot(NA,
               xlim = c(1, res$par$ntrials),
               ylim = if (isTRUE(fixed_axis)) c(-0.5, 0.5) else range(y_all, na.rm = TRUE),               xlab = "Trial",
               ylab = "V (raw)",
               xaxt = "n",
               main = paste0("Latent State ", l, " — D", d, " — Participant ", i))
          
          abline(v = 1:res$par$ntrials, col = "gray90", lty = 1)
          
          for (c in 1:res$par$C) {
            idx <- first_active:res$par$ntrials
            lines(idx, v_mat[idx, c], col = c)
          }
          
          # Belief proportion line on second y-axis
          lsb_scaled <- range(y_all, na.rm = TRUE)[1] +
            lsb_state * diff(range(y_all, na.rm = TRUE))
          lines(1:res$par$ntrials, lsb_scaled, col = "black", lty = 2, lwd = 1.5)
          axis(4, at = range(y_all, na.rm = TRUE)[1] +
                 seq(0, 1, by = 0.25) * diff(range(y_all, na.rm = TRUE)),
               labels = paste0(seq(0, 100, by = 25), "%"),
               las = 2)
          mtext("Belief (%)", side = 4, line = 10, xpd = NA)
          
          # Vertical lines at phase boundaries
          phase_ends <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
          abline(v = phase_ends[-length(phase_ends)] + 0.5, col = "gray40", lty = 2)
          
          legend(x = par("usr")[2],
                 y = par("usr")[4],
                 legend = c(res$par$base_cues, "Belief"),
                 col    = c(1:res$par$C, "black"),
                 lty    = c(rep(1, res$par$C), 2),
                 bty    = "n",
                 xpd    = TRUE)
          
          if (isTRUE(plot_label)) {
            type_seq_flat_raw <- unlist(res$type_seq)
            win_flat <- res$win_received
            arm_flat_raw <- res$arm_chosen
            
            label_cols_raw <- sapply(1:res$par$ntrials, function(t) {
              reward <- sum(win_flat[t, ])
              if (reward > 0) "darkgreen" else "red"
            })
            
            axis(1,
                 at       = 1:res$par$ntrials,
                 labels   = type_seq_flat_raw,
                 las      = 2,
                 cex.axis = 0.6,
                 col.axis = "black",
                 tick     = TRUE)
            
            for (t in 1:res$par$ntrials) {
              axis(1,
                   at       = t,
                   labels   = type_seq_flat_raw[t],
                   las      = 2,
                   cex.axis = 0.6,
                   col.axis = label_cols_raw[t],
                   tick     = FALSE)
            }
            
            if (isTRUE(plot_arm)) {
              mtext(arm_flat_raw,
                    side = 1,
                    at   = 1:res$par$ntrials,
                    line = 4.5,
                    las  = 2,
                    cex  = 0.5,
                    col  = label_cols_raw)
            }
            
          } else {
            axis(1, at = seq(0, res$par$ntrials, by = 10))
            axis(1, at = 1:res$par$ntrials, labels = FALSE, tcl = -0.3)
          }
        }
      }
    }
  }
  
  # ============================================================
  # 9. Plot V_cue: per-cue-type, mean across participants
  # ============================================================
  if (isTRUE(plot_cue)) {
    
    type_seq_flat <- unlist(all_results[[1]]$type_seq)
    
    for (ct in all_cue_types) {
      for (d in 1:par$D) {
        
        # Collect V_cue[[ct]][, d] across all participants
        raw_list <- lapply(1:n, function(i) all_results[[i]]$V_cue[[ct]][, d])
        
        # Find the expected number of rows (max across participants)
        n_rows <- max(sapply(raw_list, length))
        
        if (n_rows == 0) {
          warning(paste0("Cue type '", ct, "' is not plotted because it was never chosen."))
          next
        }
        
        # Pad shorter participants with NA so matrix is rectangular
        mat <- sapply(raw_list, function(x) {
          if (length(x) < n_rows) c(x, rep(NA_real_, n_rows - length(x))) else x
        })
        if (is.null(dim(mat))) dim(mat) <- c(length(mat) / n, n)
        
        y_mean <- rowMeans(mat, na.rm = TRUE)
        y_sd   <- apply(mat, 1, function(x) sd(x, na.rm = TRUE))
        
        # Number of participants contributing to each datapoint
        n_contrib <- apply(mat, 1, function(x) sum(!is.na(x)))
        
        y_lo <- y_mean - y_sd
        y_hi <- y_mean + y_sd
        ylim_range <- if (isTRUE(plot_sd) && n > 1) {
          sd_range <- range(c(y_lo, y_hi), na.rm = TRUE)
          if (isTRUE(fixed_axis)) c(min(-0.5, sd_range[1]), max(0.5, sd_range[2])) else sd_range
        } else if (isTRUE(fixed_axis)) c(-0.5, 0.5) else range(mat, na.rm = TRUE)
        
        # Add bottom margin for participant counts if multi-arm
        if (par$A > 1) {
          par(mar = c(7, 4, 4, 2))
        } else {
          par(mar = c(5, 4, 4, 2))
        }
        
        plot(NA,
             xlim = c(1, nrow(mat)),
             ylim = ylim_range,
             xlab = "Presentation number",
             ylab = "V (±SD)",
             xaxt = "n",
             main = paste0("V cue — D", d, " — Cue type ", ct,
                           if (n > 1) paste0(" (mean of ", n, " participants)") else ""))
        
        # Faint vertical grid lines
        abline(v = 1:nrow(mat), col = "gray90", lty = 1)
        
        # SE ribbon — trim to non-NA range before drawing
        if (isTRUE(plot_sd) && n > 1) {
          valid <- which(!is.na(y_lo) & !is.na(y_hi))
          if (length(valid) > 0) {
            xs <- c(valid, rev(valid))
            ys <- c(y_hi[valid], rev(y_lo[valid]))
            polygon(xs, ys, col = adjustcolor("black", alpha.f = 0.15), border = NA)
          }
        }
        
        # Mean line
        lines(1:nrow(mat), y_mean, col = "black", lwd = 2)
        
        # Phase boundary lines
        #find which trials each cue type appears in, considering only the chosen arm
        cue_trial_idx <- which(sapply(seq_along(type_seq_flat), function(t_idx) {
          chosen_seg <- strsplit(type_seq_flat[t_idx], "_")[[1]][all_results[[1]]$arm_chosen[t_idx]]
          chosen_seg == ct
        }))
        
        phase_ends_trial <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
        phase_boundary_pres <- sapply(phase_ends_trial[-length(phase_ends_trial)], function(end_t) {
          sum(cue_trial_idx <= end_t)
        })
        phase_boundary_pres <- phase_boundary_pres[phase_boundary_pres > 0 &
                                                     phase_boundary_pres < nrow(mat)]
        if (length(phase_boundary_pres) > 0) {
          abline(v = phase_boundary_pres + 0.5, col = "gray40", lty = 2)
        }
        
        axis(1, at = seq(0, nrow(mat), by = 10))
        axis(1, at = 1:nrow(mat), labels = FALSE, tcl = -0.3)
        
        # Participant count line below x-axis (multi-arm only)
        if (par$A > 1) {
          mtext(n_contrib,
                side = 1,
                at   = 1:nrow(mat),
                line = 4,
                las  = 2,
                cex  = 0.6)
          mtext("N",
                side = 1,
                at   = 0,
                line = 4,
                las  = 2,
                cex  = 0.6,
                font = 2)
        }
      }
    }
  }
  
  # ============================================================
  # 10. Plot mu of chosen arm
  # ============================================================
  if (isTRUE(plot_mu)) {
    
    plot_ids_mu <- if (isTRUE(plot_all)) seq_len(n) else 1L
    
    for (i in plot_ids_mu) {
      
      res        <- all_results[[i]]
      type_seq_flat_mu <- unlist(res$type_seq)
      arm_flat   <- res$arm_chosen
      
      for (d in 1:res$par$D) {
        
        # Extract mu of chosen arm for each trial: ntrials-length vector
        mu_chosen <- sapply(1:res$par$ntrials, function(t) {
          res$mu[t, d, arm_flat[t]]
        })
        
        y_all <- mu_chosen
        
        ylim_range <- if (isTRUE(fixed_axis)) c(-0.5, 0.5) else range(y_all, na.rm = TRUE)
        
        par(mar = c(5, 4, 4, 4))
        
        plot(NA,
             xlim = c(1, res$par$ntrials),
             ylim = ylim_range,
             xlab = "Trial",
             ylab = "mu (chosen arm)",
             xaxt = "n",
             main = paste0("mu chosen arm — D", d, " — Participant ", i))
        
        abline(v = 1:res$par$ntrials, col = "gray90", lty = 1)
        
        lines(1:res$par$ntrials, mu_chosen, col = "black", lwd = 2)
        
        # Phase boundary lines
        phase_ends <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
        abline(v = phase_ends[-length(phase_ends)] + 0.5, col = "gray40", lty = 2)
        
        # plot_label: trial type labels on x-axis, colored by reward
        if (isTRUE(plot_label)) {
          win_flat <- res$win_received  # ntrials x D
          
          label_cols <- sapply(1:res$par$ntrials, function(t) {
            reward <- sum(win_flat[t, ])
            if (reward > 0) "darkgreen" else "red"
          })
          
          axis(1,
               at       = 1:res$par$ntrials,
               labels   = type_seq_flat_mu,
               las      = 2,
               cex.axis = 0.6,
               col.axis = "black",
               tick     = TRUE)
          
          for (t in 1:res$par$ntrials) {
            axis(1,
                 at       = t,
                 labels   = type_seq_flat_mu[t],
                 las      = 2,
                 cex.axis = 0.6,
                 col.axis = label_cols[t],
                 tick     = FALSE)
          }
          
          if (isTRUE(plot_arm)) {
            mtext(arm_flat,
                  side = 1,
                  at   = 1:res$par$ntrials,
                  line = 4.5,
                  las  = 2,
                  cex  = 0.5,
                  col  = label_cols)
          }
          
        } else {
          axis(1, at = seq(0, res$par$ntrials, by = 10))
          axis(1, at = 1:res$par$ntrials, labels = FALSE, tcl = -0.3)
        }
      }
    }
  }
  
  # ============================================================
  # 11. Plot q (Page's statistic) across trials
  #      q is plotted after reward, before state creation/reset.
  # ============================================================
  if (isTRUE(plot_q)) {
    
    plot_ids_q <- if (isTRUE(plot_all)) seq_len(n) else 1L
    
    for (i in plot_ids_q) {
      res <- all_results[[i]]
      
      q_vec   <- res$q
      
      # simple plotting-only cap for infinities
      inf_trials <- which(is.infinite(q_vec))
      
      if (any(is.infinite(q_vec))) {
        q_vec[is.infinite(q_vec) & q_vec > 0] <-  as.numeric(eta)*2
        q_vec[is.infinite(q_vec) & q_vec < 0] <-  as.numeric(eta)*2
      }
      
      ntr     <- res$par$ntrials
      eta_val <- res$par$ls$eta
      
      # y-lim from 0 to max(q, eta), with small headroom
      q_max <- max(q_vec, eta_val, na.rm = TRUE)
      y_pad <- 0.05 * (q_max + 1e-8)
      ylim_range <- c(0, q_max + y_pad)
      
      type_seq_flat_q <- unlist(res$type_seq)
      win_flat_q      <- res$win_received
      
      par(mar = if (isTRUE(plot_label)) c(7, 4, 4, 2) else c(5, 4, 4, 2))
      
      plot(NA,
           xlim = c(1, ntr),
           ylim = ylim_range,
           xlab = "Trial",
           ylab = "q (change-point statistic)",
           xaxt = "n",
           main = paste0("q — Participant ", i))
      
      # Grid
      abline(v = 1:ntr, col = "gray90", lty = 1)
      
      # q trajectory
      lines(1:ntr, q_vec, col = "black", lwd = 2)
      
      # Mark trials where q was Inf (true value arbitrarily large, capped at 
      # eta*2 for display)
      if (length(inf_trials) > 0) {
        points(inf_trials, q_vec[inf_trials],
               pch = 17, col = "red", cex = 1.2)
        text(inf_trials, q_vec[inf_trials],
             labels = "\u221e",
             pos = 3, col = "red", cex = 0.9)
        legend("topleft",
               legend = "q = \u221e (capped at 2\u00d7\u03b7 for display)",
               pch = 17, col = "red", bty = "n", cex = 0.8)
      }
      
      # Threshold eta
      abline(h = eta_val, col = "red", lty = 2, lwd = 1.5)
      
      # Phase boundaries
      phase_ends <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
      if (length(phase_ends) > 1) {
        abline(v = phase_ends[-length(phase_ends)] + 0.5, col = "gray40", lty = 2)
      }
      
      # X‑axis labels with trial types / reward coloring, mirroring other plots
      if (isTRUE(plot_label)) {
        
        label_cols <- sapply(1:ntr, function(t) {
          reward <- sum(win_flat_q[t, ])
          if (reward > 0) "darkgreen" else "red"
        })
        
        axis(1,
             at       = 1:ntr,
             labels   = type_seq_flat_q,
             las      = 2,
             cex.axis = 0.6,
             col.axis = "black",
             tick     = TRUE)
        
        # Recolor each label individually
        for (t in 1:ntr) {
          axis(1,
               at       = t,
               labels   = type_seq_flat_q[t],
               las      = 2,
               cex.axis = 0.6,
               col.axis = label_cols[t],
               tick     = FALSE)
        }
        
        if (isTRUE(plot_arm)) {
          mtext(res$arm_chosen,
                side = 1,
                at   = 1:ntr,
                line = 4.5,
                las  = 2,
                cex  = 0.5,
                col  = label_cols)
        }
        
      } else {
        axis(1, at = seq(0, ntr, by = 10))
        axis(1, at = 1:ntr, labels = FALSE, tcl = -0.3)
      }
    }
  }
  
  all_results <- lapply(all_results, function(res) {
    res$V        <- res$V_list
    res$V_raw    <- res$V_raw_list
    res$V_list   <- NULL
    res$V_raw_list <- NULL
    res
  })
  
  return(all_results)
}

#when Mauchly's test is significant, I should use the GG sphericity correction values
#this function also checks if the significance (yes or no) matches what is expected
check_p <- function(ez, effect, sig) {
  mauchly <- ez$`Mauchly's Test for Sphericity`
  
  p <- if (!is.null(mauchly) && effect %in% mauchly$Effect &&
           mauchly %>% filter(Effect == effect) %>% pull(`p<.05`) == "*") {
    ez$`Sphericity Corrections` %>% filter(Effect == effect) %>% pull(`p[GG]`)
  } else {
    ez$ANOVA %>% filter(Effect == effect) %>% pull(p)
  }
  
  pass <<- c(pass, (p < 0.05) == sig)
  
  return(p)
}

#for non-ANOVA, this function will check if the significance of the test matches
#the intended significance (true or false)
check_p_other <- function(p, sig) {
  pass <<- c(pass, (p < 0.05) == sig)
  return(p)
}


#decide if I want to make plot_avg and maybe plot_raw_avg, each of which plot
#the mean V of each cue at every trial across participants, ignoring the fact 
#that trial order for participants differed

#also decide if I want to bother adding plot_label to plot_compare when the 
#inserted groups have identical trial sequences