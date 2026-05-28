plot_compare <- function(...,
                         plot = FALSE,
                         plot_avg = FALSE,
                         fixed_axis = FALSE,
                         plot_sd = FALSE,
                         plot_label = FALSE,
                         plot_shown = FALSE,
                         plot_cue = FALSE) {
  
  groups     <- list(...)
  group_names <- names(groups)
  if (is.null(group_names)) group_names <- paste0("group", seq_along(groups))
  
  n_groups <- length(groups)
  colors   <- seq_len(n_groups)
  
  # Derive shared structure from first group
  first   <- groups[[1]][[1]]
  par     <- first$par
  ntrials <- first$par$ntrials
  phase_ends <- cumsum(sapply(first$type_seq, function(ph) length(unlist(ph))))
  
  # Recompute phase ends from type_seq
  phase_ends <- cumsum(sapply(first$type_seq, function(ph) length(unlist(ph))))
  
  # Validate if type sequences match across all groups if using plot_label
  if (isTRUE(plot_label)) {
    ref_seq <- unlist(first$type_seq)
    for (g in seq_len(n_groups)) {
      for (j in seq_along(groups[[g]])) {
        cmp_seq <- unlist(groups[[g]][[j]]$type_seq)
        if (!identical(ref_seq, cmp_seq)) stop(paste0(
          "plot_label = TRUE requires all groups to have identical type sequences, ",
          "but group '", group_names[g], "', participant ", j, " differs from the reference."
        ))
      }
    }
  }
  
  # Helper: compute mean V across participants in one run result list
  mean_V <- function(res_list, d, c) {
    mat <- sapply(res_list, function(res) res$V[[c]][, d])
    if (is.null(dim(mat))) dim(mat) <- c(length(mat), 1)
    list(
      mean = rowMeans(mat, na.rm = TRUE),
      sd   = apply(mat, 1, function(x) sd(x, na.rm = TRUE)),
      n    = length(res_list)
    )
  }
  
  # ============================================================
  # plot=TRUE equivalent: V of participant 1 per group
  # ============================================================
  if (plot == TRUE) {
    for (d in 1:par$D) {
      
      all_vals <- unlist(lapply(seq_len(n_groups), function(g) {
        lapply(1:par$C, function(c) groups[[g]][[1]]$V[[c]][, d])
      }))
      
      ylim_range <- if (isTRUE(fixed_axis)) c(-0.5, 0.5) else range(all_vals, na.rm = TRUE)
      
      par(mar = c(5, 4, 4, 12))
      
      plot(NA,
           xlim = c(1, ntrials),
           ylim = ylim_range,
           xlab = "Trial",
           ylab = "V",
           xaxt = "n",
           main = paste0("V — D", d, " — Participant 1"))
      
      abline(v = 1:ntrials, col = "gray90", lty = 1)
      abline(v = phase_ends[-length(phase_ends)] + 0.5, col = "gray40", lty = 2)
      
      for (c in 1:par$C) {
        for (g in seq_len(n_groups)) {
          lines(1:ntrials, groups[[g]][[1]]$V[[c]][, d], col = c, lty = g, lwd = 2)
        }
      }
      
      legend(x = par("usr")[2],
             y = par("usr")[4],
             legend = c(first$par$base_cues, group_names),
             col    = c(1:par$C, rep("black", n_groups)),
             lty    = c(rep(1, par$C), seq_len(n_groups)),
             lwd    = 2,
             bty    = "n",
             xpd    = TRUE)
      
      if (isTRUE(plot_label)) {
        type_seq_flat <- unlist(first$type_seq)
        win_flat      <- first$win_received
        arm_flat      <- first$arm_chosen
        
        label_cols <- sapply(1:ntrials, function(t) {
          reward <- sum(win_flat[t, ])
          if (reward > 0) "darkgreen" else "red"
        })
        
        axis(1,
             at       = 1:ntrials,
             labels   = type_seq_flat,
             las      = 2,
             cex.axis = 0.6,
             col.axis = "black",
             tick     = TRUE)
        
        for (t in 1:ntrials) {
          axis(1,
               at       = t,
               labels   = type_seq_flat[t],
               las      = 2,
               cex.axis = 0.6,
               col.axis = label_cols[t],
               tick     = FALSE)
        }
      } else {
        axis(1, at = seq(0, ntrials, by = 10))
        axis(1, at = 1:ntrials, labels = FALSE, tcl = -0.3)
      }
    }
  }
  
  # ============================================================
  # plot_avg=TRUE: mean V across participants per group
  # ============================================================
  if (isTRUE(plot_avg)) {
    for (d in 1:par$D) {
      
      all_means <- lapply(seq_len(par$C), function(c) {
        lapply(seq_len(n_groups), function(g) mean_V(groups[[g]], d, c))
      })
      
      if (isTRUE(fixed_axis)) {
        ylim_range <- c(-0.5, 0.5)
        if (isTRUE(plot_sd)) {
          all_vals <- unlist(lapply(all_means, function(cue) lapply(cue, function(m) c(m$mean + m$sd, m$mean - m$sd))))
          sd_range <- range(all_vals, na.rm = TRUE)
          ylim_range <- c(min(ylim_range[1], sd_range[1]), max(ylim_range[2], sd_range[2]))
        }
      } else {
        all_vals <- unlist(lapply(all_means, function(cue) lapply(cue, function(m) {
          if (isTRUE(plot_sd) && m$n > 1) c(m$mean + m$sd, m$mean - m$sd) else m$mean
        })))
        ylim_range <- range(all_vals, na.rm = TRUE)
      }
      
      par(mar = c(5, 4, 4, 12))
      
      plot(NA,
           xlim = c(1, ntrials),
           ylim = ylim_range,
           xlab = "Trial",
           ylab = "V",
           xaxt = "n",
           main = paste0("V — D", d, " (group means)"))
      
      abline(v = 1:ntrials, col = "gray90", lty = 1)
      abline(v = phase_ends[-length(phase_ends)] + 0.5, col = "gray40", lty = 2)
      
      for (c in 1:par$C) {
        for (g in seq_len(n_groups)) {
          m <- all_means[[c]][[g]]
          
          if (isTRUE(plot_sd) && m$n > 1) {
            xs <- c(1:ntrials, ntrials:1)
            ys <- c(m$mean + m$sd, rev(m$mean - m$sd))
            polygon(xs, ys, col = adjustcolor(colors[g], alpha.f = 0.15), border = NA)
          }
          
          lines(1:ntrials, m$mean, col = c, lty = g, lwd = 2)
        }
      }
      
      legend(x = par("usr")[2],
             y = par("usr")[4],
             legend = c(first$par$base_cues, group_names),
             col    = c(1:par$C, rep("black", n_groups)),
             lty    = c(rep(1, par$C), seq_len(n_groups)),
             lwd    = 2,
             bty    = "n",
             xpd    = TRUE)
      
      axis(1, at = seq(0, ntrials, by = 10))
      axis(1, at = 1:ntrials, labels = FALSE, tcl = -0.3)
    }
  }
  
  # ============================================================
  # plot_shown equivalent: V per presentation of each cue
  # ============================================================
  if (isTRUE(plot_shown)) {
    
    type_seq_flat <- unlist(first$type_seq)
    base_cues     <- first$par$base_cues
    
    for (d in 1:par$D) {
      for (c in 1:par$C) {
        
        all_shown <- lapply(seq_len(n_groups), function(g) {
          mat <- sapply(groups[[g]], function(res) res$V_shown[[c]][, d])
          if (is.null(dim(mat))) dim(mat) <- c(length(mat), 1)
          list(
            mean = rowMeans(mat, na.rm = TRUE),
            sd   = apply(mat, 1, function(x) sd(x, na.rm = TRUE)),
            n    = length(groups[[g]])
          )
        })
        
        n_pres <- length(all_shown[[1]]$mean)
        
        if (isTRUE(fixed_axis)) {
          ylim_range <- c(-0.5, 0.5)
          if (isTRUE(plot_sd)) {
            all_vals <- unlist(lapply(all_shown, function(m) c(m$mean + m$sd, m$mean - m$sd)))
            sd_range <- range(all_vals, na.rm = TRUE)
            ylim_range <- c(min(ylim_range[1], sd_range[1]), max(ylim_range[2], sd_range[2]))
          }
        } else {
          all_vals <- unlist(lapply(all_shown, function(m) {
            if (isTRUE(plot_sd) && m$n > 1) c(m$mean + m$sd, m$mean - m$sd) else m$mean
          }))
          ylim_range <- range(all_vals, na.rm = TRUE)
        }
        
        par(mar = c(5, 4, 4, 12))
        
        plot(NA,
             xlim = c(1, n_pres),
             ylim = ylim_range,
             xlab = "Presentation number",
             ylab = "V (±SD)",
             xaxt = "n",
             main = paste0("V shown — D", d, " — Cue ", first$par$base_cues[c]))
        
        abline(v = 1:n_pres, col = "gray90", lty = 1)
        
        # Phase boundary presentations for this cue
        cue_trial_idx <- which(sapply(type_seq_flat, function(nm) {
          segs <- strsplit(nm, "_")[[1]]
          cue_letters <- strsplit(base_cues[c], "")[[1]]
          any(sapply(segs, function(seg) {
            seg_letters <- strsplit(seg, "")[[1]]
            all(cue_letters %in% seg_letters)
          }))
        }))
        
        phase_ends_trial <- cumsum(sapply(first$type_seq, function(ph) length(unlist(ph))))
        phase_boundary_pres <- cumsum(sapply(seq_along(phase_ends_trial), function(ph) {
          t_start <- if (ph == 1) 1 else phase_ends_trial[ph - 1] + 1
          t_end   <- phase_ends_trial[ph]
          sum(cue_trial_idx >= t_start & cue_trial_idx <= t_end)
        }))
        phase_boundary_pres <- phase_boundary_pres[-length(phase_boundary_pres)]
        phase_boundary_pres <- phase_boundary_pres[phase_boundary_pres > 0 &
                                                     phase_boundary_pres < n_pres]
        
        if (length(phase_boundary_pres) > 0) {
          abline(v = phase_boundary_pres + 0.5, col = "gray40", lty = 2)
        }
        
        for (g in seq_len(n_groups)) {
          m <- all_shown[[g]]
          
          if (isTRUE(plot_sd) && m$n > 1) {
            xs <- c(1:n_pres, n_pres:1)
            ys <- c(m$mean + m$sd, rev(m$mean - m$sd))
            polygon(xs, ys, col = adjustcolor(colors[g], alpha.f = 0.15), border = NA)
          }
          
          lines(1:n_pres, m$mean, col = colors[g], lwd = 2)
        }
        
        legend(x = par("usr")[2],
               y = par("usr")[4],
               legend = group_names,
               col    = colors,
               lty    = 1,
               lwd    = 2,
               bty    = "n",
               xpd    = TRUE)
        
        axis(1, at = seq(0, n_pres, by = 10))
        axis(1, at = 1:n_pres, labels = FALSE, tcl = -0.3)
      }
    }
  }
  
  # ============================================================
  # plot_cue equivalent: V_cue per presentation of each cue type
  # ============================================================
  if (isTRUE(plot_cue)) {
    
    all_cue_types <- names(first$V_cue)
    
    for (ct in all_cue_types) {
      for (d in 1:par$D) {
        
        all_cue <- lapply(seq_len(n_groups), function(g) {
          raw_list <- lapply(groups[[g]], function(res) res$V_cue[[ct]][, d])
          n_rows   <- max(sapply(raw_list, length))
          if (n_rows == 0) return(NULL)
          mat <- sapply(raw_list, function(x) {
            if (length(x) < n_rows) c(x, rep(NA_real_, n_rows - length(x))) else x
          })
          if (is.null(dim(mat))) dim(mat) <- c(n_rows, 1)
          list(
            mean = rowMeans(mat, na.rm = TRUE),
            sd   = apply(mat, 1, function(x) sd(x, na.rm = TRUE)),
            n    = length(groups[[g]])
          )
        })
        
        if (any(sapply(all_cue, is.null))) next
        
        n_pres <- length(all_cue[[1]]$mean)
        
        if (isTRUE(fixed_axis)) {
          ylim_range <- c(-0.5, 0.5)
          if (isTRUE(plot_sd)) {
            all_vals <- unlist(lapply(all_cue, function(m) c(m$mean + m$sd, m$mean - m$sd)))
            sd_range <- range(all_vals, na.rm = TRUE)
            ylim_range <- c(min(-0.5, sd_range[1]), max(0.5, sd_range[2]))
          }
        } else {
          all_vals <- unlist(lapply(all_cue, function(m) c(m$mean + m$sd, m$mean - m$sd)))
          ylim_range <- range(all_vals, na.rm = TRUE)
        }
        
        par(mar = c(5, 4, 4, 12))
        
        plot(NA,
             xlim = c(1, n_pres),
             ylim = ylim_range,
             xlab = "Presentation number",
             ylab = "V (±SD)",
             xaxt = "n",
             main = paste0("V cue — D", d, " — Cue type ", ct))
        
        type_seq_flat <- unlist(first$type_seq)
        cue_trial_idx <- which(sapply(seq_along(type_seq_flat), function(t_idx) {
          chosen_seg <- strsplit(type_seq_flat[t_idx], "_")[[1]][first$arm_chosen[t_idx]]
          chosen_seg == ct
        }))
        phase_ends_trial <- cumsum(sapply(first$type_seq, function(ph) length(unlist(ph))))
        phase_boundary_pres <- cumsum(sapply(seq_along(phase_ends_trial), function(ph) {
          t_start <- if (ph == 1) 1 else phase_ends_trial[ph - 1] + 1
          t_end   <- phase_ends_trial[ph]
          sum(cue_trial_idx >= t_start & cue_trial_idx <= t_end)
        }))
        phase_boundary_pres <- phase_boundary_pres[-length(phase_boundary_pres)]
        phase_boundary_pres <- phase_boundary_pres[phase_boundary_pres > 0 &
                                                     phase_boundary_pres < n_pres]
        if (length(phase_boundary_pres) > 0) {
          abline(v = phase_boundary_pres + 0.5, col = "gray40", lty = 2)
        }
        
        abline(v = 1:n_pres, col = "gray90", lty = 1)
        
        for (g in seq_len(n_groups)) {
          m <- all_cue[[g]]
          
          if (isTRUE(plot_sd) && m$n > 1) {
            valid <- which(!is.na(m$mean + m$sd))
            if (length(valid) > 0) {
              xs <- c(valid, rev(valid))
              ys <- c((m$mean + m$sd)[valid], rev((m$mean - m$sd)[valid]))
              polygon(xs, ys, col = adjustcolor(colors[g], alpha.f = 0.15), border = NA)
            }
          }
          
          lines(1:n_pres, m$mean, col = colors[g], lwd = 2)
        }
        
        legend(x = par("usr")[2],
               y = par("usr")[4],
               legend = group_names,
               col    = colors,
               lty    = 1,
               lwd    = 2,
               bty    = "n",
               xpd    = TRUE)
        
        axis(1, at = seq(0, n_pres, by = 10))
        axis(1, at = 1:n_pres, labels = FALSE, tcl = -0.3)
      }
    }
  }
  
  invisible(NULL)
}















#manual world creator
make_world <- function(trial_patterns, reward, L = 10, features = NULL) {
  
  if (!is.character(trial_patterns)) stop(
    "trial_patterns must be a character vector of cue patterns, for example c(\"A\", \"B\", \"A\") or c(\"AX_B\", \"B_C\")."
  )
  if (any(grepl("[^A-Za-z_0]", trial_patterns))) stop(
    "trial_patterns may only contain letters, underscores, and 0. ",
    "Examples include \"A\", \"AX_B\", \"A_0\", and \"A_BC_X\"."
  )
  if (any(grepl("^_|_$|__", trial_patterns))) stop(
    "trial_patterns contains invalid underscore placement. ",
    "Underscores must separate arms and cannot appear at the start, end, or consecutively."
  )
  
  ntrials <- length(trial_patterns)
  
  # Determine number of arms from first trial
  A <- length(strsplit(trial_patterns[1], "_")[[1]])
  
  # Validate all trials have same number of arms
  arm_counts <- sapply(trial_patterns, function(tp) length(strsplit(tp, "_")[[1]]))
  if (any(arm_counts != A)) stop(paste0(
    "All trial patterns must have the same number of arms (", A, "), ",
    "but trials ", paste(which(arm_counts != A), collapse = ", "),
    " have ", arm_counts[arm_counts != A], " arm(s)."
  ))
  
  # Determine C from unique letters across all patterns. [_0] removes both _ and 0
  # Determine all unique arm segments (non-nil) across all trial patterns
  all_arm_segs <- unique(Filter(function(s) s != "0",
                                unlist(lapply(trial_patterns, function(tp) strsplit(tp, "_")[[1]]))
  ))
  
  # Valid atomic letters that appear in at least one arm segment
  valid_atoms <- sort(unique(unlist(strsplit(all_arm_segs, ""))))
  
  if (!is.null(features)) {
    
    if (!is.character(features)) stop(
      "features must be NULL or a character vector."
    )
    if (length(features) == 0L) stop(
      "features must be NULL or a non-empty character vector."
    )
    if (anyDuplicated(features)) {
      warning("features contains duplicate entries; duplicates will be ignored.")
      features <- unique(features)
    }
    
    # Validate each feature: all its letters must exist, and its letters must
    # be a subset of at least one arm segment's letters
    invalid_features <- character(0)
    for (f in features) {
      f_letters <- strsplit(f, "")[[1]]
      if (!all(f_letters %in% valid_atoms)) {
        invalid_features <- c(invalid_features, f)
        next
      }
      appears <- any(sapply(all_arm_segs, function(seg) {
        seg_letters <- strsplit(seg, "")[[1]]
        all(f_letters %in% seg_letters)
      }))
      if (!appears) invalid_features <- c(invalid_features, f)
    }
    if (length(invalid_features) > 0L) stop(
      "features contains entries whose letters never co-occur in any arm segment: ",
      paste(invalid_features, collapse = ", "), "."
    )
    
    # Validate coverage: every arm segment must have at least one feature active
    for (seg in all_arm_segs) {
      seg_letters <- strsplit(seg, "")[[1]]
      covered <- any(sapply(features, function(f) {
        f_letters <- strsplit(f, "")[[1]]
        all(f_letters %in% seg_letters)
      }))
      if (!covered) stop(
        "features does not cover arm segment '", seg, "'. ",
        "At least one feature whose letters are a subset of {",
        paste(seg_letters, collapse = ", "), "} must be included."
      )
    }
    
    base_cues <- features
    
  } else {
    
    # Default: individual letters + full multi-letter arm segments
    config_cues <- sort(setdiff(all_arm_segs, valid_atoms))
    base_cues   <- c(valid_atoms, config_cues)
    
  }
  
  C <- length(base_cues)
  
  # Normalize reward: allow plain vector (D=1) or list of vectors
  if (is.numeric(reward)) {
    reward <- list(reward)
  }
  if (!is.list(reward)) stop(
    "reward must be a numeric vector (for D=1) or a list of numeric vectors (one per dimension)."
  )
  
  # Infer D from number of reward vectors
  D <- length(reward)
  
  # Validate each reward vector has length ntrials
  for (d in 1:D) {
    if (!is.numeric(reward[[d]])) stop(paste0(
      "reward[[", d, "]] must be a numeric vector."
    ))
    if (length(reward[[d]]) != ntrials) stop(paste0(
      "reward[[", d, "]] must have length equal to ntrials (", ntrials, "), ",
      "but has length ", length(reward[[d]]), "."
    ))
  }
  
  # Build c_vec array
  c_vec <- array(0, dim = c(D, C, L, A, ntrials))
  
  for (t in 1:ntrials) {
    arms <- strsplit(trial_patterns[t], "_")[[1]]
    for (a in 1:A) {
      # leave c_vec of nil arms as all zeros
      if (arms[a] == "0") next  
      seg_letters <- strsplit(arms[a], "")[[1]]
      # A feature is active on this arm if all its letters are present in the segment
      cue_idx <- which(sapply(base_cues, function(f) {
        f_letters <- strsplit(f, "")[[1]]
        all(f_letters %in% seg_letters)
      }))
      for (d in 1:D) {
        for (l in 1:L) {
          c_vec[d, cue_idx, l, a, t] <- 1
        }
      }
    }
  }
  
  # Build win array: D x A x ntrials
  win <- array(0, dim = c(D, A, ntrials))
  for (d in 1:D) {
    for (t in 1:ntrials) {
      win[d, , t] <- reward[[d]][t]
    }
  }
  
  return(list(
    c_vec = c_vec,
    win = win,
    type_seq  = list(trial_patterns) #wrapped in list to match TaskInit format
  ))
}
