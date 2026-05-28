# Example of phase_def usage. Any cue or phase that you do not need can
# be omitted, no need to manually set it to 0.
# trial type names use underscore to separate arms. E.g. "A_B" means arm1 gets 
# cue A, arm2 gets cue B. "AB_C" means arm1 gets cues A+B, arm2 gets cue C.
# "0" is a reserved keyword meaning no cues on that arm (all-zero cue vector). 
# E.g. "B_0" = forced choice, only arm1 is present.
# Single-arm names (no underscore) remain valid when A=1.
phase_def <- list(
  
  phase1 = list(
    trials = list(
      A_B   = 20,
      B_A   = 10,
      AB_C  = 12,
      C_AB  = 6,
      B_0   = 4
    ),
    rewards = list(
      D1 = c(A = 0.80, B = 0.20, AB = 0.50)
    )
  ),
  
  phase2 = list(
    trials = list(
      A_B  = 16,
      B_A  = 16,
      AB_C = 8,
      C_AB = 8
    ),
    rewards = list(
      D1 = c(A = 0.10, B = 0.50, C = 0.00, AB = 0.25)
    )
  )
)


TaskInit <- function(
    phase_def,
    pseudo = 10,
    par_L = 1,
    center = 0.5,
    ncop = 10,
    r0 = 1,
    contextShift = NULL,
    block = NULL,
    probabilistic = FALSE,
    features = NULL
) {
  
  # ============================================================
  # 1. Extract experiment details
  # ============================================================
  
  # infer number of arms per trial from the phase_def
  all_trial_names_raw <- unique(unlist(lapply(phase_def, function(ph) names(ph$trials))))
  arm_counts <- sapply(all_trial_names_raw, function(nm) length(strsplit(nm, "_")[[1]]))
  if (length(unique(arm_counts)) != 1) stop(
    "Inconsistent number of arms across trial types: every trial type name must have ",
    "the same number of underscore-separated segments."
  )
  A <- unique(arm_counts)
  
  # Helper to split a trial type name into per-arm cue strings.
  parse_trial_name <- function(nm) strsplit(nm, "_")[[1]]
  
  # error pops if user accidentally specifies reward probability
  # for a cue type that doesn't exist as an arm segment in the given phase.
  for (phase_name in names(phase_def)) {
    phase <- phase_def[[phase_name]]
    # collect every distinct non-nil arm segment in this phase's trial names
    phase_arm_segs <- unique(Filter(function(s) s != "0",
                                    unlist(lapply(names(phase$trials), parse_trial_name))
    ))
    stopifnot("all reward cues must appear in trials" =
                all(unlist(lapply(phase$rewards, names)) %in% phase_arm_segs))  # [CHANGED-BUG1]
  }
  
  # Validate that reward probability × trial count is a whole number for
  # every phase, dimension, and cue, so no rounding occurs during sampling.
  for (phase_name in names(phase_def)) {
    phase <- phase_def[[phase_name]]
    for (d_name in names(phase$rewards)) {
      for (base_cue in names(phase$rewards[[d_name]])) {
        # count trials where base_cue appears as a complete arm segment
        n <- sum(sapply(names(phase$trials), function(nm) { 
          arm_cues <- parse_trial_name(nm)
          appearances <- sum(arm_cues == base_cue)
          appearances * phase$trials[[nm]]
        }))
        p <- phase$rewards[[d_name]][[base_cue]]
        if (!probabilistic && round(n * p, 10) %% 1 != 0) stop(paste0(
          phase_name, ", ", d_name, ", cue '", base_cue, "': ",
          n, " trials × ", p, " reward prob = ", n * p, " (not a whole number)."
        ))
      }
    }
  }
  
  # Validate and expand the block argument into a per-phase vector.
  # If block is a single number, replicate it for every phase.
  # If block is a vector, it must have exactly one entry per phase.
  if (!is.null(block)) {
    if (length(block) == 1) {
      block <- rep(block, length(phase_def))
    } else if (length(block) != length(phase_def)) {
      stop(paste0(
        "block must be either a single number or a vector with one entry per phase ",
        "(expected length ", length(phase_def), ", got length ", length(block), ")."
      ))
    }
    
    # For each phase, verify that the block length evenly divides every
    # trial count AND that each block can contain an identical set of trial types
    for (phase_i in seq_along(phase_def)) {
      phase_trials <- phase_def[[phase_i]]$trials
      block_len    <- block[phase_i]
      nphase_i     <- sum(unlist(phase_trials))
      
      # Total trials in the phase must be divisible by the block length
      if (nphase_i %% block_len != 0) {
        stop(paste0(
          "Phase ", phase_i, ": total trials (", nphase_i, ") is not divisible by ",
          "block length (", block_len, "). All blocks must be identical."
        ))  
      }
      
      n_blocks <- nphase_i / block_len  # [NEW] number of blocks in this phase
      
      # Each cue's trial count must be divisible by n_blocks so that
      # every block receives the same number of that cue type
      for (cue in names(phase_trials)) {
        if (phase_trials[[cue]] %% n_blocks != 0) {
          stop(paste0(
            "Phase ", phase_i, ", cue '", cue, "': trial count (", phase_trials[[cue]],
            ") is not divisible by the number of blocks (", n_blocks, "). ",
            "It is impossible to create identical blocks with this combination of ",
            "trial counts and block length. Please adjust your trial counts or block length."
          ))  
        }
      }
    }
  }
  
  #extract number of reward dimensions from phase_def
  D <- max(sapply(phase_def, function(phase) length(phase$rewards)))
  
  
  # ============================================================
  # 2. Cue types and which are present
  # ============================================================

  # Collect all distinct non-nil arm segment strings from all
  # trial type names across all phases (e.g., "A", "B", "AB").
  all_cue_names <- unique(Filter(function(s) s != "0",
                                 unlist(lapply(phase_def, function(phase) {
                                   unlist(lapply(names(phase$trials), parse_trial_name))
                                 }))
  ))
  
  # Valid atomic letters that actually appear in at least one arm segment
  valid_atoms <- sort(unique(unlist(strsplit(all_cue_names, ""))))
  
  # Valid multi-letter segments that appear as full arm segments
  valid_segments <- all_cue_names
  
  # -----------------------------
  # 2a. If features is user-specified
  # -----------------------------
  if (!is.null(features)) {
    
    # Basic type and length checks
    if (!is.character(features)) {
      stop("features must be NULL or a character vector, but is of type ",
           typeof(features), ".")
    }
    if (length(features) == 0L) {
      stop("features must be NULL or a non-empty character vector.")
    }
    
    # Remove duplicates, warn if any
    if (anyDuplicated(features)) {
      warning("features contains duplicate entries; duplicates will be ignored.")
      features <- unique(features)
    }
    
    # Validate each feature:
    # - all letters must be among valid_atoms
    # - the feature string must appear as a contiguous substring of at least
    #   one arm segment in all_cue_names
    invalid_features <- character(0)
    
    for (f in features) {
      letters_f <- strsplit(f, "")[[1]]
      
      if (anyDuplicated(letters_f)) {
        stop(
          "features contains '", f, "', which has duplicate letters. ",
          "Each feature must contain every letter at most once."
        )
      }
      
      # letters must exist somewhere in the experiment
      if (!all(letters_f %in% valid_atoms)) {
        invalid_features <- c(invalid_features, f)
        next
      }
      
      # f must appear as a substring in at least one cue segment
      appears <- any(vapply(all_cue_names, function(seg) {
        seg_letters <- strsplit(seg, "")[[1]]
        f_letters   <- strsplit(f,   "")[[1]]
        all(f_letters %in% seg_letters)
      }, logical(1)))
      
      if (!appears) {
        invalid_features <- c(invalid_features, f)
      }
    }
    
    if (length(invalid_features) > 0L) {
      stop(
        "features contains cue types that never appear in the experiment as substrings: ",
        paste(invalid_features, collapse = ", "),
        ".\n\nValid letters are: ",
        if (length(valid_atoms) > 0) paste(valid_atoms, collapse = ", ") else "<none>",
        ".\nValid substrings must occur within at least one of these arm segments: ",
        if (length(valid_segments) > 0) paste(valid_segments, collapse = ", ") else "<none>",
        "."
      )
    }
    
    # base_cues: exactly the features the user wants to value
    base_cues <- features
    
    # Build cue_map. For each arm-segment string (e.g., "AB", "ABC"):
    # - single-letter features active if their letter is in that segment
    # - multi-letter features active if their string is a substring of that segment
    cue_map <- setNames(vector("list", length(all_cue_names)), all_cue_names)
    
    for (cue in all_cue_names) {
      letters_cue <- strsplit(cue, "")[[1]]
      
      idx <- integer(0)
      
      for (i in seq_along(base_cues)) {
        f         <- base_cues[i]
        f_letters <- strsplit(f, "")[[1]]
        if (all(f_letters %in% letters_cue)) {
          idx <- c(idx, i)
        }
      }
      
      idx <- unique(idx)
      
      if (length(idx) == 0L) {
        stop(
          "features does not include any feature that can represent cue '", cue, "'. ",
          "For example, a trial contains arm segment '", cue,
          "', but none of the specified features (single letters or substrings) ",
          "match it. Please adjust features so every arm segment activates at least one feature."
        )
      }
      
      cue_map[[cue]] <- idx
    }
    
    # -----------------------------
    # 2b. Default: hybrid atom+config encoding (Cochran–Cisler style)
    # -----------------------------
  } else {
    
    # Elemental cues = individual letters (A, B, C, ...)
    atom_cues <- valid_atoms
    
    # Configural cues = multi-letter arm segments (AB, AD, etc.) as full segments
    config_cues <- sort(setdiff(all_cue_names, atom_cues))
    
    # base_cues = all cues (elemental + configural)
    base_cues <- c(atom_cues, config_cues)
    
    # Build cue_map:
    # - For a single-letter cue "A": active features = ["A"]
    # - For a configural cue "AB": active features = ["A","B","AB"]
    # (config_cues only fire when the entire segment equals them)
    cue_map <- setNames(
      lapply(all_cue_names, function(cue) {
        letters   <- strsplit(cue, "")[[1]]
        idx_atoms <- match(letters, base_cues, nomatch = NA_integer_)
        idx_atoms <- idx_atoms[!is.na(idx_atoms)]
        idx_conf  <- which(base_cues == cue)
        unique(c(idx_atoms, idx_conf))
      }),
      all_cue_names
    )
    
    # every cue must have at least one active feature
    empty_cues <- names(cue_map)[vapply(cue_map, length, integer(1)) == 0L]
    if (length(empty_cues) > 0L) {
      stop(
        "Internal error: some cues have no active features under the default encoding: ",
        paste(empty_cues, collapse = ", "), "."
      )
    }
  }
  
  # Number of feature dimensions
  C <- length(base_cues)
  
  # ============================================================
  # 3. Trial counts per phase (from trials, not from rewards)
  # ============================================================

  # vector automatically sizes itself to however many phases are provided
  nphase <- integer(length(phase_def))
  for (i in seq_along(phase_def)) {
    nphase[i] <- sum(unlist(phase_def[[i]]$trials))
  }
  
  ntrials <- sum(nphase)
  
  # ============================================================
  # 4. Params and world allocation
  # ============================================================

  par <- list()
  par$A       <- A
  par$D       <- D
  par$C       <- C
  par$L       <- par_L
  par$center  <- center
  par$r0      <- r0
  par$context <- contextShift
  par$ntrials = ntrials
  par$base_cues <- base_cues
  par$ls$ncop   <- ncop
  
  L <- par$L * par$ls$ncop
  
  world <- list()
  world$c_vec <- array(0, dim = c(par$D, par$C, L, par$A, ntrials))
  world$win   <- array(0, dim = c(par$D, par$A, ntrials))
  
  
  # ============================================================
  # 5. Randomization helper
  # ============================================================

  randomize <- function(types, pseudo_local = pseudo) {
    repeat {
      seq <- sample(types)
      if (max(rle(seq)$lengths) <= pseudo_local) return(seq)
    }
  }
  
  # Block-based randomization helper.
  # Divides the full trial list into identical blocks (each block contains
  # the same proportional mix of cue types), shuffles each block independently,
  # and concatenates them. The pseudo argument is ignored when this is used.
  randomize_blocks <- function(types, block_len) {
    n_blocks  <- length(types) / block_len  
    block_template <- rep(names(table(types)), times = table(types) / n_blocks)
    # Shuffle each block independently and concatenate
    result <- unlist(lapply(seq_len(n_blocks), function(i) sample(block_template))) 
    return(result)  # [NEW]
  }
  
  
  # ============================================================
  # 6. Phase loop 
  # ============================================================

  type_seq_all <- vector("list", length(phase_def))  # [NEW]
  
  # loop over all phases present in phase_def
  phase_end <- 0
  
  for (phase_i in seq_along(phase_def)) {  # [NEW]
    
    phase      <- phase_def[[phase_i]]     # [NEW]
    nphase_i   <- nphase[phase_i]          # [NEW]
    
    # Create trial-type sequence
    type_seq <- unlist(lapply(names(phase$trials), function(cue) {  # [NEW]
      rep(cue, phase$trials[[cue]])
    }), use.names = FALSE)
    
    # Choose randomization method based on whether block is specified.
    # When block is non-NULL, use block-based shuffling and ignore pseudo.
    # When block is NULL, use the original pseudo-randomization.
    # skip randomization when there is only one trial type as it would loop forever
    if (length(phase$trials) > 1) {
      if (!is.null(block)) {
        type_seq <- randomize_blocks(type_seq, block_len = block[phase_i])
      } else {
        type_seq <- randomize(type_seq)
      }
    }
    
    # Store randomization after it is finalized
    type_seq_all[[phase_i]] <- type_seq
    
    # Parse each trial type in type_seq into a list of per-arm
    # cue pattern vectors (C-length binary). arm_cue_patterns[[a]] is a C x
    # nphase_i matrix of 0/1 cue indicators for arm a.
    # "0" arm segments produce an all-zero column (no active cues).
    arm_cue_patterns <- lapply(1:A, function(a) {
      mat <- sapply(type_seq, function(nm) {
        arm_str <- parse_trial_name(nm)[a]
        if (arm_str == "0") return(rep(0L, C))  
        as.integer(1:C %in% cue_map[[arm_str]])
      })
      if (is.null(dim(mat)) || length(dim(mat)) < 2) {
        dim(mat) <- c(C, length(type_seq))
      }
      mat
    })
    
    # Locations of each trial type in the trial order
    type_idx <- lapply(names(phase$trials), function(cue) {
      which(type_seq == cue)
    })
    names(type_idx) <- names(phase$trials) 
    
    # Assign rewards per arm per trial. Reward probability is
    # looked up by exact arm segment string (cue type), not by base letter.
    # phase_win: D x A x nphase_i — reward outcome per arm per trial.
    phase_win <- array(0, dim = c(D, A, nphase_i))  
    
    for (d in 1:length(phase$rewards)) {
      probs <- phase$rewards[[d]]
      
      for (base_cue in names(probs)) {
        p <- probs[[base_cue]]
        for (a in 1:A) {
          # find trials where arm a has exactly base_cue as its arm segment
          idx <- which(sapply(type_seq, function(nm) {
            parse_trial_name(nm)[a] == base_cue  
          }))
          n_occ <- length(idx)
          if (n_occ > 0) {
            if (probabilistic) {
              phase_win[d, a, idx] <- rbinom(n_occ, 1, p)
            } else {
              lucky <- sample(idx, size = round(p * n_occ), replace = FALSE)
              phase_win[d, a, lucky] <- 1           
            }
          }
        }
      }
    }
    
    # Write into world
    start_i <- phase_end + 1      
    end_i   <- phase_end + nphase_i 
    
    # Write cue vectors per arm using the parsed arm cue patterns.
    # Nil arms naturally get all-zero columns from arm_cue_patterns.
    for (l in 1:L) {
      for (d in 1:D) {
        for (a in 1:A) {                                                     
          world$c_vec[d, , l, a, start_i:end_i] <- arm_cue_patterns[[a]]  
        }
      }
    }
    
    # Write reward vector per arm and dimension.
    for (d in 1:D) {
      for (a in 1:A) {                                                        
        world$win[d, a, start_i:end_i] <- phase_win[d, a, ]                
      }
    }
    
    phase_end <- end_i 
  }
  
  return(list(
    world = world,
    par = par,
    type_seq = type_seq_all,
    cue_names = base_cues
  ))
}
