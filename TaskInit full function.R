# Example of phase_def usage. Any cue or phase that you do not need can
# be omitted, no need to manually set it to 0.
# [CHANGED-MULTI] trial type names now use underscore to separate arms.
# "A_B" means arm1 gets cue A, arm2 gets cue B.
# "AB_C" means arm1 gets cues A+B, arm2 gets cue C.
# [CHANGED-NIL] "0" is a reserved keyword meaning no cues on that arm
# (all-zero cue vector). E.g. "B_0" = forced choice, only arm1 is present.
# Single-arm names (no underscore) remain valid when A=1.
phase_def <- list(
  
  phase1 = list(
    trials = list(
      A_B   = 20,
      B_A   = 10,
      AB_C  = 12,
      C_AB  = 6,
      B_0   = 4   # [CHANGED-NIL] forced-choice trial: only arm1 is present
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
    ncop = 10,
    center = 0.5,
    r0 = 1,
    contextShift = NULL,
    block = NULL,  # [NEW] optional block length: single number (applies to all phases)
    # or vector (one entry per phase). When non-NULL, pseudo is ignored.
    probabilistic = FALSE
) {
  
  # [CHANGED] removed the stopifnot phase cap entirely; any number of phases
  # is now supported
  
  # [CHANGED-MULTI] Infer A (number of arms) from the trial type names.
  # Split every trial type name on "_" and count the segments; all must agree.
  all_trial_names_raw <- unique(unlist(lapply(phase_def, function(ph) names(ph$trials))))
  arm_counts <- sapply(all_trial_names_raw, function(nm) length(strsplit(nm, "_")[[1]]))
  if (length(unique(arm_counts)) != 1) stop(
    "Inconsistent number of arms across trial types: every trial type name must have ",
    "the same number of underscore-separated segments."
  )
  A <- unique(arm_counts)  # [CHANGED-MULTI] A inferred, no longer a parameter
  
  # [CHANGED-MULTI] Helper to split a trial type name into per-arm cue strings.
  parse_trial_name <- function(nm) strsplit(nm, "_")[[1]]
  
  # [CHANGED-MULTI] error pops if user accidentally specifies reward probability
  # for a cue type that doesn't exist as an arm segment in the given phase.
  # [CHANGED-BUG1] now checks for exact arm-segment match rather than letter
  # membership, so "AB" in rewards requires "AB" to appear as a full arm
  # segment in trials, not just "A" and "B" as individual letters.
  for (phase_name in names(phase_def)) {
    phase <- phase_def[[phase_name]]
    # collect every distinct non-nil arm segment in this phase's trial names
    phase_arm_segs <- unique(Filter(function(s) s != "0",  # [CHANGED-NIL]
                                    unlist(lapply(names(phase$trials), parse_trial_name))
    ))  # [CHANGED-BUG1]
    stopifnot("all reward cues must appear in trials" =
                all(unlist(lapply(phase$rewards, names)) %in% phase_arm_segs))  # [CHANGED-BUG1]
  }
  
  # [NEW] Validate that reward probability × trial count is a whole number for
  # every phase, dimension, and cue, so no rounding occurs during sampling.
  # [CHANGED-MULTI] trial count for a cue type = sum of trial counts for all
  # trial types that contain that exact cue type as a complete arm segment.
  # [CHANGED-BUG2] now uses exact arm-segment matching (arm_cues == base_cue)
  # rather than letter-membership, so "AB" counts only trials where arm has
  # exactly "AB", not trials where arm has "A" or "B" separately.
  for (phase_name in names(phase_def)) {
    phase <- phase_def[[phase_name]]
    for (d_name in names(phase$rewards)) {
      for (base_cue in names(phase$rewards[[d_name]])) {
        # count trials where base_cue appears as a complete arm segment
        n <- sum(sapply(names(phase$trials), function(nm) {  # [CHANGED-BUG2]
          arm_cues <- parse_trial_name(nm)
          appearances <- sum(arm_cues == base_cue)  # [CHANGED-BUG2] exact match
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
  
  # [NEW] Validate and expand the block argument into a per-phase vector.
  # If block is a single number, replicate it for every phase.
  # If block is a vector, it must have exactly one entry per phase.
  if (!is.null(block)) {
    if (length(block) == 1) {
      block <- rep(block, length(phase_def))  # [NEW] expand scalar to per-phase vector
    } else if (length(block) != length(phase_def)) {
      stop(paste0(
        "block must be either a single number or a vector with one entry per phase ",
        "(expected length ", length(phase_def), ", got length ", length(block), ")."
      ))  # [NEW] length mismatch error
    }
    
    # [NEW] For each phase, verify that the block length evenly divides every
    # trial count AND that each block can contain an identical set of trial types.
    # Identical blocks require: for each cue, trials[[cue]] / n_blocks is an integer,
    # where n_blocks = nphase_i / block_len.
    for (phase_i in seq_along(phase_def)) {
      phase_trials <- phase_def[[phase_i]]$trials
      block_len    <- block[phase_i]
      nphase_i     <- sum(unlist(phase_trials))
      
      # [NEW] Total trials in the phase must be divisible by the block length
      if (nphase_i %% block_len != 0) {
        stop(paste0(
          "Phase ", phase_i, ": total trials (", nphase_i, ") is not divisible by ",
          "block length (", block_len, "). All blocks must be identical."
        ))  # [NEW]
      }
      
      n_blocks <- nphase_i / block_len  # [NEW] number of blocks in this phase
      
      # [NEW] Each cue's trial count must be divisible by n_blocks so that
      # every block receives the same number of that cue type
      for (cue in names(phase_trials)) {
        if (phase_trials[[cue]] %% n_blocks != 0) {
          stop(paste0(
            "Phase ", phase_i, ", cue '", cue, "': trial count (", phase_trials[[cue]],
            ") is not divisible by the number of blocks (", n_blocks, "). ",
            "It is impossible to create identical blocks with this combination of ",
            "trial counts and block length. Please adjust your trial counts or block length."
          ))  # [NEW]
        }
      }
    }
  }
  
  # ============ 1. Extract D (number of reward dimensions) ============
  
  D <- max(sapply(phase_def, function(phase) length(phase$rewards)))
  
  # ============ 2. Cue types and which are present ============
  
  # [CHANGED-MULTI] Collect all distinct non-nil arm segment strings from all
  # trial type names across all phases.
  # [CHANGED-NIL] "0" segments are excluded.
  all_cue_names <- unique(Filter(function(s) s != "0",  # [CHANGED-NIL]
                                 unlist(lapply(phase_def, function(phase) {
                                   unlist(lapply(names(phase$trials), parse_trial_name))
                                 }))
  ))  # [CHANGED-MULTI] these are per-arm cue strings, e.g. "AB", "C"
  
  # Derive the set of unique base cues by splitting each name into individual
  # letters (e.g. "ABD" → c("A","B","D")) and taking the union
  base_cues <- sort(unique(unlist(strsplit(all_cue_names, ""))))
  
  # Build cue_map automatically: each per-arm cue string maps to the integer
  # positions of its constituent letters within base_cues
  cue_map <- setNames(
    lapply(all_cue_names, function(cue) {
      which(base_cues %in% strsplit(cue, "")[[1]])
    }),
    all_cue_names
  )
  
  # C is now the number of unique base cues found in this experiment
  C <- length(base_cues)
  
  # ============ 3. Trial counts per phase (from trials, not from rewards) ============
  
  # [CHANGED] replaced integer(6) with integer(length(phase_def)) so the
  # vector automatically sizes itself to however many phases are provided
  nphase <- integer(length(phase_def))
  for (i in seq_along(phase_def)) {
    nphase[i] <- sum(unlist(phase_def[[i]]$trials))
  }
  
  ntrials <- sum(nphase)
  
  # ============ 4. Params and world allocation ============
  
  par <- list()
  par$A       <- A
  par$D       <- D
  par$C       <- C
  par$L       <- par_L
  par$center  <- center
  par$r0      <- r0
  par$context <- contextShift
  par$ls$ncop <- ncop
  par$ntrials = ntrials
  
  L <- par$L * par$ls$ncop
  
  world <- list()
  world$c_vec <- array(0, dim = c(par$D, par$C, L, par$A, ntrials))
  world$win   <- array(0, dim = c(par$D, par$A, ntrials))
  
  
  # ============ 5. Randomization helper ============
  
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
    n_blocks  <- length(types) / block_len          # [NEW] number of complete blocks
    # Build the single-block template: the proportional share of each cue type
    block_template <- rep(names(table(types)), times = table(types) / n_blocks)  # [NEW]
    # Shuffle each block independently and concatenate
    result <- unlist(lapply(seq_len(n_blocks), function(i) sample(block_template)))  # [NEW]
    return(result)  # [NEW]
  }
  
  
  # ============ 6. Phase loop ============
  
  type_seq_all <- vector("list", length(phase_def))  # [NEW]
  
  # [NEW] replaces the three hard-coded phase blocks (6, 7, 8) with a single
  # loop over all phases present in phase_def
  phase_end <- 0  # [NEW] rolling end-index; replaces end1/end2 hand-offs
  
  for (phase_i in seq_along(phase_def)) {  # [NEW]
    
    phase      <- phase_def[[phase_i]]     # [NEW]
    nphase_i   <- nphase[phase_i]          # [NEW]
    
    # Create trial-type sequence
    type_seq <- unlist(lapply(names(phase$trials), function(cue) {  # [NEW]
      rep(cue, phase$trials[[cue]])
    }), use.names = FALSE)
    
    # [NEW] Choose randomization method based on whether block is specified.
    # When block is non-NULL, use block-based shuffling and ignore pseudo.
    # When block is NULL, use the original pseudo-randomization.
    # [CHANGED] skip randomization when there is only one trial type — shuffling
    # is meaningless and pseudo-randomization would loop forever
    if (length(phase$trials) > 1) {
      if (!is.null(block)) {
        type_seq <- randomize_blocks(type_seq, block_len = block[phase_i])
      } else {
        type_seq <- randomize(type_seq)
      }
    }
    
    type_seq_all[[phase_i]] <- type_seq  # [NEW] store after randomization is finalized
    
    # [CHANGED-MULTI] Parse each trial type in type_seq into a list of per-arm
    # cue pattern vectors (C-length binary). arm_cue_patterns[[a]] is a C x
    # nphase_i matrix of 0/1 cue indicators for arm a.
    # [CHANGED-NIL] "0" arm segments produce an all-zero column (no active cues).
    arm_cue_patterns <- lapply(1:A, function(a) {  # [CHANGED-MULTI]
      mat <- sapply(type_seq, function(nm) {
        arm_str <- parse_trial_name(nm)[a]
        if (arm_str == "0") return(rep(0L, C))  # [CHANGED-NIL] all-zero for nil
        as.integer(1:C %in% cue_map[[arm_str]])
      })
      if (is.null(dim(mat)) || length(dim(mat)) < 2) {
        dim(mat) <- c(C, length(type_seq))
      }
      mat
    })
    
    # Locations of each trial type in the trial order
    type_idx <- lapply(names(phase$trials), function(cue) {  # [NEW]
      which(type_seq == cue)
    })
    names(type_idx) <- names(phase$trials)  # [NEW]
    
    # [CHANGED-MULTI] Assign rewards per arm per trial. Reward probability is
    # looked up by exact arm segment string (cue type), not by base letter.
    # phase_win: D x A x nphase_i — reward outcome per arm per trial.
    phase_win <- array(0, dim = c(D, A, nphase_i))  # [CHANGED-MULTI]
    
    for (d in 1:length(phase$rewards)) {
      probs <- phase$rewards[[d]]
      
      for (base_cue in names(probs)) {
        p <- probs[[base_cue]]
        for (a in 1:A) {
          # find trials where arm a has exactly base_cue as its arm segment
          idx <- which(sapply(type_seq, function(nm) {
            parse_trial_name(nm)[a] == base_cue  # [CHANGED-BUG2] exact match
          }))
          n_occ <- length(idx)
          if (n_occ > 0) {
            if (probabilistic) {
              phase_win[d, a, idx] <- rbinom(n_occ, 1, p)  # [CHANGED-MULTI]
            } else {
              lucky <- sample(idx, size = round(p * n_occ), replace = FALSE)
              phase_win[d, a, lucky] <- 1                   # [CHANGED-MULTI]
            }
          }
        }
      }
    }
    
    # Write into world
    start_i <- phase_end + 1         # [NEW]
    end_i   <- phase_end + nphase_i  # [NEW]
    
    # [CHANGED-MULTI] Write cue vectors per arm using the parsed arm cue patterns.
    # [CHANGED-NIL] nil arms naturally get all-zero columns from arm_cue_patterns.
    for (l in 1:L) {
      for (d in 1:D) {
        for (a in 1:A) {                                                        # [CHANGED-MULTI]
          world$c_vec[d, , l, a, start_i:end_i] <- arm_cue_patterns[[a]]       # [CHANGED-MULTI]
        }
      }
    }
    
    # [CHANGED-MULTI] Write reward vector per arm and dimension.
    for (d in 1:D) {
      for (a in 1:A) {                                                          # [CHANGED-MULTI]
        world$win[d, a, start_i:end_i] <- phase_win[d, a, ]                    # [CHANGED-MULTI]
      }
    }
    
    phase_end <- end_i  # [NEW] advance rolling index for the next phase
  }
  
  return(list(
    world = world,
    par = par,
    type_seq = type_seq_all
  ))
}

out = TaskInit(phase_def)
world <- out$world
par   <- out$par