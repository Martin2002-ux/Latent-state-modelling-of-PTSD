# =========================================================================
# LatentState: Latent-state learning approach
# Inputs:
#   par   = list of parameters
#   agent = list of variables tracked by agent (NULL for initialization)
#   world = list of variables generated from task/world
# Output:
#   agent = updated list of variables tracked by agent
# =========================================================================

LatentState <- function(par, agent, world) {
  
  # Total number of latent states (number of copies x initial # hypotheses)
  L <- par$L * par$ls$ncop
  
  # -----------------------------------------------------------------------
  # Initialization of agent
  # -----------------------------------------------------------------------
  if (is.null(agent)) {
    
    # Special case: all variables already initialized externally, initialV is a
    #vector of initial variables
    if (!is.null(par$initialV)) {
      agent <- par$initialV
      
    } else {
      agent <- list()
      
      # Initial expectations (neutral, unless otherwise specified), allows to
      #set initial V directly without setting the other initial variables
      if (!is.null(par$initialex)) {
        agent$V <- par$initialex
      } else {
        #set all initial values to 0 for no prior expectations of value.
        #C is cues, D is reward dimensions, representing how many reward 
        #outcomes an agent needs to track, and is typically 1 (e.g. tracking shock)
        agent$V <- array(0, dim = c(par$D, par$C, L))
      }
      
      # Each hypothesis (contained in par$L) is initially active, and no
      #other states are active, meaning we have 1 of each initial state
      agent$Lmax <- 1
      
      # Change point statistic
      agent$q <- 0
      
      # Validate pre-specified initial latent-state beliefs (must be a probability vector of length par$L)
      if (!is.null(par$lsb0)) {
        if (length(par$lsb0) != par$L) stop(paste0(
          "par$lsb0 must have length par$L (", par$L, "), ",
          "but has length ", length(par$lsb0), "."
        ))
        if (round(sum(par$lsb0), 10) != 1) stop(paste0(
          "par$lsb0 must sum to 1, but sums to ", round(sum(par$lsb0), 8), "."
        ))
      }
      
      # Initial log-beliefs (evenly divided among initial hypotheses)
      agent$lsb <- rep(0, L)
      if (!is.null(par$lsb0)) {
        agent$lsb[1:par$L] <- par$lsb0
      } else {
        agent$lsb[1:par$L] <- 1 / par$L
      }
      
      # Initial B matrix (identity) for cues per dimension & state.
      #Uses par$C twice to create a square matrix. There is a separate
      #square matrix for each latent state L and each reward dimension D
      agent$B <- array(0, dim = c(par$C, par$C, par$D, L))
      #set everything to an identity matrix
      for (d in 1:par$D) {
        for (l in 1:L) {
          agent$B[, , d, l] <- diag(par$C)
        }
      }
      
      # Initial variance of expected mean
      agent$sigmasq <- par$ls$sigma0^2
    }
    
    # -----------------------------------------------------------------------
    # Update of agent
    # -----------------------------------------------------------------------
  } else {
    
    # --------------------------------
    # Get variables
    # --------------------------------
    V        <- agent$V       # Associative strengths:  D x C x L
    lsb      <- agent$lsb     # Latent state beliefs:   length L
    B        <- agent$B       # Effort matrix:          C x C x D x L
    Lmax     <- agent$Lmax    # No. of active copies
    q        <- agent$q       # Change point statistic
    sigmasq0 <- agent$sigmasq # Variance
    win      <- world$win - par$center  # Outcome: D x A matrix
    c_vec    <- world$c_vec             # Cue vector: D x C x L x A array
    
    # --------------------------------
    # Context or temporal shift
    # --------------------------------
    lsb <- ContextShift(lsb, Lmax, par)
    
    # --------------------------------
    # Choice behavior
    # --------------------------------
    
    # Value function per arm and latent state: D x L x A, where each arm
    #represents a choice, which can have its own cue vector
    # we sum over dim(2) so that we get overall value for each D x L
    mu_arr <- array(0, dim = c(par$D, L, par$A))
    for (a in 1:par$A) {
      c_slice <- c_vec[, , , a]
      if (length(dim(c_slice)) < 3) {
        dim(c_slice) <- c(par$D, par$C, L)
      }
      mu_arr[, , a] <- apply(V * c_slice, c(1, 3), sum)
    }
    
    #matrix multiply V of each arm by lsb so you get a belief-weighted expected value
    mu0 <- matrix(0, nrow = par$D, ncol = par$A)
    for (a in 1:par$A) {
      mu0[, a] <- mu_arr[, , a] %*% lsb 
    }
    
    # Choose arm and update relevant variables based on chosen arm
    arm   <- SoftMaxChoice(mu0, par)
    win   <- win[, arm]          # D-length vector
    mu    <- mu_arr[, , arm]     # D x L matrix
    c_vec <- c_vec[, , , arm]
    if (length(dim(c_vec)) < 3) {
      dim(c_vec) <- c(par$D, par$C, L)
    }
    
    # Guard against dimension-dropping when par$A == 1
    #when mu only has one dimension, it becomes a vector, which causes issues later
    
    if (is.null(dim(mu)) || length(dim(mu)) < 2) {
      dim(mu) <- c(par$D, L)
    }
    #same idea but applied to c_vec which needs 3 dimensions, and R may simplify 
    #it to less
    if (is.null(dim(c_vec)) || length(dim(c_vec)) < 3) {
      dim(c_vec) <- c(par$D, par$C, L)
    }
    
    # --------------------------------
    # Updates
    # --------------------------------
    
    # Update likelihood
    LL <- LogLikelihood(win, mu, sigmasq0, Lmax, L, par$ls$ncop)
    
    # Approximate Bayesian filter of latent-state beliefs
    lsb <- LatentStateBeliefs(lsb, Lmax, LL, par)
    
    # Add new state if warranted
    ns   <- NewState(V, lsb, Lmax, c_vec, win, LL, q, par)
    V    <- ns$V
    lsb  <- ns$lsb
    Lmax <- ns$Lmax
    q    <- ns$q
    
    # Update associative strengths, effort matrix, variance
    upd     <- Update(V, B, sigmasq0, lsb, Lmax, c_vec, win, par)
    V       <- upd$V
    B       <- upd$B
    sigmasq <- upd$sigmasq
    
    # --------------------------------
    # Output updated variables
    # --------------------------------
    agent$V       <- V
    agent$lsb     <- lsb
    agent$B       <- B
    agent$sigmasq <- sigmasq
    agent$Lmax    <- Lmax
    agent$q       <- q
    agent$arm     <- arm
    agent$win     <- win
    agent$mu      <- mu0
    agent$mu_arr  <- mu_arr   # D x L x A, unweighted per-state expected value
  }
  
  return(agent)
}


# =========================================================================
# ContextShift: Apply temporal or spatial context shifts to beliefs
# =========================================================================
ContextShift <- function(lsb, Lmax, par) {
  
  # Temporal shift: beliefs decay towards uniform over time gap
  if (par$t != 1 && !is.null(par$time)) {
    prob_stay <- (1 - par$ls$gamma)^(par$time[par$t] - par$time[par$t - 1] - 1)
    ix <- 1:(Lmax * par$L)
    lsb[ix] <- lsb[ix] * prob_stay + (1 - prob_stay) / (Lmax * par$L)
  }
  
  # Spatial/visual context shift: re-initialize beliefs to uniform
  #This resets beliefs about which latent state is active. However, it does
  #not erase the latent states that have already been created. In sum, lsb
  #is reset while V, B, Lmax, sigmasq, and q are all retained
  if (!is.null(par$context) && any(par$context == par$t)) {
    ix <- 1:(par$L * Lmax)
    lsb[ix] <- 1
    lsb     <- lsb / sum(lsb)
  }
  
  return(lsb)
}


# =========================================================================
# LogLikelihood: Log-likelihood of each latent state
# Inputs:
#   win     = outcome vector, length D
#   mu      = predicted values, D x L matrix
#   sigmasq = current variance estimate (scalar)
#   Lmax    = number of active copies
#   L       = total number of latent states
#   ncop    = number of copies (par$ls$ncop)
# =========================================================================
LogLikelihood <- function(win, mu, sigmasq, Lmax, L, ncop) {
  
  LL <- rep(0, L)
  
  # Set inactive (not yet activated) states to -Inf
  # Active states: indices 1 to (L/ncop)*Lmax
  start_inactive <- (L / ncop) * (Lmax + 1) + 1
  if (start_inactive <= L) { #the 'if' serves as a guard in case Lmax=ncop, which does not produce an empty start_inactive even though we want it to
    LL[start_inactive:L] <- -Inf
  }
  
  # Index of activated states
  ix <- 1:((L / ncop) * Lmax)
  
  # Log-likelihood: sum over reward dimensions (dim 1)
  # this computes the log likelihood of each state given the discrepancy between
  #expected and observed value
  LL[ix] <- colSums(-(mu[, ix, drop = FALSE] - win)^2 / (2 * sigmasq))
  
  return(LL)
}


# =========================================================================
# LatentStateBeliefs: Approximate Bayesian filter update of latent beliefs
# =========================================================================
LatentStateBeliefs <- function(lsb, Lmax, LL, par) {
  
  ix1 <- 1:(Lmax * par$L)
  
  tmp      <- rep(0, par$L * par$ls$ncop)
  #tmp[ix1] <- (PRIOR TERM) * (LIKELIHOOD TERM)
  #gamma determines how much the agent weights previous beliefs (lsb) vs how 
  #much they uniformly distribute probability without consideration for lsb 
  #of the previous trial
  tmp[ix1] <- ((1 - par$ls$gamma) * lsb[ix1] +
                 par$ls$gamma / (Lmax * par$L)) *
    exp(LL[ix1] - max(LL[ix1]))  # exponentiate to change from log likelihood to regular likelihood, subtraction for numerical stability
  
  lsb <- tmp / sum(tmp)
  
  return(lsb)
}


# =========================================================================
# NewState: Check if a new latent state should be added (Page's algorithm)
# =========================================================================
NewState <- function(V, lsb, Lmax, c_vec, win, LL, q, par) {
  
  if (Lmax < par$ls$ncop) {
    
    ix1 <- 1:(par$L * Lmax)
    ix2 <- par$L * Lmax + (1:par$L)
    mx <- max(exp(LL[c(ix1, ix2)]))
    
    tmp <- rep(0, par$L * par$ls$ncop)
    tmp[ix1] <- ((1 - par$ls$gamma) * lsb[ix1] +
                   par$ls$gamma / (Lmax * par$L)) * exp(LL[ix1] - mx)
    l0 <- sum(tmp)
    
    l1 <- sum((1 / par$L) * exp(LL[ix2] - mx))
    
    q <- max(q + log(l1) - log(l0) - par$ls$delta, 0)
    
    if (q >= par$ls$eta) {
      
      Lmax <- Lmax + 1
      q <- 0
      
      tmp[c(ix1, ix2)] <- ((1 - par$ls$gamma) * lsb[c(ix1, ix2)] +
                             par$ls$gamma / (Lmax * par$L)) *
        exp(LL[c(ix1, ix2)])
      
      # FIXED: Force D x C matrix before rowSums
      c_l <- c_vec[, , Lmax]
      if (length(dim(c_l)) < 2) {
        dim(c_l) <- c(par$D, par$C)
      }
      
      denom <- rowSums(c_l^2)
      Vnew <- sweep(c_l, 1, win / denom, "*")
      V[, , Lmax] <- Vnew
      
      lsb <- tmp / sum(tmp)
    }
  }
  
  return(list(V = V, lsb = lsb, Lmax = Lmax, q = q))
}

#########################UPDATE

Update <- function(V, B, sigmasq, lsb, Lmax, c_vec, win, par) {
  
  L <- par$L * par$ls$ncop
  
  if (par$t == par$ntrials || is.null(par$time)) {
    nIter <- 1
  } else {
    nIter <- min(par$ls$chi, round(par$time[par$t + 1] - par$time[par$t]))
  }
  
  sigmasq0 <- sigmasq
  mu <- apply(V * c_vec, c(1, 3), sum)
  
  # GUARD: restore D x L dims if apply() dropped a dimension (e.g. when par$D == 1)
  if (is.null(dim(mu)) || length(dim(mu)) < 2) {
    dim(mu) <- c(par$D, L)
  }
  
  for (d in 1:par$D) {
    for (l1 in 1:par$L) { 
      for (l2 in 1:Lmax) {
        
        l <- (l2 - 1) * par$L + l1  # FLATTENED INDEX
        
        # Extract cue vectors for this state (must be C x 1)
        cue_row <- c_vec[d, , l]  # 1 x C → C-length vector
        if (is.null(dim(cue_row)) || length(dim(cue_row)) != 2) {
          cue_row <- matrix(cue_row, nrow = par$C, ncol = 1)
        }
        
        # Update effort matrix (C x C)
        outer_prod <- cue_row %*% t(cue_row)  # C x 1 * 1 x C → C x C
        B[, , d, l] <- B[, , d, l] + par$ls$alpha2 * (lsb[l] * outer_prod - B[, , d, l])
        
        # Associability (C x 1)
        tmp <- par$ls$alpha0 * lsb[l] * solve(B[, , d, l], cue_row)
        
        # TD error
        TD <- win[d] - mu[d, l]
        
        # Update V
        V[d, , l] <- V[d, , l] + tmp * TD
        
        sigmasq <- sigmasq + par$ls$alpha1 * lsb[l] * TD^2
        
        if (nIter > 1) {
          for (n in seq_len(nIter - 1)) {
            mu <- apply(V * c_vec, c(1, 3), sum)
            
            # GUARD: restore D x L dims after rumination recompute of mu
            if (is.null(dim(mu)) || length(dim(mu)) < 2) {
              dim(mu) <- c(par$D, L)
            }
            
            TD <- win[d] - mu[d, l]
            V[d, , l] <- V[d, , l] + tmp * TD
          }
        }
      }
    }
  }
  
  sigmasq <- sigmasq - par$ls$alpha1 * sigmasq0
  return(list(V = V, B = B, sigmasq = sigmasq))
}


# =========================================================================
# SoftMaxChoice: Select arm via softmax decision rule

SoftMaxChoice <- function(mu0, par) {
  vals  <- colSums(mu0) #sums V over all reward dimensions D
  scaled <- par$tau * vals #tau is temperature parameter, controls how deterministic softmax decision making is
  p     <- exp(scaled - max(scaled))
  p     <- p / sum(p)
  arm   <- which(runif(1) < cumsum(p))[1]
  return(arm)
}
