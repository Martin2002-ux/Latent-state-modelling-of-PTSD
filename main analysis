library(latentState)
library(ez)
library(dplyr)
library(emmeans)
library(afex)
library(tidyr)
library(truncnorm)

noise=0.2
#noise=0.05
#noise=0.1
#noise=0.3

#these are the baseline hyperparameters with variation, which were used for the
#simulation
alpha0="rbeta(1,13,50)"
alpha1="rbeta(1,20,362)"
alpha2="rbeta(1,10,182)"
gamma  = "0.01" #no variation
sigma0="rlnorm(1,meanlog = log(0.5), sdlog = 0.03)"
delta  = "0.6" #no variation
lambda="rlnorm(1,meanlog = log(1), sdlog = 0.18)"
eta="rtruncnorm(1, a = 0.7, b = Inf, mean = 1.2, sd = 1)"
eta_ptsd="rtruncnorm(1, a = 0.7, b = Inf, mean = 5, sd = 1)"

#when Mauchly's test is significant, we use the GG sphericity correction values
#this function also checks if the significance (yes or no) matches what is expected.
#Errors automatically become treated as if it was non-significant output
check_p_local <- function(ez, effect, sig) {
  if (is.null(ez)) {
    return(FALSE == sig)
  }

  mauchly <- ez$`Mauchly's Test for Sphericity`

  p <- if (!is.null(mauchly) && effect %in% mauchly$Effect &&
           mauchly %>% filter(Effect == effect) %>% pull(`p<.05`) == "*") {
    ez$`Sphericity Corrections` %>% filter(Effect == effect) %>% pull(`p[GG]`)
  } else {
    ez$ANOVA %>% filter(Effect == effect) %>% pull(p)
  }

  if (is.na(p)) {
    return(FALSE == sig)
  }

  return((p < 0.05) == sig)
}

#for non-ANOVA, this function will check if the significance of the test matches
#the intended significance (true or false)
check_p_other_local <- function(p, sig) {
  is_sig <- if (is.na(p)) FALSE else (p < 0.05)
  return(is_sig == sig)
}

#function that runs the power analysis, taking simulation details as input
run_power_analysis <- function(study_fun,
                               n_iter = 1000,
                               seed_start = 1000,
                               n_ctrl,
                               n_ptsd,
                               noise,
                               eta,
                               eta_ptsd,
                               lambda,
                               alpha0,
                               alpha1,
                               alpha2,
                               gamma,
                               delta,
                               sigma0,
                               exclude_from_joint = character(0)) {

  results_list <- lapply(seq_len(n_iter), function(i) {
    cat("Iteration", i, "/", n_iter, "\n")

    study_fun(seed = seed_start + i,
              n_ctrl = n_ctrl, n_ptsd = n_ptsd, noise = noise,
              eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
              alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
              gamma = gamma, delta = delta, sigma0 = sigma0)
  })

  results_mat <- do.call(rbind, results_list)

  marginal <- colMeans(results_mat)

  # joint: computed only over checks NOT in exclude_from_joint
  joint_cols  <- setdiff(colnames(results_mat), exclude_from_joint)
  joint_mat   <- results_mat[, joint_cols, drop = FALSE]
  joint <- mean(apply(joint_mat, 1, all))

  n_passed_per_run <- rowSums(results_mat)
  partial_dist <- table(n_passed_per_run)

  # Specific combinations of checks passed per iteration
  pass_patterns <- apply(
    results_mat,
    1,
    function(x) paste(names(x)[x], collapse = ";")
  )

  pattern_dist <- sort(table(pass_patterns), decreasing = TRUE)

  list(
    marginal = marginal,                  # per-check replication rate
    joint = joint,                        # all-checks-simultaneously replication rate
    n_passed_per_run = partial_dist,      # distribution of number of checks passed
    pass_patterns = pattern_dist,         # distribution of exact pass combinations
    raw = results_mat                     # raw TRUE/FALSE results
  )
}

#################################Physiological#################################

##########################
###Blechert et al. 2007###
##########################

study_blechert2007 <- function(seed, n_ctrl, n_ptsd, noise,
                               eta, eta_ptsd, lambda,
                               alpha0, alpha1, alpha2, gamma, delta, sigma0) {

  phase_def <- list(
    phase1 = list(
      trials  = list(A = 6, B = 6),
      rewards = list(c(A = 1))
    ),
    phase2 = list(
      trials = list(A = 6, B = 6)
    )
  )

  ctrl <- run(phase_def = phase_def, seed = seed, V_noise_sd = noise,
              n = n_ctrl, block = 6, eta = eta, lambda = lambda,
              alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
              gamma = gamma, delta = delta, sigma0 = sigma0)

  ptsd <- run(phase_def = phase_def, seed = seed + 1e6, V_noise_sd = noise,
              n = n_ptsd, block = 6, eta = eta_ptsd, lambda = lambda,
              alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
              gamma = gamma, delta = delta, sigma0 = sigma0,
              prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  df_avg <- df_long %>%
    group_by(subject, group, cue_type, block) %>%
    summarise(V = mean(V), .groups = "drop")

  ext     <- ezANOVA(data = df_avg, dv = V, wid = subject,
                     within = c(cue_type, block), between = group)
  CSplus  <- ezANOVA(data = filter(df_avg, cue_type == "A"), dv = V,
                     wid = subject, within = block, between = group)
  CSmin   <- ezANOVA(data = filter(df_avg, cue_type == "B"), dv = V,
                     wid = subject, within = block, between = group)

  c(
    extinction_ANOVA_group_main_effect   = check_p_local(ext, "group", TRUE),
    extinction_ANOVA_cue_type_main_effect       = check_p_local(ext, "cue_type", TRUE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext, "group:cue_type", TRUE),
    CS_plus_group_difference  = check_p_local(CSplus, "group", TRUE),
    CS_minus_group_difference   = check_p_local(CSmin, "group", FALSE)
  )
}

blechert2007 <- run_power_analysis(
  study_fun  = study_blechert2007,
  n_ctrl = 34, n_ptsd = 36,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "blechert2007",
  n_ctrl = 34,
  n_ptsd = 36,
  Marginal = names(blechert2007$marginal),
  Value = unname(blechert2007$marginal),
  Joint = blechert2007$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("blechert2007_prob_", noise, ".csv"))

max_passes <- length(blechert2007$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- blechert2007$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "blechert2007",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("blechert2007_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "blechert2007",
  Pass_pattern = names(blechert2007$pass_patterns),
  n_passes = sapply(
    strsplit(names(blechert2007$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(blechert2007$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("blechert2007_pattern_", noise, ".csv"))

############################
###Felmingham et al. 2018###
############################

study_felmingham2018 <- function(seed, n_ctrl, n_ptsd, noise,
                               eta, eta_ptsd, lambda,
                               alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 5,
        B = 5
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        B = 10
      )
    )
  )


  ctrl=run(seed=seed, phase_def, V_noise_sd=noise, n=n_ctrl, pseudo=2, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def, V_noise_sd=noise, n=n_ptsd, pseudo=2, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, group="ptsd", prev=n_ctrl)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) <= 5)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) <= 5))
  ))

  early = ezANOVA(data = df_long, dv = V, wid = subject,
                  within = c(cue_type, presentation), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) > 5)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) > 5))
  ))

  late = ezANOVA(data = df_long, dv = V, wid = subject,
                 within = c(cue_type, presentation), between = group)

  c(
    early_extinction_ANOVA_cue_type_main_effect = check_p_local(early, "cue_type", TRUE),
    early_extinction_ANOVA_trial_number_main_effect = check_p_local(early, "presentation", TRUE),
    early_extinction_ANOVA_group_x_trial_interaction = check_p_local(early, "group:presentation", TRUE),
    late_extinction_ANOVA_trial_number_main_effect = check_p_local(late, "presentation", TRUE),
    late_extinction_ANOVA_group_x_trial_interaction = check_p_local(late, "group:presentation", FALSE),
    late_extinction_ANOVA_stimulu_main_effects = check_p_local(late, "cue_type", FALSE)
  )
}

felmingham2018 <- run_power_analysis(
  study_fun  = study_felmingham2018,
  n_ctrl = 84, n_ptsd = 22,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "felmingham2018",
  n_ctrl = 84,
  n_ptsd = 22,
  Marginal = names(felmingham2018$marginal),
  Value = unname(felmingham2018$marginal),
  Joint = felmingham2018$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("felmingham2018_prob_", noise, ".csv"))

max_passes <- length(felmingham2018$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- felmingham2018$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "felmingham2018",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("felmingham2018_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "felmingham2018",
  Pass_pattern = names(felmingham2018$pass_patterns),
  n_passes = sapply(
    strsplit(names(felmingham2018$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(felmingham2018$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("felmingham2018_pattern_", noise, ".csv"))

###########################
###Pohlchen et al., 2020###
###########################

study_pohlechen2020 <- function(seed, n_ctrl, n_ptsd, noise,
                                 eta, eta_ptsd, lambda,
                                 alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12,
        C = 12
      ),
      rewards = list(
        c(A = 0.75, B = 0.75)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        C = 10
      )
    ),

    day2 = list(
      trials = list(
        A = 8,
        B = 8,
        C = 8
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=c(12,10,12), ezITI=list(c(56,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=c(12,10,12), ezITI=list(c(56,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]])),
    lapply(ptsd, function(r) subset(r$long[[1]]))
  ))

  #the second CS+ was not included in analyses
  df_long = df_long %>%
    filter(cue_type != "B")

  acq = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, block), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  ext = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, block), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[3]])),
    lapply(ptsd, function(r) subset(r$long[[3]]))
  ))

  df_long = df_long %>%
    filter(cue_type != "B")

  rec = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, block), between = group)

  c(
    acquisition_ANOVA_group_main_effect = check_p_local(acq, "group", FALSE),
    acquisition_ANOVA_group_x_cue_type_interaction = check_p_local(acq, "group:cue_type", FALSE),
    extinction_ANOVA_group_main_effect = check_p_local(ext, "group", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext, "group:cue_type", FALSE),

    recall_ANOVA_group_main_effect = check_p_local(rec, "group", FALSE),
    recall_ANOVA_group_x_cue_type_interaction = check_p_local(rec, "group:cue_type", FALSE)
  )
}

pohlechen2020 <- run_power_analysis(
  study_fun  = study_pohlechen2020,
  n_ctrl = 35, n_ptsd = 21,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = c("recall_ANOVA_group_main_effect", "recall_ANOVA_group_x_cue_type_interaction")
)

df_marginal <- data.frame(
  Name = "pohlechen2020",
  n_ctrl = 35,
  n_ptsd = 21,
  Marginal = names(pohlechen2020$marginal),
  Value = unname(pohlechen2020$marginal),
  Joint = pohlechen2020$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("pohlechen2020_prob_", noise, ".csv"))

max_passes <- length(pohlechen2020$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- pohlechen2020$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "pohlechen2020",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("pohlechen2020_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "pohlechen2020",
  Pass_pattern = names(pohlechen2020$pass_patterns),
  n_passes = sapply(
    strsplit(names(pohlechen2020$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(pohlechen2020$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("pohlechen2020_pattern_", noise, ".csv"))

##########################
###Norrholm et al. 2011###
##########################

study_norrholm2011 <- function(seed, n_ctrl, n_ptsd, noise,
                                eta, eta_ptsd, lambda,
                                alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 24,
        B = 24
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=8, eta=eta_ptsd, alpha1=alpha1, lambda=lambda, alpha0=alpha0, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]], as.integer(presentation) > 8)),
    lapply(ptsd, function(r) subset(r$long[[1]], as.integer(presentation) > 8))
  ))

  df_acq_avg <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V, na.rm = TRUE), .groups = "drop")

  acq = ezANOVA(data = df_acq_avg, dv = V, wid = subject,
                within = c(cue_type), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  #there are 6 blocks which are combined to early, middle and late
  df_avg <- df_long %>%
    mutate(superblock = factor(ceiling(as.integer(block) / 2))) %>%
    group_by(subject, group, cue_type, superblock) %>%
    summarise(V = mean(V), .groups = "drop")

  ext <- ezANOVA(data = df_avg, dv = V, wid = subject,
                 within = c(cue_type, superblock), between = group)

  df_long_CSminus <- df_avg %>% filter(cue_type == "B")

  ext_CSminus = ezANOVA(data = df_long_CSminus, dv = V, wid = subject,
                        within = c(superblock), between = group)

  #t tests finding out which CS+ superblock had a significant difference
  ext_sb1 <- t.test(V ~ group, data = df_avg %>% filter(superblock == 1, cue_type == "A"))
  ext_sb2 <- t.test(V ~ group, data = df_avg %>% filter(superblock == 2, cue_type == "A"))
  ext_sb3 <- t.test(V ~ group, data = df_avg %>% filter(superblock == 3, cue_type == "A"))

  c(
    extinction_ANOVA_group_main_effect = check_p_local(ext, "group", TRUE),
    extinction_block_1_t_test = check_p_other_local(ext_sb1$p.value, TRUE),
    extinction_block_2_t_test = check_p_other_local(ext_sb2$p.value, TRUE),
    acquisition_ANOVA_group_main_effect = check_p_local(acq, "group", FALSE),
    CS_minus_group_difference = check_p_local(ext_CSminus, "group", FALSE),
    extinction_block_3_t_test = check_p_other_local(ext_sb3$p.value, FALSE)
  )
}

norrholm2011 <- run_power_analysis(
  study_fun  = study_norrholm2011,
  n_ctrl = 78, n_ptsd = 49,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "norrholm2011",
  n_ctrl = 78,
  n_ptsd = 49,
  Marginal = names(norrholm2011$marginal),
  Value = unname(norrholm2011$marginal),
  Joint = norrholm2011$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("norrholm2011_prob_", noise, ".csv"))

max_passes <- length(norrholm2011$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- norrholm2011$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "norrholm2011",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("norrholm2011_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "norrholm2011",
  Pass_pattern = names(norrholm2011$pass_patterns),
  n_passes = sapply(
    strsplit(names(norrholm2011$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(norrholm2011$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("norrholm2011_pattern_", noise, ".csv"))

#######################
###Shvil et al. 2014###
#######################

study_shvil2014 <- function(seed, n_ctrl, n_ptsd, noise,
                               eta, eta_ptsd, lambda,
                               alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  #X represents context A, Y represents context B
  phase_def <- list(

    phase1_block1 = list(
      trials = list(
        AX = 8,
        CX = 8
      ),
      rewards = list(
        c(AX = 0.6)
      )
    ),

    phase1_block2 = list(
      trials = list(
        BX = 8,
        CX = 8
      ),
      rewards = list(
        c(BX = 0.6)
      )
    ),

    phase2 = list(
      trials = list(
        AY = 16,
        CY = 16
      )
    ),

    day2_block1 = list(
      trials = list(
        AY = 16,
        CY = 16
      )
    ),

    day2_block2 = list(
      trials = list(
        BY = 16,
        CY = 16
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n_ctrl, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n_ptsd, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[3]], as.integer(presentation) > 4)),
    lapply(ptsd, function(r) subset(r$long[[3]], as.integer(presentation) > 4))
  ))

  df_avg <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop")

  ext = ezANOVA(data = df_avg, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #recall analysis, comparing the two conditioned stimuli in day 2
  df_recall <- do.call(rbind, c(
    lapply(ctrl, function(r) {
      ay <- subset(r$long[[4]], cue_type == "AY" & as.integer(presentation) <= 4)
      by <- subset(r$long[[5]], cue_type == "BY" & as.integer(presentation) <= 4)
      rbind(ay, by)
    }),
    lapply(ptsd, function(r) {
      ay <- subset(r$long[[4]], cue_type == "AY" & as.integer(presentation) <= 4)
      by <- subset(r$long[[5]], cue_type == "BY" & as.integer(presentation) <= 4)
      rbind(ay, by)
    })
  ))

  df_recall_avg <- df_recall %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop")

  rec = ezANOVA(data = df_recall_avg, dv = V, wid = subject,
                within = cue_type, between = group)

  c(
    extinction_ANOVA_group_main_effect = check_p_local(ext, "group", FALSE),
    extinction_ANOVA_cue_type_main_effect = check_p_local(ext, "cue_type", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext, "group:cue_type", FALSE),
    recall_ANOVA_group_main_effect = check_p_local(rec, "group", FALSE),
    recall_ANOVA_cue_type_main_effect = check_p_local(rec, "cue_type", FALSE),
    recall_ANOVA_group_x_cue_type_interaction = check_p_local(rec, "group:cue_type", FALSE)
  )
}

shvil2014 <- run_power_analysis(
  study_fun  = study_shvil2014,
  n_ctrl = 25, n_ptsd = 31,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = c("recall_ANOVA_group_main_effect", "recall_ANOVA_cue_type_main_effect", "recall_ANOVA_group_x_cue_type_interaction")
)

df_marginal <- data.frame(
  Name = "shvil2014",
  n_ctrl = 25,
  n_ptsd = 31,
  Marginal = names(shvil2014$marginal),
  Value = unname(shvil2014$marginal),
  Joint = shvil2014$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("shvil2014_prob_", noise, ".csv"))

max_passes <- length(shvil2014$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- shvil2014$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "shvil2014",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("shvil2014_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "shvil2014",
  Pass_pattern = names(shvil2014$pass_patterns),
  n_passes = sapply(
    strsplit(names(shvil2014$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(shvil2014$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("shvil2014_pattern_", noise, ".csv"))

##########################
###Acheson et al. 2015####
##########################

study_acheson2015 <- function(seed, n_ctrl, n_ptsd, noise,
                            eta, eta_ptsd, lambda,
                            alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 8,
        B = 8
      ),
      rewards = list(
        c(A = 0.75)
      )
    ),

    phase2 = list(
      trials = list(
        A = 16,
        B = 16
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=8, eta=eta, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, lambda=lambda)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=8, eta=eta_ptsd, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, lambda=lambda, prev=n_ctrl, group="ptsd")


  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]], as.integer(presentation) > 6)),
    lapply(ptsd, function(r) subset(r$long[[1]], as.integer(presentation) > 6))
  ))

  df_avg_last_two_acq <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop")

  acq = ezANOVA(data = df_avg_last_two_acq, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #follow up t-test comparing A vs B in ctrl
  V_A <- df_avg_last_two_acq %>% filter(group == "ctrl", cue_type == "A") %>% arrange(subject) %>% pull(V)
  V_B <- df_avg_last_two_acq %>% filter(group == "ctrl", cue_type == "B") %>% arrange(subject) %>% pull(V)
  ctrl_ttest = t.test(V_A, V_B, paired = TRUE)

  #same t-test for ptsd
  V_A <- df_avg_last_two_acq %>% filter(group == "ptsd", cue_type == "A") %>% arrange(subject) %>% pull(V)
  V_B <- df_avg_last_two_acq %>% filter(group == "ptsd", cue_type == "B") %>% arrange(subject) %>% pull(V)
  ptsd_ttest = t.test(V_A, V_B, paired = TRUE)

  #the extinction test involves comparing cue A to the highest V in acq to get
  #a % conditioned fear. The following data wrangling is based on this goal
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]])),
    lapply(ptsd, function(r) subset(r$long[[1]]))
  ))

  #V can go as low as -0.5. To make % fear a sensible conversion, we need all values
  #to be positive, so we add 1 so that V = 0 becomes the lowest possible score
  df_long <- df_long %>%
    mutate(V = V + 0.5)

  #we want the V of the highest block average in acq
  df_avg_acq_block <- df_long %>%
    group_by(subject, group, cue_type, block) %>%
    summarise(V = mean(V), .groups = "drop")

  peak_acq <- df_avg_acq_block %>%
    group_by(subject) %>%
    summarise(peak_V = max(V), .groups = "drop")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  df_long <- df_long %>%
    mutate(V = V + 0.5)

  #we want block average for extinction, only looking at cue A
  df_avg_ext_block <- df_long %>%
    group_by(subject, group, cue_type, block) %>%
    summarise(V = mean(V), .groups = "drop") %>%
    filter(cue_type == "A")

  prcnt_fear_ext_block <- df_avg_ext_block %>%
    left_join(peak_acq, by = "subject") %>%
    mutate(V_pct = V / peak_V * 100) %>%
    dplyr::select(-peak_V)

  ext = ezANOVA(data = prcnt_fear_ext_block, dv = V_pct, wid = subject,
                within = block, between = group)

  #follow-up ttests were done on each block
  ext_block_1 <- t.test(V_pct ~ group, data = prcnt_fear_ext_block %>% filter(block == 1))
  ext_block_2 <- t.test(V_pct ~ group, data = prcnt_fear_ext_block %>% filter(block == 2))
  ext_block_3 <- t.test(V_pct ~ group, data = prcnt_fear_ext_block %>% filter(block == 3))
  ext_block_4 <- t.test(V_pct ~ group, data = prcnt_fear_ext_block %>% filter(block == 4))

  c(
    acquisition_ANOVA_group_x_cue_type_interaction = check_p_local(acq, "group:cue_type", TRUE),
    extinction_ANOVA_group_main_effect = check_p_local(ext, "group", TRUE),
    control_group_t_test_between_cue_type_types_during_acquisition = check_p_other_local(ctrl_ttest$p.value, TRUE),
    extinction_block_3_t_test = check_p_other_local(ext_block_3$p.value, TRUE),
    extinction_block_4_t_test = check_p_other_local(ext_block_4$p.value, TRUE),
    ptsd_group_t_test_between_cue_type_types_during_acquisition = check_p_other_local(ptsd_ttest$p.value, FALSE), #ptsd only had p value of 0.09, driven by higher CS- response in acq. There was no group differences in expectancy ratings
    extinction_block_1_t_test = check_p_other_local(ext_block_1$p.value, FALSE),
    extinction_block_2_t_test = check_p_other_local(ext_block_2$p.value, FALSE)
  )
}

acheson2015 <- run_power_analysis(
  study_fun  = study_acheson2015,
  n_ctrl = 923, n_ptsd = 42,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "acheson2015",
  n_ctrl = 923,
  n_ptsd = 42,
  Marginal = names(acheson2015$marginal),
  Value = unname(acheson2015$marginal),
  Joint = acheson2015$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("acheson2015_prob_", noise, ".csv"))

max_passes <- length(acheson2015$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- acheson2015$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "acheson2015",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("acheson2015_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "acheson2015",
  Pass_pattern = names(acheson2015$pass_patterns),
  n_passes = sapply(
    strsplit(names(acheson2015$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(acheson2015$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("acheson2015_pattern_", noise, ".csv"))

#########################
###Pineles et al. 2016###
#########################

study_pineles2016 <- function(seed, n_ctrl, n_ptsd, noise,
                              eta, eta_ptsd, lambda,
                              alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 5,
        B = 5
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        B = 10
      )
    ),

    day2 = list(
      trials = list(
        A = 5,
        B = 5
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, pseudo=3, ezITI=list(c(30,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, pseudo=3, ezITI=list(c(30,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #mean difference between last two trials of acquisition for cue A vs cue B for each participant
  acq_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[4:5], na.rm = TRUE) -
                       mean(x$V_cue$B[4:5], na.rm = TRUE))

  acq_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[4:5], na.rm = TRUE) -
                       mean(x$V_cue$B[4:5], na.rm = TRUE))

  #mean difference on trials 2 to 5 of extinction
  early_ext_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[7:10], na.rm = TRUE) -
                             mean(x$V_cue$B[7:10], na.rm = TRUE))

  early_ext_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[7:10], na.rm = TRUE) -
                             mean(x$V_cue$B[7:10], na.rm = TRUE))

  #mean difference on trials 6 to 10 of extinction
  late_ext_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[11:15], na.rm = TRUE) -
                            mean(x$V_cue$B[11:15], na.rm = TRUE))

  late_ext_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[11:15], na.rm = TRUE) -
                            mean(x$V_cue$B[11:15], na.rm = TRUE))

  #mean difference in extinction recall
  rec_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[16:20], na.rm = TRUE) -
                       mean(x$V_cue$B[16:20], na.rm = TRUE))

  rec_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[16:20], na.rm = TRUE) -
                       mean(x$V_cue$B[16:20], na.rm = TRUE))

  #calculating differences between phases, which are used as the outcome variables
  initial_change_ctrl <- early_ext_ctrl - acq_ctrl
  initial_change_ptsd <- early_ext_ptsd - acq_ptsd

  late_change_ctrl <- late_ext_ctrl - early_ext_ctrl
  late_change_ptsd <- late_ext_ptsd - early_ext_ptsd

  ext_retention_ctrl <- rec_ctrl - early_ext_ctrl
  ext_retention_ptsd <- rec_ptsd - early_ext_ptsd

  df <- data.frame(
    initial_change = c(initial_change_ctrl, initial_change_ptsd),
    late_change = c(late_change_ctrl, late_change_ptsd),
    ext_retention = c(ext_retention_ctrl, ext_retention_ptsd),
    group = factor(c(rep("ctrl", n_ctrl),
                     rep("ptsd", n_ptsd)))
  )

  fit_initial <- lm(initial_change ~ group, data = df)
  fit_late <- lm(late_change ~ group, data = df)
  fit_ext_retention <- lm(ext_retention ~ group, data = df)

  initial_p = summary(fit_initial)$coefficients["groupptsd", "Pr(>|t|)"]
  late_p = summary(fit_late)$coefficients["groupptsd", "Pr(>|t|)"]
  retention_p = summary(fit_ext_retention)$coefficients["groupptsd", "Pr(>|t|)"]
  retention_est <- summary(fit_ext_retention)$coefficients["groupptsd", "Estimate"]

  c(
    group_effect_on_regression_for_initial_fear_change = check_p_other_local(initial_p, FALSE),
    group_effect_on_regression_for_late_fear_change = check_p_other_local(late_p, FALSE),
    group_effect_on_regression_for_fear_recall = (retention_p < 0.05) && (retention_est > 0)
  )
}

pineles2016 <- run_power_analysis(
  study_fun  = study_pineles2016,
  n_ctrl = 16, n_ptsd = 16,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = "group_effect_on_regression_for_fear_recall"
)

df_marginal <- data.frame(
  Name = "pineles2016",
  n_ctrl = 16,
  n_ptsd = 16,
  Marginal = names(pineles2016$marginal),
  Value = unname(pineles2016$marginal),
  Joint = pineles2016$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("pineles2016_prob_", noise, ".csv"))

max_passes <- length(pineles2016$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- pineles2016$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "pineles2016",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("pineles2016_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "pineles2016",
  Pass_pattern = names(pineles2016$pass_patterns),
  n_passes = sapply(
    strsplit(names(pineles2016$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(pineles2016$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("pineles2016_pattern_", noise, ".csv"))

##########################
###Norrholm et al. 2015###
##########################

study_norrholm2015 <- function(seed, n_ctrl, n_ptsd, noise,
                              eta, eta_ptsd, lambda,
                              alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 24,
        B = 24
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block = 8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block = 8, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #mean V of cue A during early extinction
  early_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[13:20], na.rm = TRUE))
  early_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[13:20], na.rm = TRUE))

  #mean V of cue A during mid extinction
  mid_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[21:28], na.rm = TRUE))
  mid_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[21:28], na.rm = TRUE))

  #mean V of cue A during late extinction
  late_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[29:36], na.rm = TRUE))
  late_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[29:36], na.rm = TRUE))

  df <- data.frame(
    early = c(early_ctrl, early_ptsd),
    mid = c(mid_ctrl, mid_ptsd),
    late = c(late_ctrl, late_ptsd),
    group = factor(c(rep("ctrl", n_ctrl),
                     rep("ptsd", n_ptsd)))
  )

  early_p = t.test(early_ctrl, early_ptsd)
  mid_p = t.test(mid_ctrl, mid_ptsd)
  late_p = t.test(late_ctrl, late_ptsd)

  c(
    early_extinction_t_test_between_groups = check_p_other_local(early_p$p.value, TRUE),
    middle_extinction_t_test_between_groups = check_p_other_local(mid_p$p.value, TRUE),
    late_extinction_t_test_between_groups = check_p_other_local(late_p$p.value, FALSE)
  )
}

#This study looked at symptom severity rather than binary PTSD classification.
#As such, it doesn't have an explicit PTSD sample size. However, it does give
#the mean and sd of PSS symptom score. The PSS cutoff for PTSD is 29. Therefore,
#we can use this to calculate how many of the 269 participants we expect to be
#in the PTSD group

#Probability of X >= 29 for a normal distribution truncated at 0
prob <- 1 - ptruncnorm(
  29,
  a = 0,
  b = Inf,
  mean = 14.78,
  sd = 12.5
)

prob * 269

norrholm2015 <- run_power_analysis(
  study_fun  = study_norrholm2015,
  n_ctrl = 230, n_ptsd = 39,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "norrholm2015",
  n_ctrl = 230,
  n_ptsd = 39,
  Marginal = names(norrholm2015$marginal),
  Value = unname(norrholm2015$marginal),
  Joint = norrholm2015$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("norrholm2015_prob_", noise, ".csv"))

max_passes <- length(norrholm2015$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- norrholm2015$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "norrholm2015",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("norrholm2015_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "norrholm2015",
  Pass_pattern = names(norrholm2015$pass_patterns),
  n_passes = sapply(
    strsplit(names(norrholm2015$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(norrholm2015$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("norrholm2015_pattern_", noise, ".csv"))

#######################
###Handy et al. 2018###
#######################

study_handy2018 <- function(seed, n_ctrl, n_ptsd, noise,
                               eta, eta_ptsd, lambda,
                               alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 60
      ),
      rewards = list(
        c(A = 0.5)
      )
    ),

    phase2 = list(
      trials = list(
        A = 20
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, pseudo=3, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, pseudo=3, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #mean of last 20 acquisition trials
  acq_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[51:60], na.rm = TRUE))
  acq_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[51:60], na.rm = TRUE))

  #mean of extinction
  ext_ctrl <- sapply(ctrl, function(x) mean(x$V_cue$A[61:80], na.rm = TRUE))
  ext_ptsd <- sapply(ptsd, function(x) mean(x$V_cue$A[61:80], na.rm = TRUE))

  # Combine into one data frame with group labels
  df <- data.frame(
    group = factor(
      c(rep("ctrl", length(acq_ctrl)),
        rep("ptsd", length(acq_ptsd)))
    ),
    acq = c(acq_ctrl, acq_ptsd),
    ext = c(ext_ctrl, ext_ptsd)
  )

  #ANCOVA
  ext_test <- aov(ext ~ group * acq, data = df)
  ext_p=summary(ext_test)[[1]][[1, "Pr(>F)"]]

  c(
    extinction_group_effect_on_regression = check_p_other_local(ext_p, TRUE)
  )
}

handy2018 <- run_power_analysis(
  study_fun  = study_handy2018,
  n_ctrl = 51, n_ptsd = 15,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "handy2018",
  n_ctrl = 51,
  n_ptsd = 15,
  Marginal = names(handy2018$marginal),
  Value = unname(handy2018$marginal),
  Joint = handy2018$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("handy2018_prob_", noise, ".csv"))

max_passes <- length(handy2018$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- handy2018$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "handy2018",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("handy2018_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "handy2018",
  Pass_pattern = names(handy2018$pass_patterns),
  n_passes = sapply(
    strsplit(names(handy2018$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(handy2018$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("handy2018_pattern_", noise, ".csv"))

#########################
###Burriss et al. 2006###
#########################

study_burriss2006 <- function(seed, n_ctrl, n_ptsd, noise,
                            eta, eta_ptsd, lambda,
                            alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 60,
        B = 60
      ),
      rewards = list(
        c(A = 1) #THE LAST CUE A OF EACH BLOCK SHOULD BE UNREWARDED
      )
    ),

    phase2 = list(
      trials = list(
        A = 20,
        B = 20
      )
    )
  )

  #We cannot automatically create experiments where specifically the last cue A trial
  #of each block is unrewarded, so instead we manually create a random trial order
  #which matches the requirements separately for each participant

  build_phase1_world <- function(n_blocks = 12, block_size = 10, ncop = 10) {
    p1_trials  <- character(0)
    p1_rewards <- numeric(0)

    for (i in seq_len(n_blocks)) {
      trials  <- sample(c(rep("A", 5), rep("B", 5)))
      rewards <- ifelse(trials == "A", 1, 0)
      rewards[max(which(trials == "A"))] <- 0   #last A in block gets no reward
      p1_trials  <- c(p1_trials,  trials)
      p1_rewards <- c(p1_rewards, rewards)
    }

    p2_trials  <- sample(c(rep("A", 20), rep("B", 20)))
    p2_rewards <- rep(0, 40)

    all_trials  <- c(p1_trials, p2_trials)
    all_rewards <- c(p1_rewards, p2_rewards)

    w <- make_world(trial_patterns = all_trials, reward = all_rewards, L = ncop)
    w$type_seq <- list(p1_trials, p2_trials)
    w
  }

  #ctrl group
  ctrl <- vector("list", n_ctrl)
  for (i in seq_len(n_ctrl)) {
    set.seed(seed+i*100)
    my_world_i <- build_phase1_world()
    ctrl[[i]] <- run(
      seed=seed+i*100, phase_def = phase_def, V_noise_sd=noise, n = 1, my_world = my_world_i,
      eta = eta, lambda = lambda, alpha0 = alpha0, alpha1 = alpha1,
      alpha2 = alpha2, gamma = gamma, delta = delta, sigma0 = sigma0,
      prev = i - 1, group = "ctrl"
    )[[1]]  # run() returns a list of length n, extract the single result
  }

  #ptsd group
  ptsd <- vector("list", n_ptsd)
  for (i in seq_len(n_ptsd)) {
    set.seed(seed+i*100+1e6)
    my_world_i <- build_phase1_world()
    ptsd[[i]] <- run(
      seed=seed+i*100+1e6, phase_def = phase_def, V_noise_sd=noise, n = 1, my_world = my_world_i,
      eta = eta_ptsd, lambda = lambda, alpha0 = alpha0, alpha1 = alpha1,
      alpha2 = alpha2, gamma = gamma, delta = delta, sigma0 = sigma0,
      prev = n_ctrl + i - 1, group = "ptsd"
    )[[1]]
  }

  #add block labels, because our manual world creation doesn't include them
  ctrl <- lapply(ctrl, function(r) {
    r$long[[1]]$block <- factor(ceiling(as.integer(r$long[[1]]$presentation) / 5))
    r$long[[2]]$block <- factor(ceiling(as.integer(r$long[[2]]$presentation) / 5))
    r
  })

  ptsd <- lapply(ptsd, function(r) {
    r$long[[1]]$block <- factor(ceiling(as.integer(r$long[[1]]$presentation) / 5))
    r$long[[2]]$block <- factor(ceiling(as.integer(r$long[[2]]$presentation) / 5))
    r
  })

  #RM ANOVA acquisition
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]])),
    lapply(ptsd, function(r) subset(r$long[[1]]))
  ))

  acq = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, block), between = group)

  #RM ANOVA extinction
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  ext = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, block), between = group)

  c(
    extinction_ANOVA_block_main_effect = check_p_local(ext, "block", TRUE),
    acquistion_ANOVA_group_x_block_interaction = check_p_local(acq, "group:block", FALSE),
    extinction_ANOVA_group_x_block_interaction = check_p_local(ext, "group:block", FALSE),
    extinction_ANOVA_cue_type_x_block_interaction = check_p_local(ext, "cue_type:block", FALSE)
  )
}

burriss2006 <- run_power_analysis(
  study_fun  = study_burriss2006,
  n_ctrl = 51, n_ptsd = 29,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "burriss2006",
  n_ctrl = 51,
  n_ptsd = 29,
  Marginal = names(burriss2006$marginal),
  Value = unname(burriss2006$marginal),
  Joint = burriss2006$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("burriss2006_prob_", noise, ".csv"))

max_passes <- length(burriss2006$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- burriss2006$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "burriss2006",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("burriss2006_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "burriss2006",
  Pass_pattern = names(burriss2006$pass_patterns),
  n_passes = sapply(
    strsplit(names(burriss2006$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(burriss2006$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("burriss2006_pattern_", noise, ".csv"))

#############################
###Grillon and Morgan 1999###
#############################

study_grillon1999 <- function(seed, n_ctrl, n_ptsd, noise,
                              eta, eta_ptsd, lambda,
                              alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 10,
        B = 10
      ),
      rewards = list(
        c(A = 0.8)
      )
    ),

    phase2 = list(
      trials = list(
        A = 6,
        B = 6
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=c(10,4), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=c(10,4), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #RM ANOVA acquisition
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]])),
    lapply(ptsd, function(r) subset(r$long[[1]]))
  ))

  #first we do an omnibus ANOVA, with afex package so it's compatible with emmeans
  acq <- aov_ez(
    id      = "subject",
    dv      = "V",
    data    = df_long,
    within  = c("cue_type", "block"),
    between = "group"
  )

  #now we do a simple effects analysis using the omnibus anova, which allows us
  #to test cue_type within each group while maintaining power from the entire
  #sample. This is the approach used in the original study.
  emm_acq <- emmeans(acq, ~ cue_type | group)
  pairs(emm_acq)

  ctrl_acq = summary(pairs(emm_acq))$p.value[summary(pairs(emm_acq))$group == "ctrl"]
  ptsd_acq = summary(pairs(emm_acq))$p.value[summary(pairs(emm_acq))$group == "ptsd"]

  #RM ANOVA extinction
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  #first we do an omnibus ANOVA, with afex package so it's compatible with emmeans
  ext <- aov_ez(
    id      = "subject",
    dv      = "V",
    data    = df_long,
    within  = c("cue_type", "block"),
    between = "group"
  )

  #now we do a simple effects analysis using the omnibus anova, which allows us
  #to test cue_type within each group while maintaining power from the entire
  #sample. This is the approach used in the original study.
  emm_ext <- emmeans(ext, ~ cue_type | group)
  pairs(emm_ext)

  ctrl_ext = summary(pairs(emm_ext))$p.value[summary(pairs(emm_ext))$group == "ctrl"]
  ptsd_ext = summary(pairs(emm_ext))$p.value[summary(pairs(emm_ext))$group == "ptsd"]

  c(
    control_group_effect_of_cue_type_in_acquisition = check_p_other_local(ctrl_acq, TRUE),
    control_group_effect_of_cue_type_in_extinction = check_p_other_local(ctrl_ext, TRUE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_other_local(ext$anova_table["group:cue_type", "Pr(>F)"], TRUE),
    acquisition_ANOVA_group_x_cue_type_interaction = check_p_other_local(acq$anova_table["group:cue_type", "Pr(>F)"], FALSE),
    PTSD_group_effect_of_cue_type_in_acquisition = check_p_other_local(ptsd_acq, FALSE),
    PTSD_group_effect_of_cue_type_in_extinction = check_p_other_local(ptsd_ext, FALSE)
  )
}

grillon1999 <- run_power_analysis(
  study_fun  = study_grillon1999,
  n_ctrl = 12, n_ptsd = 12,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "grillon1999",
  n_ctrl = 12,
  n_ptsd = 12,
  Marginal = names(grillon1999$marginal),
  Value = unname(grillon1999$marginal),
  Joint = grillon1999$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("grillon1999_prob_", noise, ".csv"))

max_passes <- length(grillon1999$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- grillon1999$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "grillon1999",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("grillon1999_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "grillon1999",
  Pass_pattern = names(grillon1999$pass_patterns),
  n_passes = sapply(
    strsplit(names(grillon1999$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(grillon1999$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("grillon1999_pattern_", noise, ".csv"))

#######################
###Milad et al. 2009###
#######################

study_milad2009 <- function(seed, n_ctrl, n_ptsd, noise,
                              eta, eta_ptsd, lambda,
                              alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  #X represents context A, Y represents context B
  phase_def <- list(

    phase1_block1 = list(
      trials = list(
        AX = 8,
        CX = 8
      ),
      rewards = list(
        c(AX = 0.6)
      )
    ),

    phase1_block2 = list(
      trials = list(
        BX = 8,
        CX = 8
      ),
      rewards = list(
        c(BX = 0.6)
      )
    ),

    phase2 = list(
      trials = list(
        AY = 16,
        CY = 16
      )
    ),

    day2_block1 = list(
      trials = list(
        AY = 16,
        CY = 16
      )
    ),

    day2_block2 = list(
      trials = list(
        BY = 16,
        CY = 16
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n_ctrl, lambda=lambda, eta=eta, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3) #I guessed pseudo since they only said "pseudorandomization" without a specific number
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n_ptsd, lambda=lambda, eta=eta_ptsd, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[3]], as.integer(presentation) > 4)),
    lapply(ptsd, function(r) subset(r$long[[3]], as.integer(presentation) > 4))
  ))

  df_avg <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop")

  ext = ezANOVA(data = df_avg, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #recall analysis, comparing the two conditioned stimuli in day 2
  df_recall <- do.call(rbind, c(
    lapply(ctrl, function(r) {
      ay <- subset(r$long[[4]], cue_type == "AY" & as.integer(presentation) <= 4)
      by <- subset(r$long[[5]], cue_type == "BY" & as.integer(presentation) <= 4)
      rbind(ay, by)
    }),
    lapply(ptsd, function(r) {
      ay <- subset(r$long[[4]], cue_type == "AY" & as.integer(presentation) <= 4)
      by <- subset(r$long[[5]], cue_type == "BY" & as.integer(presentation) <= 4)
      rbind(ay, by)
    })
  ))

  df_recall_avg <- df_recall %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop")

  rec = ezANOVA(data = df_recall_avg, dv = V, wid = subject,
                within = cue_type, between = group)

  c(
    extinction_ANOVA_group_main_effect = check_p_local(ext, "group", FALSE),
    extinction_ANOVA_cue_type_main_effect = check_p_local(ext, "cue_type", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext, "group:cue_type", FALSE),
    recall_ANOVA_group_x_cue_type_interaction = check_p_local(rec, "group:cue_type", TRUE)
  )

}

milad2009 <- run_power_analysis(
  study_fun  = study_milad2009,
  n_ctrl = 19, n_ptsd = 20,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = "recall_ANOVA_group_x_cue_type_interaction"
)

df_marginal <- data.frame(
  Name = "milad2009",
  n_ctrl = 19,
  n_ptsd = 20,
  Marginal = names(milad2009$marginal),
  Value = unname(milad2009$marginal),
  Joint = milad2009$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("milad2009_prob_", noise, ".csv"))

max_passes <- length(milad2009$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- milad2009$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "milad2009",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("milad2009_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "milad2009",
  Pass_pattern = names(milad2009$pass_patterns),
  n_passes = sapply(
    strsplit(names(milad2009$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(milad2009$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("milad2009_pattern_", noise, ".csv"))

#####################
###Zuj et al. 2017###
#####################

study_zuj2017 <- function(seed, n_ctrl, n_ptsd, noise,
                            eta, eta_ptsd, lambda,
                            alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 5,
        B = 5
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        B = 10
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=10, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=10, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]])),
    lapply(ptsd, function(r) subset(r$long[[1]]))
  ))

  acq = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, presentation), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) <= 5)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) <= 5))
  ))

  early = ezANOVA(data = df_long, dv = V, wid = subject,
                  within = c(cue_type, presentation), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) > 5)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) > 5))
  ))

  late = ezANOVA(data = df_long, dv = V, wid = subject,
                 within = c(cue_type, presentation), between = group)

  #comparing reaction to the first vs second extinguished cue A presentation
  df_1v2 <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) <= 2)),
    lapply(ptsd, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) <= 2))
  ))

  ctrl_1v2 <- t.test(
    x = subset(df_1v2, group == "ctrl" & as.integer(presentation) == 1)$V,
    y = subset(df_1v2, group == "ctrl" & as.integer(presentation) == 2)$V,
    paired = TRUE
  )

  ptsd_1v2 <- t.test(
    x = subset(df_1v2, group == "ptsd" & as.integer(presentation) == 1)$V,
    y = subset(df_1v2, group == "ptsd" & as.integer(presentation) == 2)$V,
    paired = TRUE
  )

  #comparing reaction to the second vs third extinguished cue A presentation
  df_2v3 <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) %in% c(2, 3))),
    lapply(ptsd, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) %in% c(2, 3)))
  ))

  ctrl_2v3 <- t.test(
    x = subset(df_2v3, group == "ctrl" & as.integer(presentation) == 2)$V,
    y = subset(df_2v3, group == "ctrl" & as.integer(presentation) == 3)$V,
    paired = TRUE
  )

  ptsd_2v3 <- t.test(
    x = subset(df_2v3, group == "ptsd" & as.integer(presentation) == 2)$V,
    y = subset(df_2v3, group == "ptsd" & as.integer(presentation) == 3)$V,
    paired = TRUE
  )

  c(
    early_extinction_ANOVA_cue_type_main_effect = check_p_local(early, "cue_type", TRUE),
    early_extinction_ANOVA_trial_number_main_effect = check_p_local(early, "presentation", TRUE),
    early_extinction_ANOVA_group_x_cue_type_interaction = check_p_local(early, "group:cue_type", TRUE),
    late_extinction_ANOVA_trial_number_main_effect = check_p_local(late, "presentation", TRUE),
    control_difference_between_trials_1_and_2_of_extinction = check_p_other_local(ctrl_1v2$p.value, TRUE),
    PTSD_difference_between_trials_2_and_3_of_extinction = check_p_other_local(ptsd_2v3$p.value, TRUE),
    acquisition_ANOVA_group_main_effect = check_p_local(acq, "group", FALSE),
    group_extinction_group_main_effect = check_p_local(early, "group", FALSE),
    late_extinction_group_main_effect = check_p_local(late, "group", FALSE),
    late_extinction_cue_type_main_effect = check_p_local(late, "cue_type", FALSE),
    late_extinction_group_x_cue_type_interaction = check_p_local(late, "group:cue_type", FALSE),
    PTSD_difference_between_trials_1_and_2_of_extinction = check_p_other_local(ptsd_1v2$p.value, FALSE),
    control_difference_between_trials_2_and_3_of_extinction = check_p_other_local(ctrl_2v3$p.value, FALSE)
  )
}

zuj2017 <- run_power_analysis(
  study_fun  = study_zuj2017,
  n_ctrl = 33, n_ptsd = 21,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "zuj2017",
  n_ctrl = 33,
  n_ptsd = 21,
  Marginal = names(zuj2017$marginal),
  Value = unname(zuj2017$marginal),
  Joint = zuj2017$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("zuj2017_prob_", noise, ".csv"))

max_passes <- length(zuj2017$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- zuj2017$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "zuj2017",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("zuj2017_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "zuj2017",
  Pass_pattern = names(zuj2017$pass_patterns),
  n_passes = sapply(
    strsplit(names(zuj2017$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(zuj2017$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("zuj2017_pattern_", noise, ".csv"))

###########################
###Garfinkel et al. 2014###
###########################

study_garfinkel2014 <- function(seed, n_ctrl, n_ptsd, noise,
                          eta, eta_ptsd, lambda,
                          alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        AX = 8,
        BX = 8,
        CX = 16
      ),
      rewards = list(
        c(AX = 5/8, BX = 5/8)
      )
    ),

    phase2 = list(
      trials = list(
        AY = 16,
        CY = 16
      )
    ),

    day2_recall = list(
      trials = list(
        AY = 8,
        BY = 8,
        CY = 16
      )
    ),

    day2_renewal = list(
      trials = list(
        AX = 8,
        BX = 8,
        CX = 16
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, features = c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), pseudo=3, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, features = c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), pseudo=3, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) > 8)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) > 8))
  ))

  ext = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) rbind(
      subset(r$long[[3]], cue_type == "AY" & as.integer(presentation) <= 4),
      subset(r$long[[3]], cue_type == "CY" & as.integer(presentation) <= 8)
    )),
    lapply(ptsd, function(r) rbind(
      subset(r$long[[3]], cue_type == "AY" & as.integer(presentation) <= 4),
      subset(r$long[[3]], cue_type == "CY" & as.integer(presentation) <= 8)
    ))
  ))

  rec = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) rbind(
      subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4),
      subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)
    )),
    lapply(ptsd, function(r) rbind(
      subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4),
      subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)
    ))
  ))

  ren = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #compute renewal magnitude (AX - CX) per subject using day2_renewal (phase 4)
  renewal_ctrl <- sapply(ctrl, function(r) {
    ax <- mean(subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4)$V)
    cx <- mean(subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)$V)
    ax - cx
  })

  renewal_ptsd <- sapply(ptsd, function(r) {
    ax <- mean(subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4)$V)
    cx <- mean(subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)$V)
    ax - cx
  })

  t_renewal <- t.test(
    x = renewal_ctrl,
    y = renewal_ptsd,
    paired = FALSE
  )

  #t test for cue discrimination in extinction
  ctrl_ay_ext <- sapply(ctrl, function(r) {
    mean(subset(r$long[[2]], cue_type == "AY" & as.integer(presentation) > 8)$V)
  })

  ctrl_cy_ext <- sapply(ctrl, function(r) {
    mean(subset(r$long[[2]], cue_type == "CY" & as.integer(presentation) > 8)$V)
  })

  t_ctrl_ext <- t.test(
    x = ctrl_ay_ext,
    y = ctrl_cy_ext,
    paired = TRUE
  )

  ptsd_ay_ext <- sapply(ptsd, function(r) {
    mean(subset(r$long[[2]], cue_type == "AY" & as.integer(presentation) > 8)$V)
  })

  ptsd_cy_ext <- sapply(ptsd, function(r) {
    mean(subset(r$long[[2]], cue_type == "CY" & as.integer(presentation) > 8)$V)
  })

  t_ptsd_ext <- t.test(
    x = ptsd_ay_ext,
    y = ptsd_cy_ext,
    paired = TRUE
  )

  #t test for cue discrimination in recall
  ctrl_ay_rec <- sapply(ctrl, function(r) {
    mean(subset(r$long[[3]], cue_type == "AY" & as.integer(presentation) <= 4)$V)
  })

  ctrl_cy_rec <- sapply(ctrl, function(r) {
    mean(subset(r$long[[3]], cue_type == "CY" & as.integer(presentation) <= 8)$V)
  })

  t_ctrl_rec <- t.test(
    x = ctrl_ay_rec,
    y = ctrl_cy_rec,
    paired = TRUE
  )

  ptsd_ay_rec <- sapply(ptsd, function(r) {
    mean(subset(r$long[[3]], cue_type == "AY" & as.integer(presentation) <= 4)$V)
  })

  ptsd_cy_rec <- sapply(ptsd, function(r) {
    mean(subset(r$long[[3]], cue_type == "CY" & as.integer(presentation) <= 8)$V)
  })

  t_ptsd_rec <- t.test(
    x = ptsd_ay_rec,
    y = ptsd_cy_rec,
    paired = TRUE
  )

  #t test for cue discrimination in renewal
  ctrl_ax_ren <- sapply(ctrl, function(r) {
    mean(subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4)$V)
  })

  ctrl_cx_ren <- sapply(ctrl, function(r) {
    mean(subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)$V)
  })

  t_ctrl_ren <- t.test(
    x = ctrl_ax_ren,
    y = ctrl_cx_ren,
    paired = TRUE
  )

  ptsd_ax_ren <- sapply(ptsd, function(r) {
    mean(subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4)$V)
  })

  ptsd_cx_ren <- sapply(ptsd, function(r) {
    mean(subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)$V)
  })

  t_ptsd_ren <- t.test(
    x = ptsd_ax_ren,
    y = ptsd_cx_ren,
    paired = TRUE
  )

  c(
    extinction_ANOVA_group_main_effect = check_p_local(ext, "group", FALSE),
    control_t_test_between_stimuli_in_extinction = check_p_other_local(t_ctrl_ext$p.value, FALSE),
    PTSD_t_test_between_stimuli_in_extinction = check_p_other_local(t_ptsd_ext$p.value, FALSE),

    recall_ANOVA_group_x_cue_type_interaction = check_p_local(rec, "group:cue_type", TRUE),
    renewal_ANOVA_group_x_cue_type_interaction = check_p_local(ren, "group:cue_type", TRUE),
    PTSD_t_test_between_stimuli_in_recall = check_p_other_local(t_ptsd_rec$p.value, TRUE),
    control_t_test_between_stimuli_in_renewal = check_p_other_local(t_ctrl_ren$p.value, TRUE),
    t_test_between_groups_in_renewal = check_p_other_local(t_renewal$p.value, TRUE),
    control_t_test_between_stimuli_in_recall = check_p_other_local(t_ctrl_rec$p.value, FALSE),
    PTSD_t_test_between_stimuli_in_renewal = check_p_other_local(t_ptsd_ren$p.value, FALSE)
  )
}

garfinkel2014 <- run_power_analysis(
  study_fun  = study_garfinkel2014,
  n_ctrl = 14, n_ptsd = 14,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = c(
    "recall_ANOVA_group_x_cue_type_interaction",
    "renewal_ANOVA_group_x_cue_type_interaction",
    "PTSD_t_test_between_stimuli_in_recall",
    "control_t_test_between_stimuli_in_renewal",
    "t_test_between_groups_in_renewal",
    "control_t_test_between_stimuli_in_recall",
    "PTSD_t_test_between_stimuli_in_renewal"
  )
)

df_marginal <- data.frame(
  Name = "garfinkel2014",
  n_ctrl = 14,
  n_ptsd = 14,
  Marginal = names(garfinkel2014$marginal),
  Value = unname(garfinkel2014$marginal),
  Joint = garfinkel2014$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("garfinkel2014_prob_", noise, ".csv"))

max_passes <- length(garfinkel2014$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- garfinkel2014$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "garfinkel2014",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("garfinkel2014_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "garfinkel2014",
  Pass_pattern = names(garfinkel2014$pass_patterns),
  n_passes = sapply(
    strsplit(names(garfinkel2014$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(garfinkel2014$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("garfinkel2014_pattern_", noise, ".csv"))

###########################
###Jovanovic et al. 2013###
###########################

study_jovanovic2013 <- function(seed, n_ctrl, n_ptsd, noise,
                                eta, eta_ptsd, lambda,
                                alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 24,
        B = 24
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=8, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #ANOVA comparing response to cue A vs cue B separately for each group
  df_ctrl_acq <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]]))
  ))

  ctrl_acq = ezANOVA(data = df_ctrl_acq, dv = V, wid = subject,
                     within = c(cue_type))

  df_ptsd_acq <- do.call(rbind, c(
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  ptsd_acq = ezANOVA(data = df_ptsd_acq, dv = V, wid = subject,
                     within = c(cue_type))

  #ANOVAs checking whether response to cue A changed with block
  df_ctrl_ext <- do.call(rbind,
                         lapply(ctrl, function(r) subset(r$long[[2]], cue_type == "A"))
  )

  ext_ctrl <- ezANOVA(data = df_ctrl_ext, dv = V, wid = subject,
                      within = block)

  df_ptsd_ext <- do.call(rbind,
                         lapply(ptsd, function(r) subset(r$long[[2]], cue_type == "A"))
  )

  ext_ptsd <- ezANOVA(data = df_ptsd_ext, dv = V, wid = subject,
                      within = block)

  c(
    control_acquistion_cue_type_main_effect = check_p_local(ctrl_acq, "cue_type", TRUE),
    control_extinction_cue_type_main_effect = check_p_local(ext_ctrl, "block", TRUE),
    PTSD_acquistion_cue_type_main_effect = check_p_local(ptsd_acq, "cue_type", FALSE),
    PTSD_extinction_cue_type_main_effect = check_p_local(ext_ptsd, "block", FALSE)

  )
}

jovanovic2013 <- run_power_analysis(
  study_fun  = study_jovanovic2013,
  n_ctrl = 12, n_ptsd = 12,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "jovanovic2013",
  n_ctrl = 12,
  n_ptsd = 12,
  Marginal = names(jovanovic2013$marginal),
  Value = unname(jovanovic2013$marginal),
  Joint = jovanovic2013$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("jovanovic2013_prob_", noise, ".csv"))

max_passes <- length(jovanovic2013$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- jovanovic2013$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "jovanovic2013",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("jovanovic2013_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "jovanovic2013",
  Pass_pattern = names(jovanovic2013$pass_patterns),
  n_passes = sapply(
    strsplit(names(jovanovic2013$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(jovanovic2013$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("jovanovic2013_pattern_", noise, ".csv"))

#########################
###Wicking et al. 2016###
#########################

study_wicking2016 <- function(seed, n_ctrl, n_ptsd, noise,
                                eta, eta_ptsd, lambda,
                                alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  #X is context 1, Y is context 2, Z is context 3
  phase_def <- list(

    phase1 = list(
      trials = list(
        AX = 30,
        BX = 30
      ),
      rewards = list(
        c(AX = 1)
      )
    ),

    phase2_day2 = list(
      trials = list(
        AY = 30,
        BY = 30
      )
    ),

    renewal_day9 = list(
      trials = list(
        AZ = 30,
        BZ = 30
      )
    )
  )

  #no need to put ezITI between the first and second phase since there is only
  #one latent state after the first phase
  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, pseudo=3, features=c("A","B","X","Y","Z"), ezITI=list(c(120,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, pseudo=3, features=c("A","B","X","Y","Z"), ezITI=list(c(120,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #second half of acquisition
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]], as.integer(presentation) > 15)),
    lapply(ptsd, function(r) subset(r$long[[1]], as.integer(presentation) > 15))
  ))

  acq = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #second half of extinction
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) > 15)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) > 15))
  ))

  ext = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #first half of renewal
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[3]], as.integer(presentation) <= 15)),
    lapply(ptsd, function(r) subset(r$long[[3]], as.integer(presentation) <= 15))
  ))

  ren = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #the t tests use D-score, which is mean CS+ minus mean CS- for the renewal phase
  df_dscore <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop") %>%
    pivot_wider(names_from = cue_type, values_from = V) %>%
    mutate(Dscore = AZ - BZ)

  #within-group one-sample t-tests on D-score vs 0
  ren_ttest_ctrl <- t.test(df_dscore$Dscore[df_dscore$group == "ctrl"],
                           mu = 0, alternative = "greater")
  ren_ttest_ptsd <- t.test(df_dscore$Dscore[df_dscore$group == "ptsd"],
                           mu = 0, alternative = "greater")

  c(
    acquisition_ANOVA_group_x_cue_type_interaction = check_p_local(acq, "group:cue_type", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext, "group:cue_type", FALSE),

    recall_ANOVA_group_x_cue_type_interaction = check_p_local(ren, "group:cue_type", TRUE),
    PTSD_t_test_between_stimuli_in_recall = check_p_other_local(ren_ttest_ptsd$p.value, TRUE),
    control_t_test_between_stimuli_in_recall = check_p_other_local(ren_ttest_ctrl$p.value, FALSE)
  )
}

wicking2016 <- run_power_analysis(
  study_fun  = study_wicking2016,
  n_ctrl = 18, n_ptsd = 36,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = c("recall_ANOVA_group_x_cue_type_interaction", "PTSD_t_test_between_stimuli_in_recall", "control_t_test_between_stimuli_in_recall")
)

df_marginal <- data.frame(
  Name = "wicking2016",
  n_ctrl = 18,
  n_ptsd = 36,
  Marginal = names(wicking2016$marginal),
  Value = unname(wicking2016$marginal),
  Joint = wicking2016$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("wicking2016_prob_", noise, ".csv"))

max_passes <- length(wicking2016$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- wicking2016$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "wicking2016",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("wicking2016_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "wicking2016",
  Pass_pattern = names(wicking2016$pass_patterns),
  n_passes = sapply(
    strsplit(names(wicking2016$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(wicking2016$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("wicking2016_pattern_", noise, ".csv"))

################################Expectancy######################################

##########################
###Blechert et al. 2007###
##########################

study_blechert2007_x <- function(seed, n_ctrl, n_ptsd, noise,
                               eta, eta_ptsd, lambda,
                               alpha0, alpha1, alpha2, gamma, delta, sigma0) {

  phase_def <- list(
    phase1 = list(
      trials  = list(A = 6, B = 6),
      rewards = list(c(A = 1))
    ),
    phase2 = list(
      trials = list(A = 6, B = 6)
    )
  )

  ctrl <- run(phase_def = phase_def, seed = seed, V_noise_sd = noise,
              n = n_ctrl, block = 6, eta = eta, lambda = lambda,
              alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
              gamma = gamma, delta = delta, sigma0 = sigma0)

  ptsd <- run(phase_def = phase_def, seed = seed + 1e6, V_noise_sd = noise,
              n = n_ptsd, block = 6, eta = eta_ptsd, lambda = lambda,
              alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
              gamma = gamma, delta = delta, sigma0 = sigma0,
              prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) == 6)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) == 6))
  ))

  ext = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type), between = group)

  #post-hoc single-CS comparisons use simple between-group t-tests, since there
  #is no within-subject factor left once block is removed
  CSplus <- t.test(V ~ group, data = filter(df_long, cue_type == "A"))
  CSmin  <- t.test(V ~ group, data = filter(df_long, cue_type == "B"))

  c(
    extinction_ANOVA_group_main_effect   = check_p_local(ext, "group", TRUE),
    extinction_ANOVA_cue_type_main_effect       = check_p_local(ext, "cue_type", TRUE),
    CS_plus_group_difference  = check_p_other_local(CSplus$p.value, TRUE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext, "group:cue_type", FALSE),
    CS_minus_group_difference   = check_p_other_local(CSmin$p.value, FALSE)
  )
}

blechert2007_x <- run_power_analysis(
  study_fun  = study_blechert2007_x,
  n_ctrl = 34, n_ptsd = 36,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "blechert2007_x",
  n_ctrl = 34,
  n_ptsd = 36,
  Marginal = names(blechert2007_x$marginal),
  Value = unname(blechert2007_x$marginal),
  Joint = blechert2007_x$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("blechert2007_x_prob_", noise, ".csv"))

max_passes <- length(blechert2007_x$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- blechert2007_x$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "blechert2007_x",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("blechert2007_x_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "blechert2007_x",
  Pass_pattern = names(blechert2007_x$pass_patterns),
  n_passes = sapply(
    strsplit(names(blechert2007_x$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(blechert2007_x$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("blechert2007_x_pattern_", noise, ".csv"))

###########################
###Pohlchen et al., 2020###
###########################

study_pohlechen2020_x <- function(seed, n_ctrl, n_ptsd, noise,
                                eta, eta_ptsd, lambda,
                                alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12,
        C = 12
      ),
      rewards = list(
        c(A = 0.75, B = 0.75)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        C = 10
      )
    ),

    day2 = list(
      trials = list(
        A = 8,
        B = 8,
        C = 8
      )
    ),

    #this phase is not part of the experiment, and is not analyzed, except to allow
    #us to see the V of each cue at the very end of the experiment. Because V
    #stores the value of each cue BEFORE the outcome is shown, this bonus trial
    #actually represents the V at the end of the final block of day2, which is
    #needed for analyzing threat expectancy ratings

    bonus = list(
      trials = list(
        A = 1
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=c(12,10,12,1), ezITI=list(c(56,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=c(12,10,12,1), ezITI=list(c(56,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #Now we're going to do the threat expectancy analysis. Since V records
  #associative value on a given trial before the outcome is shown, we want to
  #index the trial directly after each block. This is why we created the bonus
  #phase. We do a separate ANOVA for each phase, using only the "end of block"
  #trials
  df_acq <- do.call(rbind, c(
    lapply(seq_along(ctrl), function(i) {
      r <- ctrl[[i]]
      rbind(
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "A", block = factor(1:3), V = r$V[["A"]][c(13, 25, 37), 1]),
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "C", block = factor(1:3), V = r$V[["C"]][c(13, 25, 37), 1])
      )
    }),
    lapply(seq_along(ptsd), function(i) {
      r <- ptsd[[i]]
      rbind(
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "A", block = factor(1:3), V = r$V[["A"]][c(13, 25, 37), 1]),
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "C", block = factor(1:3), V = r$V[["C"]][c(13, 25, 37), 1])
      )
    })
  ))


  acq_x <- ezANOVA(data = df_acq, dv = V, wid = subject,
                   within = c(cue_type, block), between = group)

  df_rec <- do.call(rbind, c(
    lapply(seq_along(ctrl), function(i) {
      r <- ctrl[[i]]
      rbind(
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "A", block = factor(1:2), V = r$V[["A"]][c(69, 81), 1]),
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "C", block = factor(1:2), V = r$V[["C"]][c(69, 81), 1])
      )
    }),
    lapply(seq_along(ptsd), function(i) {
      r <- ptsd[[i]]
      rbind(
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "A", block = factor(1:2), V = r$V[["A"]][c(69, 81), 1]),
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "C", block = factor(1:2), V = r$V[["C"]][c(69, 81), 1])
      )
    })
  ))


  rec_x <- ezANOVA(data = df_rec, dv = V, wid = subject,
                   within = c(cue_type, block), between = group)

  #Pohlechen also looked at threat expectancy, which involved rating each cue
  #after each block. Since V records the associative value on a given trial
  #before the outcome is shown, we want to index the trial directly after
  #each block. However, ezITI means that beliefs normalize on trial 57 BEFORE
  #the associative value is recorded, so V of trial 57 with the "true" experiment
  #setup is not actually the V directly after trial 56. We therefore make a new
  #version of the experiment that, instead of ezITI, has one last phase 2 trial
  #whose outcome is irrelevant, so that we can record V directly after trial 56

  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12,
        C = 12
      ),
      rewards = list(
        c(A = 0.75, B = 0.75)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        C = 10
      )
    ),

    #this bonus phase is just so we can index the V of the trial directly after
    #phase 2, which represents the V at the end of phase 2
    bonus = list(
      trials = list(
        A = 1
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=c(12,10,1), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=c(12,10,1), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #threat expectancy analysis using the end of each block (i.e. the V of the
  #trial directly after each block) for the extinction phase
  df_ext <- do.call(rbind, c(
    lapply(seq_along(ctrl), function(i) {
      r <- ctrl[[i]]
      rbind(
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "A", block = factor(1:2), V = r$V[["A"]][c(47, 57), 1]),
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "C", block = factor(1:2), V = r$V[["C"]][c(47, 57), 1])
      )
    }),
    lapply(seq_along(ptsd), function(i) {
      r <- ptsd[[i]]
      rbind(
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "A", block = factor(1:2), V = r$V[["A"]][c(47, 57), 1]),
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "C", block = factor(1:2), V = r$V[["C"]][c(47, 57), 1])
      )
    })
  ))


  ext_x <- ezANOVA(data = df_ext, dv = V, wid = subject,
                   within = c(cue_type, block), between = group)

  c(
    acquisition_ANOVA_group_main_effect = check_p_local(acq_x, "group", FALSE),
    extinction_ANOVA_group_main_effect = check_p_local(ext_x, "group", FALSE),
    acquisition_ANOVA_group_x_cue_type_interaction = check_p_local(acq_x, "group:cue_type", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext_x, "group:cue_type", FALSE),

    recall_ANOVA_group_main_effect = check_p_local(rec_x, "group", FALSE),
    recall_ANOVA_group_x_cue_type_interaction = check_p_local(rec_x, "group:cue_type", FALSE)
  )
}

pohlechen2020_x <- run_power_analysis(
  study_fun  = study_pohlechen2020_x,
  n_ctrl = 35, n_ptsd = 21,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = c("recall_ANOVA_group_main_effect", "recall_ANOVA_group_x_cue_type_interaction")
)

df_marginal <- data.frame(
  Name = "pohlechen2020_x",
  n_ctrl = 35,
  n_ptsd = 21,
  Marginal = names(pohlechen2020_x$marginal),
  Value = unname(pohlechen2020_x$marginal),
  Joint = pohlechen2020_x$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("pohlechen2020_x_prob_", noise, ".csv"))

max_passes <- length(pohlechen2020_x$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- pohlechen2020_x$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "pohlechen2020_x",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("pohlechen2020_x_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "pohlechen2020_x",
  Pass_pattern = names(pohlechen2020_x$pass_patterns),
  n_passes = sapply(
    strsplit(names(pohlechen2020_x$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(pohlechen2020_x$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("pohlechen2020_x_pattern_", noise, ".csv"))

##########################
###Norrholm et al. 2011###
##########################

study_norrholm2011_x <- function(seed, n_ctrl, n_ptsd, noise,
                                  eta, eta_ptsd, lambda,
                                  alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 12,
        B = 12
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 24,
        B = 24
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=8, eta=eta_ptsd, alpha1=alpha1, lambda=lambda, alpha0=alpha0, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #for expectancy ratings, participants could only rate "danger", "safety", and
  #"uncertain". Therefore, we divide V into three categories of equal size
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]], as.integer(presentation) > 8)),
    lapply(ptsd, function(r) subset(r$long[[1]], as.integer(presentation) > 8))
  ))

  df_long$expectancy <- cut(
    df_long$V,
    breaks = c(-0.5, -1/6, 1/6, 0.5),
    include.lowest = TRUE,
    labels = c(-1, 0, 1)
  )

  df_long$expectancy <- as.numeric(as.character(df_long$expectancy))

  df_long$subject  <- factor(df_long$subject)
  df_long$group    <- factor(df_long$group)
  df_long$cue_type <- factor(df_long$cue_type)

  df_acq_avg <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(
      expectancy = mean(expectancy, na.rm = TRUE),
      .groups = "drop"
    )

  acq_x <- tryCatch(
    ezANOVA(data = df_acq_avg, dv = expectancy, wid = subject,
            within = c(cue_type), between = group),
    error = function(e) NULL
  )

  #same thing for extinction
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  df_long$expectancy <- cut(
    df_long$V,
    breaks = c(-0.5, -1/6, 1/6, 0.5),
    include.lowest = TRUE,
    labels = c(-1, 0, 1)
  )

  df_long$expectancy <- as.numeric(as.character(df_long$expectancy))

  df_long$subject  <- factor(df_long$subject)
  df_long$group    <- factor(df_long$group)
  df_long$cue_type <- factor(df_long$cue_type)
  df_long$superblock <- factor(ceiling(as.integer(df_long$block) / 2))

  ext_x <- tryCatch(
    ezANOVA(
      data = df_long,
      dv = expectancy,
      wid = subject,
      within = .(cue_type, superblock),
      between = .(group)
    ),
    error = function(e) NULL
  )

  c(
    acqustion_ANOVA_group_main_effect = check_p_local(acq_x, "group", FALSE),
    acqustion_ANOVA_group_x_cue_type_interaction = check_p_local(acq_x, "group:cue_type", FALSE),
    extinction_ANOVA_group_main_effect = check_p_local(ext_x, "group", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext_x, "group:cue_type", FALSE)
  )
}

norrholm2011_x <- run_power_analysis(
  study_fun  = study_norrholm2011_x,
  n_ctrl = 78, n_ptsd = 49,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "norrholm2011_x",
  n_ctrl = 78,
  n_ptsd = 49,
  Marginal = names(norrholm2011_x$marginal),
  Value = unname(norrholm2011_x$marginal),
  Joint = norrholm2011_x$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("norrholm2011_x_prob_", noise, ".csv"))

max_passes <- length(norrholm2011_x$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- norrholm2011_x$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "norrholm2011_x",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("norrholm2011_x_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "norrholm2011_x",
  Pass_pattern = names(norrholm2011_x$pass_patterns),
  n_passes = sapply(
    strsplit(names(norrholm2011_x$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(norrholm2011_x$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("norrholm2011_x_pattern_", noise, ".csv"))

##########################
###Acheson et al. 2015####
##########################

study_acheson2015_x <- function(seed, n_ctrl, n_ptsd, noise,
                                 eta, eta_ptsd, lambda,
                                 alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 8,
        B = 8
      ),
      rewards = list(
        c(A = 0.75)
      )
    ),

    phase2 = list(
      trials = list(
        A = 16,
        B = 16
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=8, eta=eta, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, lambda=lambda)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=8, eta=eta_ptsd, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, lambda=lambda, prev=n_ctrl, group="ptsd")


  #for expectancy ratings, participants could only rate "expect", "don't expect",
  #and "unsure". Therefore, we divide V into three categories of equal size
  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]], as.integer(presentation) > 6)),
    lapply(ptsd, function(r) subset(r$long[[1]], as.integer(presentation) > 6))
  ))

  df_avg_last_two_acq <- df_long %>%
    group_by(subject, group, cue_type) %>%
    summarise(V = mean(V), .groups = "drop")

  df_avg_last_two_acq$expectancy <- cut(
    df_avg_last_two_acq$V,
    breaks = c(-0.5, -1/6, 1/6, 0.5),
    include.lowest = TRUE,
    labels = c(-1, 0, 1)
  )

  df_avg_last_two_acq$expectancy <- as.numeric(as.character(df_avg_last_two_acq$expectancy))

  df_avg_last_two_acq$subject  <- factor(df_avg_last_two_acq$subject)
  df_avg_last_two_acq$group    <- factor(df_avg_last_two_acq$group)
  df_avg_last_two_acq$cue_type <- factor(df_avg_last_two_acq$cue_type)

  #because there are only 3 categories which is very coarse, in some analyses
  #there will be no differentiation between groups and the ANOVA will fail,
  #we can be considered a non-significant difference
  acq_x <- tryCatch(
    ezANOVA(data = df_avg_last_two_acq, dv = expectancy, wid = subject,
            within = c(cue_type), between = group),
    error = function(e) NULL
  )

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]])),
    lapply(ptsd, function(r) subset(r$long[[2]]))
  ))

  #no averaging for extinction expectancy ratings, though only cue A is considered
  df_ext <- df_long %>%
    filter(cue_type == "A")

  df_ext$expectancy <- cut(
    df_ext$V,
    breaks = c(-0.5, -1/6, 1/6, 0.5),
    include.lowest = TRUE,
    labels = c(-1, 0, 1)
  )

  df_ext$expectancy <- as.numeric(as.character(df_ext$expectancy))

  df_ext$subject  <- factor(df_ext$subject)
  df_ext$group    <- factor(df_ext$group)
  df_ext$cue_type <- factor(df_ext$cue_type)

  ext_x <- tryCatch(
    ezANOVA(data = df_ext, dv = expectancy, wid = subject,
            within = c(presentation), between = group),
    error = function(e) NULL
  )

  c(
    acquisition_ANOVA_group_main_effect = check_p_local(acq_x, "group", FALSE),
    acquisition_ANOVA_group_x_cue_type_interaction = check_p_local(acq_x, "group:cue_type", FALSE),
    extinction_ANOVA_group_main_effect = check_p_local(ext_x, "group", FALSE),
    extinction_ANOVA_group_x_trial_number_interaction = check_p_local(ext_x, "group:presentation", FALSE)
  )
}

acheson2015_x <- run_power_analysis(
  study_fun  = study_acheson2015_x,
  n_ctrl = 923, n_ptsd = 42,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "acheson2015_x",
  n_ctrl = 923,
  n_ptsd = 42,
  Marginal = names(acheson2015_x$marginal),
  Value = unname(acheson2015_x$marginal),
  Joint = acheson2015_x$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("acheson2015_x_prob_", noise, ".csv"))

max_passes <- length(acheson2015_x$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- acheson2015_x$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "acheson2015_x",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("acheson2015_x_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "acheson2015_x",
  Pass_pattern = names(acheson2015_x$pass_patterns),
  n_passes = sapply(
    strsplit(names(acheson2015_x$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(acheson2015_x$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("acheson2015_x_pattern_", noise, ".csv"))

#####################
###Zuj et al. 2017###
#####################

study_zuj2017_x <- function(seed, n_ctrl, n_ptsd, noise,
                          eta, eta_ptsd, lambda,
                          alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  phase_def <- list(

    phase1 = list(
      trials = list(
        A = 5,
        B = 5
      ),
      rewards = list(
        c(A = 1)
      )
    ),

    phase2 = list(
      trials = list(
        A = 10,
        B = 10
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, block=10, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, block=10, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[1]])),
    lapply(ptsd, function(r) subset(r$long[[1]]))
  ))

  acq = ezANOVA(data = df_long, dv = V, wid = subject,
                within = c(cue_type, presentation), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) <= 5)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) <= 5))
  ))

  early = ezANOVA(data = df_long, dv = V, wid = subject,
                  within = c(cue_type, presentation), between = group)

  df_long <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], as.integer(presentation) > 5)),
    lapply(ptsd, function(r) subset(r$long[[2]], as.integer(presentation) > 5))
  ))

  late = ezANOVA(data = df_long, dv = V, wid = subject,
                 within = c(cue_type, presentation), between = group)

  #comparing reaction to the first vs second extinguished cue A presentation
  df_1v2 <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) <= 2)),
    lapply(ptsd, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) <= 2))
  ))

  ctrl_1v2 <- t.test(
    x = subset(df_1v2, group == "ctrl" & as.integer(presentation) == 1)$V,
    y = subset(df_1v2, group == "ctrl" & as.integer(presentation) == 2)$V,
    paired = TRUE
  )

  ptsd_1v2 <- t.test(
    x = subset(df_1v2, group == "ptsd" & as.integer(presentation) == 1)$V,
    y = subset(df_1v2, group == "ptsd" & as.integer(presentation) == 2)$V,
    paired = TRUE
  )

  #comparing reaction to the second vs third extinguished cue A presentation
  df_2v3 <- do.call(rbind, c(
    lapply(ctrl, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) %in% c(2, 3))),
    lapply(ptsd, function(r) subset(r$long[[2]], cue_type == "A" & as.integer(presentation) %in% c(2, 3)))
  ))

  ctrl_2v3 <- t.test(
    x = subset(df_2v3, group == "ctrl" & as.integer(presentation) == 2)$V,
    y = subset(df_2v3, group == "ctrl" & as.integer(presentation) == 3)$V,
    paired = TRUE
  )

  ptsd_2v3 <- t.test(
    x = subset(df_2v3, group == "ptsd" & as.integer(presentation) == 2)$V,
    y = subset(df_2v3, group == "ptsd" & as.integer(presentation) == 3)$V,
    paired = TRUE
  )

  c(
    early_extinction_cue_type_x_trial_number_interaction = check_p_local(early, "cue_type:presentation", TRUE),
    late_extinction_group_main_effect = check_p_local(late, "group", TRUE),
    late_extinction_cue_type_main_effect = check_p_local(late, "cue_type", TRUE),
    late_extinction_trial_number_main_effect = check_p_local(late, "presentation", TRUE),
    early_extinction_group_x_cue_type_interaction = check_p_local(early, "group:cue_type", FALSE),
    late_extinction_group_x_cue_type_interaction = check_p_local(late, "group:cue_type", FALSE)
  )
}

zuj2017_x <- run_power_analysis(
  study_fun  = study_zuj2017_x,
  n_ctrl = 33, n_ptsd = 21,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0
)

df_marginal <- data.frame(
  Name = "zuj2017_x",
  n_ctrl = 33,
  n_ptsd = 21,
  Marginal = names(zuj2017_x$marginal),
  Value = unname(zuj2017_x$marginal),
  Joint = zuj2017_x$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("zuj2017_x_prob_", noise, ".csv"))

max_passes <- length(zuj2017_x$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- zuj2017_x$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "zuj2017_x",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("zuj2017_x_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "zuj2017_x",
  Pass_pattern = names(zuj2017_x$pass_patterns),
  n_passes = sapply(
    strsplit(names(zuj2017_x$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(zuj2017_x$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("zuj2017_x_pattern_", noise, ".csv"))

#########################
###Wicking et al. 2016###
#########################

study_wicking2016_x <- function(seed, n_ctrl, n_ptsd, noise,
                            eta, eta_ptsd, lambda,
                            alpha0, alpha1, alpha2, gamma, delta, sigma0) {
  #X is context 1, Y is context 2, Z is context 3
  phase_def <- list(

    phase1 = list(
      trials = list(
        AX = 30,
        BX = 30
      ),
      rewards = list(
        c(AX = 1)
      )
    ),

    phase2_day2 = list(
      trials = list(
        AY = 30,
        BY = 30
      )
    ),

    renewal_day9 = list(
      trials = list(
        AZ = 30,
        BZ = 30
      )
    )
  )

  #no need to put ezITI between the first and second phase since there is only
  #one latent state after the first phase
  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, pseudo=3, features=c("A","B","X","Y","Z"), ezITI=list(c(120,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, pseudo=3, features=c("A","B","X","Y","Z"), ezITI=list(c(120,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  #Now we're going to do the threat expectancy analysis. Since V records
  #associative value on a given trial before the outcome is shown, we want to
  #index the trial directly after each block. This is why we created the bonus
  #phase. We do a separate ANOVA for each phase, using only the "end of block"
  #trials
  df_long <- do.call(rbind, c(
    lapply(seq_along(ctrl), function(i) {
      r <- ctrl[[i]]
      rbind(
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "AX",
                   V = r$V[["A"]][61, 1] + r$V[["X"]][61, 1]),
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "BX",
                   V = r$V[["B"]][61, 1] + r$V[["X"]][61, 1])
      )
    }),
    lapply(seq_along(ptsd), function(i) {
      r <- ptsd[[i]]
      rbind(
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "AX",
                   V = r$V[["A"]][61, 1] + r$V[["X"]][61, 1]),
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "BX",
                   V = r$V[["B"]][61, 1] + r$V[["X"]][61, 1])
      )
    })
  ))

  #for expectancy ratings, participants could only rate on a 9 point Likert scale.
  #Therefore, we divide V into nine categories of equal size
  #the bounds are extra wide in case there's overshoot
  df_long$expectancy <- cut(
    df_long$V,
    breaks = c(-1, -7/18, -5/18, -3/18, -1/18, 1/18, 3/18, 5/18, 7/18, 1),
    include.lowest = TRUE,
    labels = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  )

  df_long$expectancy <- as.numeric(as.character(df_long$expectancy))

  df_long$subject  <- factor(df_long$subject)
  df_long$group    <- factor(df_long$group)
  df_long$cue_type <- factor(df_long$cue_type)

  acq_x = tryCatch(ezANOVA(data = df_long, dv = expectancy, wid = subject,
                           within = c(cue_type), between = group),
                   error = function(e) NULL
  )


  #now we do the same for renewal
  df_long <- do.call(rbind, c(
    lapply(seq_along(ctrl), function(i) {
      r <- ctrl[[i]]
      rbind(
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "AZ",
                   V = r$V[["A"]][151, 1] + r$V[["Z"]][151, 1]),
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "BZ",
                   V = r$V[["B"]][151, 1] + r$V[["Z"]][151, 1])
      )
    }),
    lapply(seq_along(ptsd), function(i) {
      r <- ptsd[[i]]
      rbind(
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "AZ",
                   V = r$V[["A"]][151, 1] + r$V[["Z"]][151, 1]),
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "BZ",
                   V = r$V[["B"]][151, 1] + r$V[["Z"]][151, 1])
      )
    })
  ))

  #the bounds are extra wide in case there's overshoot
  df_long$expectancy <- cut(
    df_long$V,
    breaks = c(-1, -7/18, -5/18, -3/18, -1/18, 1/18, 3/18, 5/18, 7/18, 1),
    include.lowest = TRUE,
    labels = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  )

  df_long$expectancy <- as.numeric(as.character(df_long$expectancy))

  df_long$subject  <- factor(df_long$subject)
  df_long$group    <- factor(df_long$group)
  df_long$cue_type <- factor(df_long$cue_type)

  ren_x = tryCatch(
    ezANOVA(
      data = df_long,
      dv = expectancy,
      wid = subject,
      within = c(cue_type),
      between = group
    ),
    error = function(e)
      NULL
  )

  #the t tests use D-score, which is CS+ minus mean CS- for the renewal phase
  df_dscore <- df_long %>%
    select(subject, group, cue_type, expectancy) %>%
    group_by(subject, group) %>%
    pivot_wider(names_from = cue_type, values_from = expectancy) %>%
    mutate(Dscore = AZ - BZ)

  #within-group one-sample t-tests on D-score vs 0
  ren_ttest_ctrl_x <- tryCatch(
    t.test(df_dscore$Dscore[df_dscore$group == "ctrl"],
           mu = 0, alternative = "greater"),
    error = function(e) list(p.value = NA)
  )

  ren_ttest_ptsd_x <- tryCatch(
    t.test(df_dscore$Dscore[df_dscore$group == "ptsd"],
           mu = 0, alternative = "greater"),
    error = function(e) list(p.value = NA)
  )


  #Wicking also looked at threat expectancy, which involved rating each cue
  #after each block. Since V records the associative value on a given trial
  #before the outcome is shown, we want to index the trial directly after
  #each block. However, ezITI means that beliefs normalize on trial 121 BEFORE
  #the associative value is recorded, so V of trial 121 with the "true" experiment
  #setup is not actually the V directly after trial 120 We therefore make a new
  #version of the experiment that, instead of ezITI, has one last phase 2 trial
  #whose outcome is irrelevant, so that we can record V directly after trial 56

  phase_def <- list(

    phase1 = list(
      trials = list(
        AX = 30,
        BX = 30
      ),
      rewards = list(
        c(AX = 1)
      )
    ),

    phase2_day2 = list(
      trials = list(
        AY = 30,
        BY = 30
      )
    ),

    #this bonus phase is just so we can index the V of the trial directly after
    #phase 2, which represents the V at the end of phase 2
    bonus = list(
      trials = list(
        AY = 1
      )
    )
  )

  ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n_ctrl, pseudo=3, features=c("A","B","X","Y"), ezITI=list(c(120,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
  ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n_ptsd, pseudo=3, features=c("A","B","X","Y"), ezITI=list(c(120,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n_ctrl, group = "ptsd")

  df_long <- do.call(rbind, c(
    lapply(seq_along(ctrl), function(i) {
      r <- ctrl[[i]]
      rbind(
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "AY",
                   V = r$V[["A"]][121, 1] + r$V[["Y"]][121, 1]),
        data.frame(subject = factor(i),     group = "ctrl", cue_type = "BY",
                   V = r$V[["B"]][121, 1] + r$V[["Y"]][121, 1])
      )
    }),
    lapply(seq_along(ptsd), function(i) {
      r <- ptsd[[i]]
      rbind(
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "AY",
                   V = r$V[["A"]][121, 1] + r$V[["Y"]][121, 1]),
        data.frame(subject = factor(i + n_ctrl), group = "ptsd", cue_type = "BY",
                   V = r$V[["B"]][121, 1] + r$V[["Y"]][121, 1])
      )
    })
  ))

  #the bounds are extra wide in case there's overshoot
  df_long$expectancy <- cut(
    df_long$V,
    breaks = c(-1, -7/18, -5/18, -3/18, -1/18, 1/18, 3/18, 5/18, 7/18, 1),
    include.lowest = TRUE,
    labels = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  )

  df_long$expectancy <- as.numeric(as.character(df_long$expectancy))

  df_long$subject  <- factor(df_long$subject)
  df_long$group    <- factor(df_long$group)
  df_long$cue_type <- factor(df_long$cue_type)

  ext_x = tryCatch(ezANOVA(data = df_long, dv = expectancy, wid = subject,
                           within = c(cue_type), between = group),
                   error = function(e) NULL
  )

  c(
    acquistion_ANOVA_group_x_cue_type_interaction = check_p_local(acq_x, "group:cue_type", FALSE),
    extinction_ANOVA_group_x_cue_type_interaction = check_p_local(ext_x, "group:cue_type", FALSE),

    recall_ANOVA_group_x_cue_type_interaction = check_p_local(ren_x, "group:cue_type", FALSE),
    PTSD_t_test_between_stimuli_in_recall = check_p_other_local(ren_ttest_ptsd_x$p.value, TRUE),
    control_t_test_between_stimuli_in_recall = check_p_other_local(ren_ttest_ctrl_x$p.value, FALSE)
  )

}

wicking2016_x <- run_power_analysis(
  study_fun  = study_wicking2016_x,
  n_ctrl = 18, n_ptsd = 36,
  noise = noise, eta = eta, eta_ptsd = eta_ptsd, lambda = lambda,
  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
  gamma = gamma, delta = delta, sigma0 = sigma0,
  exclude_from_joint = c("recall_ANOVA_group_x_cue_type_interaction", "PTSD_t_test_between_stimuli_in_recall", "control_t_test_between_stimuli_in_recall")
)

df_marginal <- data.frame(
  Name = "wicking2016_x",
  n_ctrl = 18,
  n_ptsd = 36,
  Marginal = names(wicking2016_x$marginal),
  Value = unname(wicking2016_x$marginal),
  Joint = wicking2016_x$joint,
  stringsAsFactors = FALSE
)

write.csv(df_marginal, paste0("wicking2016_x_prob_", noise, ".csv"))

max_passes <- length(wicking2016_x$marginal)
pass_counts <- integer(max_passes + 1)
pass_table <- wicking2016_x$n_passed_per_run
pass_counts[as.integer(names(pass_table)) + 1] <- as.vector(pass_table)

df_passes <- data.frame(
  Name = "wicking2016_x",
  Passes = 0:max_passes,
  Count = pass_counts,
  stringsAsFactors = FALSE
)

write.csv(df_passes, paste0("wicking2016_x_pass_", noise, ".csv"))

df_patterns <- data.frame(
  Name = "wicking2016_x",
  Pass_pattern = names(wicking2016_x$pass_patterns),
  n_passes = sapply(
    strsplit(names(wicking2016_x$pass_patterns), ";"),
    function(x) ifelse(x[1] == "", 0, length(x))
  ),
  Count = as.vector(wicking2016_x$pass_patterns),
  stringsAsFactors = FALSE
)

write.csv(df_patterns, paste0("wicking2016_x_pattern_", noise, ".csv"))
