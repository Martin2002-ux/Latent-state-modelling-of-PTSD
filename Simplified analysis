library(latentState)
library(afex)
library(dplyr)
library(emmeans)
library(ez)
library(tidyr)
library(truncnorm)

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

  return((p < 0.05) == sig)
}

check_p <- function(ez, effect, sig) {
  if (is.null(ez)) {
    pass <<- c(pass, FALSE == sig)
    p=NA
  } else{
    mauchly <- ez$`Mauchly's Test for Sphericity`

    p <- if (!is.null(mauchly) && effect %in% mauchly$Effect &&
             mauchly %>% filter(Effect == effect) %>% pull(`p<.05`) == "*") {
      ez$`Sphericity Corrections` %>% filter(Effect == effect) %>% pull(`p[GG]`)
    } else {
      ez$ANOVA %>% filter(Effect == effect) %>% pull(p)
    }
    if (is.na(p)) {
      pass <<- c(pass, FALSE == sig)
    } else{
      pass <<- c(pass, (p < 0.05) == sig)
    }
  }

  p_values <<- c(p_values, p)

  return(p)
}

#for non-ANOVA, this function will check if the significance of the test matches
#the intended significance (true or false)
check_p_other <- function(p, sig) {
  if (is.na(p)) {
    pass <<- c(pass, FALSE==sig)
  } else{
    pass <<- c(pass, (p < 0.05) == sig)
  }
  p_values <<- c(p_values, p)
  return(p)
}

#the check_p equivalent for non-physiological expectancy analyses
check_p_x <- function(ez, effect, sig) {
  if (is.null(ez)) {
    pass_x <<- c(pass_x, FALSE == sig)
    p=NA
  } else{
    mauchly <- ez$`Mauchly's Test for Sphericity`

    p <- if (!is.null(mauchly) && effect %in% mauchly$Effect &&
             mauchly %>% filter(Effect == effect) %>% pull(`p<.05`) == "*") {
      ez$`Sphericity Corrections` %>% filter(Effect == effect) %>% pull(`p[GG]`)
    } else {
      ez$ANOVA %>% filter(Effect == effect) %>% pull(p)
    }
    if (is.na(p)) {
      pass_x <<- c(pass_x, FALSE == sig)
    } else{
      pass_x <<- c(pass_x, (p < 0.05) == sig)
    }

  }
  p_values_x <<- c(p_values_x, p)

  return(p)
}

#the check_p_other equivalent for non-physiological expectancy analyses
check_p_other_x <- function(p, sig) {
  if (is.na(p)) {
    pass_x <<- c(pass_x, FALSE==sig)
  } else{
    pass_x <<- c(pass_x, (p < 0.05) == sig)
  }
  p_values_x <<- c(p_values_x, p)
  return(p)
}

#These are the original default values from the Cochran Cisler (2019) paper
alpha0 = "0.05"
alpha1 = "0.05"
alpha2 = "0.05"
gamma  = "0.05"
eta    = "0.2"
delta  = "0.6"
sigma0 = "0.5"
tau = "10"
lambda = "0"

#these are the baseline hyperparameters
alpha0 = "0.2"
alpha1 = "0.05"
alpha2 = "0.05"
gamma  = "0.01"
eta    = "1.2"
delta  = "0.6"
sigma0 = "0.5"
lambda = "1"

eta_ptsd = "5"


pass=c()
p_values=c()
pass_x=c()
p_values_x=c()
n=100
noise=0.2

seed=42

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


##########################
###Blechert et al. 2007###
##########################

phase_def <- list(

  phase1 = list(
    trials = list(
      A = 6,
      B = 6
    ),
    rewards = list(
      c(A = 1)
    )
  ),

  phase2 = list(
    trials = list(
      A = 6,
      B = 6
    )
  )
)

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=6, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=6, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, fixed_axis=T)


df_long <- do.call(rbind, c(
  lapply(ctrl, function(r) subset(r$long[[2]])),
  lapply(ptsd, function(r) subset(r$long[[2]]))
))

df_avg <- df_long %>%
  group_by(subject, group, cue_type, block) %>%
  summarise(V = mean(V), .groups = "drop")

ext = ezANOVA(data = df_avg, dv = V, wid = subject,
              within = c(cue_type, block), between = group)

df_CSplus <- df_avg %>%
  filter(cue_type == "A")

CSplus = ezANOVA(data = df_CSplus, dv = V, wid = subject,
                 within = block, between = group)

df_CSmin <- df_avg %>%
  filter(cue_type == "B")

CSmin = ezANOVA(data = df_CSmin, dv = V, wid = subject,
                within = block, between = group)

pass=c(pass, "Blechert 2007")
p_values=c(p_values, "Blechert 2007")

#significant
check_p(ext, "group", TRUE)
check_p(ext, "group:cue_type", TRUE)
check_p(ext, "cue_type", TRUE)
check_p(CSplus, "group", TRUE)

#not significant
check_p(CSmin, "group", FALSE)


#now we're gonna do expectancy ratings. For For extinction, one rating was taken
#on the last extinction trial.
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

pass_x=c(pass_x, "Blechert 2007")
p_values_x=c(p_values_x, "Blechert 2007")

#significant
check_p_x(ext, "group", TRUE)
check_p_x(ext, "cue_type", TRUE)
check_p_other(CSplus$p.value, TRUE)

#not significant
check_p_x(ext, "group:cue_type", FALSE)
check_p_other(CSmin$p.value, FALSE)

############################
###Felmingham et al. 2018###
############################

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


ctrl=run(seed=seed, phase_def, V_noise_sd=noise, n=n, pseudo=2, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def, V_noise_sd=noise, n=n, pseudo=2, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, group="ptsd", prev=n)

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, fixed_axis = T, plot_cue=T)

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

pass=c(pass, "Felmingham 2018")
p_values=c(p_values, "Felmingham 2018")

#significant
check_p(early, "cue_type", TRUE)
check_p(early, "presentation", TRUE)
check_p(early, "group:presentation", TRUE)
check_p(late, "presentation", TRUE)

#not significant
check_p(late, "cue_type", FALSE)
check_p(late, "group:presentation", FALSE)

###########################
###Pohlchen et al., 2020###
###########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=c(12,10,12,1), ezITI=list(c(56,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=c(12,10,12,1), ezITI=list(c(56,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, fixed_axis = T, plot_cue=T)

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

pass=c(pass, "Pohlechen 2020")
p_values=c(p_values, "Pohlechen 2020")

#not significant
check_p(acq, "group", FALSE)
check_p(ext, "group", FALSE)
check_p(rec, "group", FALSE)
check_p(acq, "group:cue_type", FALSE)
check_p(ext, "group:cue_type", FALSE)
check_p(rec, "group:cue_type", FALSE)


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
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "A", block = factor(1:3), V = r$V[["A"]][c(13, 25, 37), 1]),
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "C", block = factor(1:3), V = r$V[["C"]][c(13, 25, 37), 1])
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
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "A", block = factor(1:2), V = r$V[["A"]][c(69, 81), 1]),
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "C", block = factor(1:2), V = r$V[["C"]][c(69, 81), 1])
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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=c(12,10,1), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=c(12,10,1), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

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
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "A", block = factor(1:2), V = r$V[["A"]][c(47, 57), 1]),
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "C", block = factor(1:2), V = r$V[["C"]][c(47, 57), 1])
    )
  })
))


ext_x <- ezANOVA(data = df_ext, dv = V, wid = subject,
               within = c(cue_type, block), between = group)

pass_x=c(pass_x, "Pohlechen 2020")
p_values_x=c(p_values_x, "Pohlechen 2020")

#not significant
check_p_x(acq_x, "group", FALSE)
check_p_x(ext_x, "group", FALSE)
check_p_x(rec_x, "group", FALSE)
check_p_x(acq_x, "group:cue_type", FALSE)
check_p_x(ext_x, "group:cue_type", FALSE)
check_p_x(rec_x, "group:cue_type", FALSE)

##########################
###Norrholm et al. 2011###
##########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=8, eta=eta_ptsd, alpha1=alpha1, lambda=lambda, alpha0=alpha0, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, fixed_axis = T, plot=T)

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

pass=c(pass, "Norrholm 2011")
p_values=c(p_values, "Norrholm 2011")

#significant
check_p(ext, "group", TRUE)
check_p_other(ext_sb1$p.value, TRUE)
check_p_other(ext_sb2$p.value, TRUE)

#not significant
check_p(acq, "group", FALSE)
check_p(ext_CSminus, "group", FALSE)
check_p_other(ext_sb3$p.value, FALSE)


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

pass_x=c(pass_x, "Norrholm 2011")
p_values_x=c(p_values_x, "Norrholm 2011")

#not significant
check_p_x(acq_x, "group", FALSE)
check_p_x(ext_x, "group", FALSE)

#######################
###Shvil et al. 2014###
#######################

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

ctrl=run(seed=seed, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, plot_cue=T, fixed_axis = T)

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


pass=c(pass, "Shvil 2014")
p_values=c(p_values, "Shvil 2014")

#not significant
check_p(ext, "cue_type", FALSE)
check_p(ext, "group", FALSE)
check_p(ext, "group:cue_type", FALSE)

check_p(rec, "cue_type", FALSE)
check_p(rec, "group", FALSE)
check_p(rec, "group:cue_type", FALSE) #original study had p value of 0.06, trend for significance

#there was also a trend for better extinction recall index in controls (0.07),
#which we did not recreate in code.

##########################
###Acheson et al. 2015####
##########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=8, eta=eta, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, lambda=lambda)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=8, eta=eta_ptsd, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, lambda=lambda, prev=n, group="ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, fixed_axis = T)

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

pass=c(pass, "Acheson 2015")
p_values=c(p_values, "Acheson 2015")

#significant
check_p(acq, "group:cue_type", TRUE)
check_p_other(ctrl_ttest$p.value, TRUE)
check_p(ext, "group", TRUE)
check_p_other(ext_block_3$p.value, TRUE)
check_p_other(ext_block_4$p.value, TRUE)

#not significant
check_p_other(ptsd_ttest$p.value, FALSE) #ptsd only had p value of 0.09, driven by higher CS- response in acq. There was no group differences in expectancy ratings
check_p_other(ext_block_1$p.value, FALSE)
check_p_other(ext_block_2$p.value, FALSE)


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


pass_x=c(pass_x, "Acheson 2015")
p_values_x=c(p_values_x, "Acheson 2015")

#not significant
check_p_x(acq_x, "group", FALSE)
check_p_x(acq_x, "group:cue_type", FALSE)
check_p_x(ext_x, "group", FALSE)
check_p_x(ext_x, "group:presentation", FALSE)

#########################
###Pineles et al. 2016###
#########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, pseudo=3, ezITI=list(c(30,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, pseudo=3, ezITI=list(c(30,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, fixed_axis = T, plot=T)

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
  group = factor(c(rep("ctrl", n),
                   rep("ptsd", n)))
)

fit_initial <- lm(initial_change ~ group, data = df)
fit_late <- lm(late_change ~ group, data = df)
fit_ext_retention <- lm(ext_retention ~ group, data = df)

initial_p = summary(fit_initial)$coefficients["groupptsd", "Pr(>|t|)"]
late_p = summary(fit_late)$coefficients["groupptsd", "Pr(>|t|)"]
retention_p = summary(fit_ext_retention)$coefficients["groupptsd", "Pr(>|t|)"]
retention_est <- summary(fit_ext_retention)$coefficients["groupptsd", "Estimate"]

pass=c(pass, "Pineles 2016")
p_values=c(p_values, "Pineles 2016")

#not significant
check_p_other(initial_p, FALSE)
check_p_other(late_p, FALSE) #although the study got a p value of 0.06, indicating marginal significance

#significant (specifically, extinction retention is worse for PTSD, as they show more day 2 fear)
pass <<- c(pass, retention_p < 0.05 && retention_est > 0)
if (retention_est > 0){
  p_values=c(p_values, retention_p)
} else {
  p_values=c(p_values, "Wrong direction")
}

##########################
###Norrholm et al. 2015###
##########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block = 8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block = 8, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, fixed_axis = T)


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
  group = factor(c(rep("ctrl", n),
                   rep("ptsd", n)))
)

early_p = t.test(early_ctrl, early_ptsd)
mid_p = t.test(mid_ctrl, mid_ptsd)
late_p = t.test(late_ctrl, late_ptsd)

pass=c(pass, "Norrholm 2015")
p_values=c(p_values, "Norrholm 2015")

#significant
check_p_other(early_p$p.value, TRUE)
check_p_other(mid_p$p.value, TRUE)

#not significant
check_p_other(late_p$p.value, FALSE)

#######################
###Handy et al. 2018###
#######################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, pseudo=3, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, pseudo=3, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl,ptsd=ptsd, plot=T)

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

pass=c(pass, "Handy 2018")
p_values=c(p_values, "Handy 2018")

#significant
check_p_other(ext_p, TRUE)

#########################
###Burriss et al. 2006###
#########################

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

n_participants <- n

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
ctrl <- vector("list", n_participants)
for (i in seq_len(n_participants)) {
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
ptsd <- vector("list", n_participants)
for (i in seq_len(n_participants)) {
  set.seed(seed+i*100+1e6)
  my_world_i <- build_phase1_world()
  ptsd[[i]] <- run(
    seed=seed+i*100+1e6, phase_def = phase_def, V_noise_sd=noise, n = 1, my_world = my_world_i,
    eta = eta_ptsd, lambda = lambda, alpha0 = alpha0, alpha1 = alpha1,
    alpha2 = alpha2, gamma = gamma, delta = delta, sigma0 = sigma0,
    prev = n_participants + i - 1, group = "ptsd"
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

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T)

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

pass=c(pass, "Burriss 2006")
p_values=c(p_values, "Burriss 2006")

#significant
check_p(ext, "block", TRUE)

#not significant
check_p(acq, "group:block", FALSE)
check_p(ext, "group:block", FALSE)
check_p(ext, "cue_type:block", FALSE)

#############################
###Grillon and Morgan 1999###
#############################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=c(10,4), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=c(10,4), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl,ptsd=ptsd, plot=T)

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

pass=c(pass, "Grillon and Morgan 1999")
p_values=c(p_values, "Grillon and Morgan 1999")

#significant
check_p_other(ctrl_acq, TRUE)
check_p_other(ctrl_ext, TRUE)
check_p_other(ext$anova_table["group:cue_type", "Pr(>F)"], TRUE)

#not significant
check_p_other(acq$anova_table["group:cue_type", "Pr(>F)"], FALSE)
check_p_other(ptsd_acq, FALSE)
check_p_other(ptsd_ext, FALSE)

#######################
###Milad et al. 2009###
#######################

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

ctrl=run(seed=seed, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n, lambda=lambda, eta=eta, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3) #I guessed pseudo since they only said "pseudorandomization" without a specific number
ptsd=run(seed=seed + 1e6, phase_def=phase_def, features=c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), V_noise_sd=noise, n=n, lambda=lambda, eta=eta_ptsd, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, pseudo=3, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T, fixed_axis = T)

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


pass=c(pass, "Milad 2009")
p_values=c(p_values, "Milad 2009")

#significant
check_p(rec, "group:cue_type", TRUE)

#not significant
check_p(ext, "cue_type", FALSE)
check_p(ext, "group", FALSE)
check_p(ext, "group:cue_type", FALSE)

#####################
###Zuj et al. 2017###
#####################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=10, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=10, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl,ptsd=ptsd, plot=T)

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


pass=c(pass, "Zuj 2017")
p_values=c(p_values, "Zuj 2017")

#significant
check_p(early, "cue_type", TRUE)
check_p(early, "presentation", TRUE)
check_p(early, "group:cue_type", TRUE)
check_p(late, "presentation", TRUE)
check_p_other(ctrl_1v2$p.value, TRUE)
check_p_other(ptsd_2v3$p.value, TRUE)

#not significant
check_p(acq, "group", FALSE)
check_p(late, "group", FALSE)
check_p(late, "cue_type", FALSE)
check_p(late, "group:cue_type", FALSE)
check_p_other(ptsd_1v2$p.value, FALSE)
check_p_other(ctrl_2v3$p.value, FALSE)


pass_x=c(pass_x, "Zuj 2017")
p_values_x=c(p_values_x, "Zuj 2017")

#significant
check_p_x(early, "cue_type:presentation", TRUE)
check_p_x(late, "group", TRUE)
check_p_x(late, "cue_type", TRUE)
check_p_x(late, "presentation", TRUE)

#not significant
check_p_x(early, "group:cue_type", FALSE)
check_p_x(late, "group:cue_type", FALSE)

###########################
###Garfinkel et al. 2014###
###########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, features = c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), pseudo=3, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, features = c("A", "B", "C", "X", "Y"), ezITI=list(c(64,100000)), pseudo=3, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl,ptsd=ptsd, plot=T)

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

ptsd_ay_ext <- sapply(ctrl, function(r) {
  mean(subset(r$long[[2]], cue_type == "AY" & as.integer(presentation) > 8)$V)
})

ptsd_cy_ext <- sapply(ctrl, function(r) {
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

ptsd_ay_rec <- sapply(ctrl, function(r) {
  mean(subset(r$long[[3]], cue_type == "AY" & as.integer(presentation) <= 4)$V)
})

ptsd_cy_rec <- sapply(ctrl, function(r) {
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

ptsd_ax_ren <- sapply(ctrl, function(r) {
  mean(subset(r$long[[4]], cue_type == "AX" & as.integer(presentation) <= 4)$V)
})

ptsd_cx_ren <- sapply(ctrl, function(r) {
  mean(subset(r$long[[4]], cue_type == "CX" & as.integer(presentation) <= 8)$V)
})

t_ptsd_ren <- t.test(
  x = ptsd_ax_ren,
  y = ptsd_cx_ren,
  paired = TRUE
)


pass=c(pass, "Garfinkel 2014")
p_values=c(p_values, "Garfinkel 2014")

#significant
check_p(rec, "group:cue_type", TRUE)
check_p(ren, "group:cue_type", TRUE)
check_p_other(t_ptsd_rec$p.value, TRUE)
check_p_other(t_ctrl_ren$p.value, TRUE)
check_p_other(t_renewal$p.value, TRUE)

#not significant
check_p(ext, "group", FALSE)
check_p_other(t_ctrl_ext$p.value, FALSE)
check_p_other(t_ptsd_ext$p.value, FALSE)

check_p_other(t_ctrl_rec$p.value, FALSE)
check_p_other(t_ptsd_ren$p.value, FALSE)

###########################
###Jovanovic et al. 2013###
###########################

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

ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, block=8, eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, block=8, eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T)

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


pass=c(pass, "Jovanovic 2013")
p_values=c(p_values, "Jovanovic 2013")

#significant
check_p(ctrl_acq, "cue_type", TRUE)
check_p(ext_ctrl, "block", TRUE)

#not significant
check_p(ptsd_acq, "cue_type", FALSE)
check_p(ext_ptsd, "block", FALSE)

#########################
###Wicking et al. 2016###
#########################

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
ctrl=run(seed=seed, phase_def=phase_def, V_noise_sd=noise, n=n, pseudo=3, features=c("A","B","X","Y","Z"), ezITI=list(c(120,100000)), eta=eta, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0)
ptsd=run(seed=seed + 1e6, phase_def=phase_def, V_noise_sd=noise, n=n, pseudo=3, features=c("A","B","X","Y","Z"), ezITI=list(c(120,100000)), eta=eta_ptsd, lambda=lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma = gamma, delta=delta, sigma0=sigma0, prev = n, group = "ptsd")

plot_compare(ctrl=ctrl, ptsd=ptsd, plot=T)

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

pass=c(pass, "Wicking 2016")
p_values=c(p_values, "Wicking 2016")

#significant
check_p(ren, "group:cue_type", TRUE)
check_p_other(ren_ttest_ptsd$p.value, TRUE)

#not significant
check_p(acq, "group:cue_type", FALSE)
check_p(ext, "group:cue_type", FALSE)
check_p_other(ren_ttest_ctrl$p.value, FALSE)


#Now we're going to do the threat expectancy analysis. Since V records
#associative value on a given trial before the outcome is shown, we want to
#index the trial directly after each block. This is why we created the bonus
#phase. We do a separate ANOVA for each phase, using only the "end of block"
#trials. We also want to combine individual cues to create the V for a given
#cue + context
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
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "AX",
                 V = r$V[["A"]][61, 1] + r$V[["X"]][61, 1]),
      data.frame(subject = factor(i + n), group = "ptsd", cue_type = "BX",
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


pass_x=c(pass_x, "Wicking 2016")
p_values_x=c(p_values_x, "Wicking 2016")

#significant
check_p_other_x(ren_ttest_ptsd_x$p.value, TRUE)

#not significant
check_p_x(acq_x, "group:cue_type", FALSE)
check_p_x(ext_x, "group:cue_type", FALSE)
check_p_x(ren_x, "group:cue_type", FALSE)
check_p_other_x(ren_ttest_ctrl_x$p.value, FALSE)


