library(dplyr)
library(latentState)
alpha0 = "0.2"
alpha1 = "0.05"
alpha2 = "0.05"
gamma  = "0.01"
eta    = "1.2"
delta  = "0.6"
sigma0 = "0.5"
lambda = "1"
success = c()

#all run functions use the default n value of 1

##############
###Blocking###
##############

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      AB=8,
      CD=8
    ),
    rewards = list(
      c(AB = 1, CD = 1)
    )
  )
)

my_world = make_world(
  trial_patterns = c(rep("A", 20), rep(c("AB", "CD"), 8)),
  reward = rep(1, 36)
)


ctrl=run(phase_def=phase_def, my_world=my_world,plot_raw=T, plot=T)
upd=run(plot_raw=T,phase_def=phase_def, my_world=my_world,alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=T, plot_sd=T, fixed_axis=T,plot_label = T)

#check if B and C have different V at the end of the experiment
B = round(last(upd[[1]][["V_cue"]][["B"]]), 2)
C = round(last(upd[[1]][["V_cue"]][["C"]]), 2)
success = c(success, B<C)

#####################
###Overexpectation###
#####################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=50,
      B=50
    ),
    rewards = list(
      c(A = 1,
        B = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      AB=10
    ),
    rewards = list(
      c(AB = 1)
    )
  )
)


my_world = make_world(trial_patterns = c(rep(c("A", "B"), 50), rep("AB", 10)),
                      reward = rep(1, 110))


ctrl=run(phase_def=phase_def, my_world=my_world,plot_raw=F, plot=T)
upd=run(plot=T,phase_def=phase_def, my_world=my_world,alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=T, plot_sd=T, fixed_axis=T)

#check if A at the end of phase 1 is different than A at the end of phase 2
A1 = round(upd[[1]][["V_cue"]][["A"]][50], 2)
A2 = round(upd[[1]][["V_cue"]][["A"]][60], 2)
success = c(success, A1>A2)

############################
###Conditioned Inhibition###
############################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=100,
      AB=100
    ),
    rewards = list(
      c(A = 1, #this is 100% reward rate of 1 reward
        AB = 1) #this is actually 100% reward rate of 0.5 reward
    )
  )
)

my_world = make_world(trial_patterns = c(rep(c("A", "AB"), 100)),
                      reward = rep(c(1, 0.5), 100))

ctrl=run(phase_def=phase_def, my_world=my_world,plot_raw=F, plot=T)
upd=run(phase_def=phase_def, my_world=my_world, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=T, plot_sd=T, fixed_axis=T)

#check if B has negative value
B = round(last(upd[[1]][["V_cue"]][["B"]]), 2)
success = c(success, B<0)

########################
###Backwards blocking###
########################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      AB=20
    ),
    rewards = list(
      c(AB = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      A=16
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

ctrl=run(phase_def=phase_def, my_world=NULL,plot_raw=F, plot=T)
upd=run(phase_def=phase_def, my_world=NULL, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=T, plot_sd=T, fixed_axis=T)

#check if B at the end of phase 1 is different than B at the end of phase 2
B1 = round(upd[[1]][["V"]][["B"]][20], 2)
B2 = round(upd[[1]][["V"]][["B"]][36], 2)
success = c(success, B1>B2)

#################################
###Rescorla 2000 Experiment 1A###
#################################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=100,
      C=100,
      X=100,
      BX=100,
      DX=100
    ),
    rewards = list(
      c(A = 1,
        C = 1,
        X = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      AB=8
    ),
    rewards = list(
      c(AB=1)
    )
  )
)

ctrl=run(phase_def=phase_def, my_world=NULL,plot_raw=F, plot=F)
upd=run(phase_def=phase_def, my_world=NULL, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=T, plot_sd=T, fixed_axis=F)

#change in associative strength from the end of Stage 1 to the end of Stage 2
#for A vs B
Adif = round(upd[[1]][["V_cue"]][["A"]][100] - upd[[1]][["V_cue"]][["A"]][108], 2)
Bdif = round(abs(upd[[1]][["V_cue"]][["B"]][100] - upd[[1]][["V_cue"]][["B"]][108]), 2)
success = c(success, Adif < Bdif)


#################################
###Rescorla 2000 Experiment 1B###
#################################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=100,
      C=100,
      X=100,
      BX=100,
      DX=100
    ),
    rewards = list(
      c(A = 1,
        C = 1,
        X = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      AB=8
    )
  )
)

ctrl=run(phase_def=phase_def, my_world=NULL,plot_raw=T, plot=T)
upd=run(phase_def=phase_def, my_world=NULL, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)
plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=T, plot_sd=T, fixed_axis=F)

#change in associative strength from the end of Stage 1 to the end of Stage 2
#for A vs B
Adif = round(upd[[1]][["V_cue"]][["A"]][100] - upd[[1]][["V_cue"]][["A"]][108], 2)
Bdif = round(abs(upd[[1]][["V_cue"]][["B"]][100] - upd[[1]][["V_cue"]][["B"]][108]), 2)
success = c(success, Adif > Bdif)

########################
###Wilson et al. 1992###
########################
#group E
phase_def <- list(
  
  phase1 = list(
    trials = list(
      AB=10
    ),
    rewards = list(
      c(AB = 1) #actually alternates between rewarded and unrewarded
    )
  ),
  
  phase2 = list(
    trials = list(
      AB=20,
      A=20
    ),
    rewards = list(
      c(AB = 1) 
    )
  ),
  
  phase3 = list(
    trials = list(
      A=10
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

my_world = make_world(trial_patterns = c(rep(c("AB"), 10), rep(c("AB", "A"), 20), 
                                         rep("A", 10)),
                      reward = c(rep(c(1, 0), 25), rep(1, 10)))

ctrl_E=run(phase_def=phase_def, my_world=my_world,plot_raw=T, plot=T)
upd_E=run(phase_def=phase_def, my_world=my_world,alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)


#group C
phase_def <- list(
  
  phase1 = list(
    trials = list(
      AB=50
    ),
    rewards = list(
      c(AB = 1) #actually alternates between rewarded and unrewarded
    )
  ),
  
  phase2 = list(
    trials = list(
      A=10
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

my_world = make_world(trial_patterns = c(rep(c("AB"), 50), rep("A", 10)),
                      reward = c(rep(c(1, 0), 25), rep(1, 10)))

ctrl_C=run(phase_def=phase_def, my_world=my_world,plot_raw=T, plot=F)
upd_C=run(phase_def=phase_def, my_world=my_world, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl_C=ctrl_C, ctrl_E=ctrl_E, plot=T, plot_avg=F, plot_sd=T, fixed_axis=F)
plot_compare(upd_C=upd_C, upd_E=upd_E, plot=T, plot_avg=F, plot_sd=T, fixed_axis=F)

#does cue A of group E end with higher V than cue A of group C
EA = round(last(upd_E[[1]][["V_cue"]][["A"]]), 2)
CA = round(last(upd_C[[1]][["V_cue"]][["A"]]), 2)
success = c(success, EA>CA)

#######################
###PREE Experiment 1###
#######################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1) #alternating group alternates reward/not, continuous group always rewarded
    )
  ),
  
  phase2 = list(
    trials = list(
      A=20
    )
  )
)

#alternating group
my_world = make_world(trial_patterns = rep(c("A"), 40),
                      reward = c(rep(c(1, 0), 10), rep(0, 20)))

ctrl_alt=run(phase_def=phase_def, my_world=my_world,plot_raw=F, plot=T)
upd_alt=run(phase_def=phase_def, my_world=my_world, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

#continuous group
my_world = make_world(trial_patterns = rep(c("A"), 40),
                      reward = c(rep(1, 20), rep(0, 20)))

ctrl_con=run(phase_def=phase_def, my_world=my_world)
upd_con=run(phase_def=phase_def, my_world=my_world, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl_alt=ctrl_alt, ctrl_con=ctrl_con, plot=T)
plot_compare(upd_alt=upd_alt, upd_con=upd_con, plot=T)

#is cue A higher at the start of extinction for the alternating group than
#for the continuous group. Since Cochran Cisler did not specify exactly which
#trials should be checked, we will look at the average of trials 1 to 10
alt = round(mean(upd_alt[[1]][["V_cue"]][["A"]][21:30]), 2)
con = round(mean(upd_con[[1]][["V_cue"]][["A"]][21:30]), 2)
success = c(success, alt>con)

#######################
###PREE Experiment 2###
#######################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1) #alternating group alternates reward/not, continuous group always rewarded
    )
  ),
  
  phase2 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1) 
    )
  ),
  
  phase3 = list(
    trials = list(
      A=20
    )
  )
)

#alternating group
my_world = make_world(trial_patterns = rep(c("A"), 60),
                      reward = c(rep(c(1, 0), 10), rep(1, 20), rep(0, 20)))


ctrl_alt=run(phase_def=phase_def, my_world=my_world,plot_raw=F, plot=F)
upd_alt=run(phase_def=phase_def, my_world=my_world,alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

#continuous group
my_world = make_world(trial_patterns = rep(c("A"), 60),
                      reward = c(rep(1, 40), rep(0, 20)))

ctrl_con=run(phase_def=phase_def, my_world=my_world)
upd_con=run(phase_def=phase_def, my_world=my_world, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl_alt=ctrl_alt, ctrl_con=ctrl_con, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)
plot_compare(upd_alt=upd_alt, upd_con=upd_con, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)

#is cue A higher at the start of extinction for the alternating group than
#for the continuous group
alt = round(mean(upd_alt[[1]][["V_cue"]][["A"]][41:50]), 2)
con = round(mean(upd_con[[1]][["V_cue"]][["A"]][41:50]), 2)
success = c(success, alt>con)

############################
###Renewal (rapid return)###
############################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      A=20
    )
  ),
  
  phase3 = list(
    trials = list(
      A=10
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

ctrl=run(phase_def=phase_def, plot_raw=T)
upd=run(phase_def=phase_def, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl=ctrl, upd=upd, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)

#the Cochran Cisler paper used the difference between acquisition trial 2 and 
#renewal trial 2 to measure strength of renewal. However, I am only interested
#in whether the effect can be reproduced. Therefore, my test quantity will be
#whether renewal has higher average value than acquisition, which could not be
#achieved if renewal phase began at V=-0.5
#I exclude trial 41 from renewal as V is calculated before the agent sees reward
#and so they behave as is it was still extinction
acq = round(mean(upd[[1]][["V_cue"]][["A"]][1:20]), 2)
ren = round(mean(upd[[1]][["V_cue"]][["A"]][41:50]), 2)
success = c(success, ren>acq)

############################
###Renewal (w/ context)#####
############################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      A=20
    )
  ),
  
  phase3 = list(
    trials = list(
      A=10
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

ctrl_ctx=run(phase_def=phase_def, contextShift = c(21, 41))
upd_ctx=run(phase_def=phase_def, contextShift = c(21, 41), alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

ctrl_non=run(phase_def=phase_def)
upd_non=run(phase_def=phase_def, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl_ctx=ctrl_ctx, ctrl_non=ctrl_non, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)
plot_compare(upd_ctx=upd_ctx, upd_non=upd_non, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)

#test whether renewal is superior when there is context change
ctx = round(upd_ctx[[1]][["V_cue"]][["A"]][41], 2)
non = round(upd_non[[1]][["V_cue"]][["A"]][41], 2)

success = c(success, ctx > non)

##########################################################
###Renewal (w/ context) using cues to represent context###
##########################################################

#This is a repeat of the above simulation, but using cues rather than the
#contextShift mechanism. This justifies my choice to use cues to represent
#context in my main simulations.
phase_def <- list(
  
  phase1 = list(
    trials = list(
      AB=20
    ),
    rewards = list(
      c(AB = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      AC=20
    )
  ),
  
  phase3 = list(
    trials = list(
      AB=10
    ),
    rewards = list(
      c(AB = 1)
    )
  )
)

ctrl_ctx=run(phase_def=phase_def, features=c("A", "B", "C"))
upd_ctx=run(phase_def=phase_def, features=c("A", "B", "C"), alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      A=20
    )
  ),
  
  phase3 = list(
    trials = list(
      A=10
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

ctrl_non=run(phase_def=phase_def)
upd_non=run(plot=T,phase_def=phase_def, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)


#test whether renewal is superior when there is context change
ctx = round(upd_ctx[[1]][["V_shown"]][["AB"]][21], 2)
non = round(upd_non[[1]][["V_shown"]][["A"]][41], 2)

success = c(success, ctx > non) 

##########################
###Spontaneous recovery###
##########################

phase_def <- list(
  
  phase1 = list(
    trials = list(
      A=20
    ),
    rewards = list(
      c(A = 1)
    )
  ),
  
  phase2 = list(
    trials = list(
      A=20
    )
  ),
  
  phase3 = list(
    trials = list(
      A=10
    ),
    rewards = list(
      c(A = 1)
    )
  )
)

# Baseline timestamps: 1, 2, ..., 50 (no extra temporal gap)
ITI_vec <- 1:50

# Add a 200-unit temporal gap between trial 40 and 41:
# gap = ITI[41] - ITI[40] - 1 = 200
ITI_vec[41:50] <- ITI_vec[41:50] + 199

# Wrap in a list because run() expects a list of length n
ITI_list <- list(ITI_vec)

ctrl_iti=run(phase_def=phase_def, ITI = ITI_list)
upd_iti=run(phase_def=phase_def, ITI = ITI_list, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

ctrl_non=run(phase_def=phase_def)
upd_non=run(phase_def=phase_def, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl_iti=ctrl_iti, ctrl_non=ctrl_non, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)
plot_compare(upd_iti=upd_iti, upd_non=upd_non, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)

#test whether renewal is superior when there is context change
iti = round(upd_iti[[1]][["V_cue"]][["A"]][41], 2)
non = round(upd_non[[1]][["V_cue"]][["A"]][41], 2)

success = c(success, iti > non) #the difference is only 0.02, but even with ctrl parameters it's only 0.04


#########################
###Memory modification###
#########################

phase_def <- list(
  #acquisition
  phase1 = list(
    trials = list(
      A=3
    ),
    rewards = list(
      c(A = 1)
    )
  ),
  #non-reinforced retrieval
  phase2 = list(
    trials = list(
      A=1
    )
  ),
  #extinction
  phase3 = list(
    trials = list(
      A=18
    )
  ),
  #non-reinforced test (2 trials)
  phase4 = list(
    trials = list(
      A=2
    )
  )
)

ntrials <- 24

# Baseline timestamps (consecutive trials = gap of 1)
ITI_vec_5_delay <- 1:ntrials

# Gap of 21 at t=3→4 (nIter fires at t=3: time(4)-time(3)=21)
ITI_vec_5_delay[4:ntrials] <- ITI_vec_5_delay[4:ntrials] + 20

# Gap of 6 at t=4→5 (nIter fires at t=4: time(5)-time(4)=6)
ITI_vec_5_delay[5:ntrials] <- ITI_vec_5_delay[5:ntrials] + 5

# Gap of 201 at t=23→24 (nIter fires at t=23: time(24)-time(23)=201)
ITI_vec_5_delay[24:ntrials] <- ITI_vec_5_delay[24:ntrials] + 200

# Wrap for run()
ITI_list_5_delay <- list(ITI_vec_5_delay)


ctrl_5 <- run(phase_def=phase_def, ITI=ITI_list_5_delay, chi="5")
upd_5 <- run(phase_def=phase_def, ITI=ITI_list_5_delay, chi="5", alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

# Baseline timestamps
ITI_vec_1_delay <- 1:ntrials

# Gap of 21 at t=3→4
ITI_vec_1_delay[4:ntrials] <- ITI_vec_1_delay[4:ntrials] + 20

# Gap at t=4→5 stays 1
ITI_vec_1_delay[5:ntrials] <- ITI_vec_1_delay[5:ntrials] + 1

# Gap of 201 at t=23→24
ITI_vec_1_delay[24:ntrials] <- ITI_vec_1_delay[24:ntrials] + 200

# Wrap for run()
ITI_list_1_delay <- list(ITI_vec_1_delay)


ctrl_1 <- run(phase_def=phase_def, ITI=ITI_list_1_delay, chi="5")
upd_1 <- run(phase_def=phase_def, ITI=ITI_list_1_delay, chi="5", alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

# Baseline timestamps
ITI_vec_100_delay <- 1:ntrials

# Gap of 21 at t=3→4
ITI_vec_100_delay[4:ntrials] <- ITI_vec_100_delay[4:ntrials] + 20

# Gap of 101 at t=4→5 
ITI_vec_100_delay[5:ntrials] <- ITI_vec_100_delay[5:ntrials] + 100

# Gap of 201 at t=23→24
ITI_vec_100_delay[24:ntrials] <- ITI_vec_100_delay[24:ntrials] + 200

# Wrap for run()
ITI_list_100_delay <- list(ITI_vec_100_delay)


ctrl_100 <- run(phase_def=phase_def, ITI=ITI_list_100_delay, chi="5")
upd_100 <- run(phase_def=phase_def, ITI=ITI_list_100_delay, chi="5", alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, gamma=gamma, eta=eta, delta=delta, sigma0=sigma0, lambda=lambda)

plot_compare(ctrl_1=ctrl_1, ctrl_5=ctrl_5, ctrl_100=ctrl_100, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)
plot_compare(upd_1=upd_1, upd_5=upd_5, upd_100=upd_100, plot=T, plot_avg=F, plot_sd=T, fixed_axis=T)

del1   <- round(upd_1[[1]][["V_cue"]][["A"]][24],   2)
del5   <- round(upd_5[[1]][["V_cue"]][["A"]][24],   2)
del100 <- round(upd_100[[1]][["V_cue"]][["A"]][24],   2)

success=c(success, del1 > del100 & del100 > del5)

success
