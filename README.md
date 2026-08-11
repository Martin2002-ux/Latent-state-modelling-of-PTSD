# Latent-state-modelling-of-PTSD
“recreating_cochran_cisler_simulations” contains the code to recreate all analyses performed in Cochran and Cisler (2019) using the updated latent state model and the new baseline hyperparameters. It shows that the updated model and hyperparameters can recreate all the phenomena tested in Cochran and Cisler (2019) except for the “memory modification” effect.

“main_analysis” contains the code we used to simulate all experiments and generate the results we put into the various tables. This code runs each simulated experiment 1000 times.

“simplified_analysis” contains the same core code as “main_analysis” but is configured to run a single experiment rather than 1,000 simulations. This made it useful for exploratory hyperparameter tuning, as results could be quickly inspected using plotting functions. All reported results are based on the main analysis script. The simplified version is retained for posterity and for convenient exploration of how hyperparameter choices affect the results.
