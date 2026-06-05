# Latent-state-modelling-of-PTSD
Code works as is. Documentation coming soon.

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007331#sec002
This is the original article that proposed this latent state model, by Dr. Amy Cochran and Dr. Josh Cisler (2019). Model architecture has been updated compared to the original. Absent cues are now initiated at their mean value in all other latent states, rather than at 0. An additional hyperparameter lambda has been added. The model now uses entropy to calculate the informativeness of each trial. The likelihood is exponentiated by lambda x normalized entropy, allowing informative trials to more quickly change belief in latent state, proportional to lambda.
