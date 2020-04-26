# setwd("final_proj")
library(purrr)
source("scripts/functions.R")

# run jags models

n <-  100 # number of simulations to run

# fitting hierarhical models

K <- 5 # iterations

# fitting correct model

# looping so save intermittently in case of failure (repeatedly failed)
for (k in 1:K) {
  jags_dats <- map(1:n, function(i) make_jags_data())
  re_fits <- map2(jags_dats, 1:n, function(x, i) {
    print(i)
    me_jags(x)
  })

  saveRDS( 500, paste0("data/jags_dats_re", k, ".rds"))
  saveRDS(re_fits, paste0("data/jags_re_summaries", k, ".rds"))
}



# fitting incorrect models

for (k in 1:K) {
    jags_dats <- map(1:n, function(i) make_jags_data_inc())
    re_fits <- map2(jags_dats, 1:n, function(x, i) {
        print(i)
        inc_jags(x)
    })
    
    saveRDS( 500, paste0("data/jags_dats_inc", k, ".rds"))
    saveRDS(re_fits, paste0("data/jags_inc_summaries", k, ".rds"))
}
