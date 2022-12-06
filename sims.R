options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, scales, data.table, r2r, checkmate, rlist)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

#number of cores to use
SEED = 1984
set.seed(SEED)
nC = 8
ns = c(48, 112, 256)
nX = 1
ps = c(1, 2, 5, 10)


#data parameters
bbeta_vec = c(0.5, 0.5, rep(0, max(ps) - 2))
response_params = list(
  "continuous" = list(
    beta0 = 0,
    betaT = 1,
    betas = bbeta_vec,
    sigma = 1
  ),
  "binary" = list(
    beta0 = 0,
    betaT = 1,
    betas = bbeta_vec   
  ),
  "count" = list(
    beta0 = 0,
    betaT = 1,
    betas = bbeta_vec     
  ),
  "survival" = list(
    beta0 = 0,
    betaT = 1,
    betas = bbeta_vec,
    shape = 1
  ),
  "proportion" = list(
    beta0 = 0,
    betaT = 1,
    betas = bbeta_vec,
    shape = 1
  )  
)
rm(bbeta_vec)
response_types = names(response_params)

expit = function(u){1 / (1 + exp(-u))}

survival_scale = 0.5
response_mean_functions = list(
  "continuous" = function(n, X, bbeta, w, betaT){
    X %*% bbeta + w * betaT
  },
  "binary" = function(n, X, bbeta, w, betaT){
    expit(X %*% bbeta + w * betaT)
  },
  "count" = function(n, X, bbeta, w, betaT){
    exp(X %*% bbeta + w * betaT)
  },
  "survival" = function(n, X, bbeta, w, betaT){
    exp(X %*% bbeta + w * betaT) * gamma(1 + 1 / survival_scale)
  },
  "proportion" = function(n, X, bbeta, w, betaT){
    expit(X %*% bbeta + w * betaT)
  }
)

continuous_sd = 2
proportion_phi = 2
response_dgp_functions = list(
  "continuous" = function(n, mu){
    rnorm(n, mean = mu, sd = continuous_sd)
  },
  "binary" = function(n, mu){
    rbinom(n, size = 1, prob = mu)
  },
  "count" = function(n, mu){
    rpois(n, lambda = mu)
  },
  "survival" = function(n, mu){
    rweibull(n, shape = mu / gamma(1 + 1 / survival_scale)\beta_2, scale = survival_scale)
  },
  "proportion" = function(n, mu){
    rbeta(n, shape1 = mu * proportion_phi, shape2 = (1 - mu) * proportion_phi)
  }
)

covariate_params = list(
  "continuous" = list(
    dgp = rnorm,
    param_1 = 0,
    param_2 = 1
  ),
  "binary" = list(
    dgp = rnorm,
    param_1 = 0,
    param_2 = 1
  ),
  "count" = list(
    dgp = rnorm,
    param_1 = 0,
    param_2 = 0.2    
  ),
  "survival" = list(
    dgp = rnorm,
    param_1 = 0,
    param_2 = 0.1
  ),
  "proportion" = list(
    dgp = rnorm,
    param_1 = 0,
    param_2 = 1  
  )  
)

Xalls = list()
for (response_type in response_types){
  Xalls[[response_type]] = list()
  for (nx in 1 : nX){
    Xalls[[response_type]][[nx]] = hashmap()
    for (n in ns){
      Xalls[[response_type]][[nx]][[n]] = hashmap()
      dgp = covariate_params[[response_type]]$dgp
      param1 = covariate_params[[response_type]]$param_1
      param2 = covariate_params[[response_type]]$param_2
      for (p in ps){
        Xtemp = cbind(1, matrix(dgp(n * p, param1, param2), nrow = n, ncol = p))
        #WLOG order the matrix by first covariate
        Xtemp = Xtemp[order(Xtemp[,2]), ]
        if (p > 1){
          #WLOG order the quarters of the matrix by the second covariate
          X1 = Xtemp[1 : (n/4), ]
          X2 = Xtemp[(n/4 + 1) : (n/2), ]
          X3 = Xtemp[(n/2 + 1) : (3*n/4), ]
          X4 = Xtemp[(3*n/4 + 1) : n, ]
          X1 = X1[order(X1[,3]), ]
          X2 = X2[order(X2[,3]), ]
          X3 = X3[order(X3[,3]), ]
          X4 = X4[order(X4[,3]), ]
          Xalls[[response_type]][[nx]][[n]][[p]] = rbind(X1, X2, X3, X4)
        } else {
          Xalls[[response_type]][[nx]][[n]][[p]] = Xtemp
        }
      }
    }
  }
}
rm(covariate_params, dgp, param1, param2, Xtemp, X1, X2, X3, X4, n, p, response_type, nx)

gen_block_designs = function(n, B, nR){
  n_over_2B = n / (2 * B)
  assertCount(n_over_2B) #balance within the blocks
  
  dummy_block = c(rep(1, n_over_2B), rep(-1, n_over_2B))
  
  Ws = list()
  for (b in 1 : B){
    Ws[[b]] = matrix(NA, nrow = nR, ncol = n_over_2B * 2)
    for (nr in 1 : nR){
      Ws[[b]][nr, ] = sample(dummy_block)
    }
  }
  list.cbind(Ws)
}


#design parameters
designs = c("BCRD", "PB", "PM", "B2", "B4", "B8", "R", "KK")
objective = "mahal_dist"
rerand_threshold = 0.01
nR = 3e5

#now generate all designs for all designs
if (!exists("all_ws")){
  all_ws = hashmap()
}

for (n in ns){
  if (is.null(all_ws[[n]])){
    all_ws[[n]] = hashmap()
  }
  

  cat("generating bcrd designs...\n")
  w_bcrds = complete_randomization_with_forced_balanced(n, nR, form = "pos_one_min_one")          
  
  #for blocking, we already have the first and second feature ordered
  cat("generating blocking designs...\n")
  w_two_blocks =   gen_block_designs(n, B = 2, nR)
  w_four_blocks =  gen_block_designs(n, B = 4, nR)
  w_eight_blocks = gen_block_designs(n, B = 8, nR)
  
  for (response_type in response_types){
    
    if (is.null(all_ws[[n]][[response_type]])){
      all_ws[[n]][[response_type]] = hashmap()
    }
    
    for (nx in 1 : nX){
      if (is.null(all_ws[[n]][[response_type]][[nx]])){
        all_ws[[n]][[response_type]][[nx]] = hashmap()
      }
      
      for (p in ps){
        cat("n", n, "response", response_type, "nx", nx, "p", p, "\n")
        
        if (is.null(all_ws[[n]][[response_type]][[nx]][[p]])){
          all_ws[[n]][[response_type]][[nx]][[p]] = hashmap()
        }
        
        
        #already done and standard for all response types, nx, p
        all_ws[[n]][[response_type]][[nx]][[p]][["BCRD"]] = w_bcrds
        all_ws[[n]][[response_type]][[nx]][[p]][["B2"]] = w_two_blocks
        all_ws[[n]][[response_type]][[nx]][[p]][["B4"]] = w_four_blocks
        all_ws[[n]][[response_type]][[nx]][[p]][["B8"]] = w_eight_blocks
        
        #for all other designs, we need the subject-specific covariates
        X = Xalls[[response_type]][[nx]][[n]][[p]]
        Xwithoutint = X[, -1, drop = FALSE]

        if (is.null(all_ws[[n]][[response_type]][[nx]][[p]][["PM"]])){
          cat("generating PM designs...\n") #not parallelized...
          PM_match_structure = computeBinaryMatchStructure(Xwithoutint, mahal_match = TRUE)
          PM_des = initBinaryMatchExperimentalDesignSearch(PM_match_structure, num_cores = nC, wait = TRUE, max_designs = min(nR, 2^(n/2)), seed = SEED)
          w_PM_res = resultsBinaryMatchSearch(PM_des, form = "pos_one_min_one")
          all_ws[[n]][[response_type]][[nx]][[p]][["PM"]] = w_PM_res
        }
        
        #rerand
        
        if (is.null(all_ws[[n]][[response_type]][[nx]][[p]][["R"]])){
          cat("generating rerand designs...\n")
          rerand_des = initRerandomizationExperimentalDesignObject(Xwithoutint, objective = objective, obj_val_cutoff_to_include = Inf, num_cores = nC, wait = TRUE, max_designs = nR, seed = SEED)
          rerand_des_res = resultsRerandomizationSearch(rerand_des, include_assignments = TRUE, form = "pos_one_min_one")
          w_rerands = rerand_des_res$ending_indicTs[order(rerand_des_res$obj_vals)[1 : (nR * rerand_threshold)], ]
          all_ws[[n]][[response_type]][[nx]][[p]][["R"]] = w_rerands
        }
        
        if (is.null(all_ws[[n]][[response_type]][[nx]][[p]][["KK"]])){
          cat("generating greedy KK designs + perfect balance design...\n")
          KK_des = initGreedyExperimentalDesignObject(Xwithoutint, objective = objective, num_cores = nC, wait = TRUE, max_designs = nR, seed = SEED)
          KK_des_res = resultsGreedySearch(KK_des, max_vectors = nR, form = "pos_one_min_one")
          
          all_ws[[n]][[response_type]][[nx]][[p]][["KK"]] = KK_des_res$ending_indicTs
          #the perfect balance design is defined as the best greedy design
          all_ws[[n]][[response_type]][[nx]][[p]][["PB"]] = KK_des_res$ending_indicTs[1, , drop = FALSE]         
        }
        
        rm(PM_match_structure, PM_des, w_PM_res, rerand_des, rerand_des_res, w_rerands, KK_des, KK_des_res)
        gc()
      }
    }
  }
}
all_ws
rm(w_bcrds, w_four_blocks, gen_block_designs, objective, rerand_threshold, nR, X, Xwithoutint, n, p, response_type, nx)

#sim params
nY = 2000 #the number of y vectors to simulate
num_w_max = 1000 #number of vectors to average over for each design
filename = "all_response_model_sims"

exp_settings = expand.grid(
  nx = 1 : nX,
  design = designs,
  n = ns,
  p = ps,
  nsim = 1 : nY
)


# cl = makeCluster(nC)
# # registerDoParallel(cl)
# registerDoSEQ(cl)
# clusterExport(cl, list(
#   "exp_settings", "Xalls", "response_dgp_functions", "response_mean_functions", "response_params", "num_w_max"
# ), envir = environment())

# for (n_setting in 1 : nrow(exp_settings)){
res = foreach(n_setting = 1 : nrow(exp_settings), .inorder = FALSE, .combine = rbind) %do% {
  inner_res = data.table(
    nx = numeric(), 
    response_type = character(), 
    n = numeric(), 
    p = numeric(), 
    design = character(), 
    nsim = numeric(), 
    nW = numeric(),
    tau = numeric(), 
    avg_tauhat = numeric(),
    avg_sq_err = numeric()
  )
  
  for (response_type in response_types){
    #get this run's settings
    nx = exp_settings$nx[n_setting]
    n = exp_settings$n[n_setting]
    p = exp_settings$p[n_setting]
    design = as.character(exp_settings$design[n_setting])
    nsim = exp_settings$nsim[n_setting]
    
    if (p > 2 & design %in% c("B2", "B4", "B8")){
      next
    }
    
    #get this run's data
    X = Xalls[[response_type]][[nx]][[n]][[p]]
    
    #get this run's tau parameter by simulating the potential outcomes
    betaT = response_params[[response_type]]$betaT
    bbeta = c(response_params[[response_type]]$beta0, response_params[[response_type]]$betas)
    mean_func = response_mean_functions[[response_type]]
    muT = mean_func(n, X, bbeta[1 : (p + 1)], rep(+1, n), betaT)
    muC = mean_func(n, X, bbeta[1 : (p + 1)], rep(-1, n), betaT)
    response_dgp = response_dgp_functions[[response_type]]
    
    if (response_type == "continuous"){
      yT = response_dgp(n, muT)
      z = yT - muT #let zT = zC for the continuous response dgp  
      yC = muC + z
    } else {
      yT = response_dgp(n, muT)
      yC = response_dgp(n, muC)
    }
    tau = mean(yT - yC)

    
    #get the w's for this run
    W = all_ws[[n]][[response_type]][[nx]][[p]][[design]] #all designs
    # mean(abs(W%*%X[,2]))
    if (is.null(W)){
      cat("W not calculated\n")
      stop()
    }
    
    if (design != "PB"){
      if (nrow(W) < num_w_max){
        cat("not enough vectors in W\n")
        stop()
      }
      W = W[sample(1 : nrow(W), num_w_max), ]
    }
    nW = nrow(W)
    
    tauhats = array(NA, nW)
    for (i_w in 1 : nW){
      w = W[i_w, ]
      tauhats[i_w] = mean(yT[w == 1]) - mean(yC[w == -1])
    }
    avg_tauhat = mean(tauhats)
    avg_sq_err = mean((tauhats - tau)^2)
    cat("run nsim", nsim, "nx", nx, "n", n, "p", p, "response_type", response_type, "design", design, " nW", nW, "tau", tau, "tauhat", mean(tauhats), "\n")
    
    inner_res = rbind(inner_res, data.table(
      nx = nx, 
      response_type = response_type, 
      n = n, 
      p = p, 
      design = design, 
      nsim = nsim, 
      nW = nW,
      tau = tau, 
      avg_tauhat = avg_tauhat,
      avg_sq_err = avg_sq_err
    ))
  }
  

  # save(inner_res, file = paste0(filename, "_", n_setting, ".RData"))
  inner_res
}
# stopCluster(cl)

save(res, file = paste0(filename, ".RData"))




pacman::p_load(tidyverse, magrittr, data.table, scales)
# load(paste0(filename, ".RData"))

res[, n := factor(n)]
res[, ase := mean(avg_sq_err), by = .(n, p, response_type, design)]
res[, q95 := quantile(avg_sq_err, .95), by = .(n, p, response_type, design)]

res_low_n = res[n == 48]
res_low_n

# design_cols = hue_pal()(6)


ggplot(res_low_n) + 
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), linetype = "dotted") + 
  geom_vline(aes(xintercept = q95, color = design)) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n) + 
  geom_density(aes(x = log(avg_sq_err), fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = log(ase), color = design), linetype = "dotted") + 
  geom_vline(aes(xintercept = log(q95), color = design)) +   
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "continuous" & p == 1  & design %in% c("PM", "PB", "BCRD")]) + 
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), linetype = "dotted") + 
  geom_vline(aes(xintercept = q95, color = design)) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "continuous" & p == 1]) + 
  xlim(-20, 0.5) +
  geom_density(aes(x = log(avg_sq_err), fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = log(ase), color = design), linetype = "dotted") + 
  geom_vline(aes(xintercept = log(q95), color = design)) +   
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))


res[, design := factor(design, levels = c("BCRD", "B2", "B4", "B8", "PM", "PB", "R", "KK"))]
res_high_n = res[n == 256]
res_high_n[order(ase)]
res_high_n


ggplot(res_high_n) + 
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dotted") + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n) + 
  xlim(-10, -1) + 
  geom_density(aes(x = log(avg_sq_err), fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = log(ase), color = design), linetype = "dotted") + 
  geom_vline(aes(xintercept = log(q95), color = design)) +   
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))



###### p = 1

ggplot(res_high_n[response_type == "continuous" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 0.23) + 
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "binary" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 0.0062) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "count" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 0.05) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "proportion" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 0.002) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "survival" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 1.1) +
  xlab("Average Squared Error") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

###### p = 5

ggplot(res_high_n[response_type == "continuous" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.24) + 
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "binary" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.0062) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "count" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.05) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "proportion" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.0019) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "survival" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 1.15) +
  xlab("Average Squared Error") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))


###### p = 1

ggplot(res_high_n[response_type == "continuous" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.04, 0.08) + 
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "binary" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.001, 0.002) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "count" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.005, 0.016) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "proportion" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.0002, 0.00075) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_high_n[response_type == "survival" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0, 0.9) +
  xlab("Average Squared Error") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))




res_high_n[, .(avg_sq_err = mean(avg_sq_err), q = mean(q95)), by = .(n, p, design, response_type)][order(n, response_type, design, p)][1:100]








######################### Supplementary Material


res_low_n = res[n == 48]
res_low_n[order(ase)]
res_low_n

###### p = 1

ggplot(res_low_n[response_type == "continuous" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 1.3) + 
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "binary" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, .029) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "count" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 0.28) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "proportion" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 0.011) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "survival" & p == 1 & design %in% c("B8", "PB", "BCRD")]) + 
  xlim(0, 7.1) +
  xlab("Average Squared Error") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design), size = 1.5) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~paste("response:", .), p = ~paste("p:", .), .multi_line = FALSE))

###### p = 5

ggplot(res_low_n[response_type == "continuous" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 1.3) + 
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "binary" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.0362) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "count" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.25) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "proportion" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 0.01) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "survival" & p == 5 & design %in% c("PM", "PB", "BCRD")]) + 
  xlim(0, 7.5) +
  xlab("Average Squared Error") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))


###### p = 1

ggplot(res_low_n[response_type == "continuous" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.14, 0.51) + 
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "binary" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.002, 0.012) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "count" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.005, 0.105) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "proportion" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0.001, 0.00495) +
  xlab("") + 
  ylab("") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))

ggplot(res_low_n[response_type == "survival" & p == 1 & design %in% c("PM", "BCRD", "B2", "B4", "B8")]) + 
  xlim(0, 7.5) +
  xlab("Average Squared Error") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_density(aes(x = avg_sq_err, fill = design), alpha = 0.5) +
  geom_vline(aes(xintercept = ase, color = design)) + 
  geom_vline(aes(xintercept = q95, color = design), linetype = "dashed", size = 2) + 
  facet_wrap(response_type ~ p, ncol = 4, scales = "free", labeller = labeller(response_type = ~ ., p = ~paste("p:", .), .multi_line = FALSE))
