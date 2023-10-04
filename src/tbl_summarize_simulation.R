library(kableExtra)

x <- list(
  n900s1 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N900SIGMA21.rds"),
  n400s1 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N400SIGMA21.rds"),
  n225s1 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N225SIGMA21.rds"),
  
  
  n900s01 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N900SIGMA20.1.rds"),
  n400s01 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N400SIGMA20.1.rds"),
  n225s01 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N225SIGMA20.1.rds")
  
)


summ_results <- lapply(x,FUN = function(y){ round(apply(y$full_rslts,2,mean),4)})
summ_results_sd <- lapply(x,FUN = function(y){ round(apply(y$full_rslts,2,sd),4)})

a<- sapply(1:6, FUN = function(i){paste0( summ_results[[i]], " (",summ_results_sd[[i]], ")")})


res_tbl <- t(data.frame(a))
colnames(res_tbl) <- c("beta0","beta1","sigma2eta","phieta")
rownames(res_tbl) <- c("N=900","N=400","N=225","N=900","N=400","N=225")


kable(res_tbl, "latex",booktabs=TRUE,caption="Simulation Results") %>% pack_rows("sigma2eta=1",1,3) %>% pack_rows("sigma1eta=0.1",4,6)



# add making table when you sample
b0sigma1 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N900SIGMA21SAMPSIZE450.rds")
b0sigma0.1 <- readRDS("~/Documents/phd/pf-ext-reg/outputs/full_set_sim/sim-no-sampling-PHIX0.5N900SIGMA20.1SAMPSIZE450.rds")


summb0_results <- lapply(list(b0sigma1,b0sigma0.1),FUN = function(y){ round(apply(y$full_rslts,2,mean),4)})
summb0_results_sd <- lapply(list(b0sigma1,b0sigma0.1),FUN = function(y){ round(apply(y$full_rslts,2,sd),4)})

a2<- sapply(1:2, FUN = function(i){paste0( summb0_results[[i]], " (",summb0_results_sd[[i]], ")")})


res_tbl2 <- t(data.frame(a2))
colnames(res_tbl2) <- c("beta0","beta1","sigma2eta","phieta")
rownames(res_tbl2) <- c("N=900, SIZE=450","N=900, SIZE=450")


kable(res_tbl2, "latex",booktabs=TRUE,caption="Simulation Results, one population of 900 stations and sample size 50% (450 stations)") %>% pack_rows("sigma2eta=1",1,1) %>% pack_rows("sigma1eta=0.1",2,2)