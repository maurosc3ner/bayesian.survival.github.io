# v2 get_legend() added
get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

df_gk_ <- tibble(
  locations = c(
    0.991455371120813, 		
    0.949107912342759, 	 
    0.864864423359769, 	
    0.741531185599394, 	 
    0.586087235467691, 	
    0.405845151377397, 	 
    0.207784955007898, 	
    0.000000000000000
  ),
  weights = c(
    0.022935322010529,
    0.063092092629979,
    0.104790010322250,
    0.140653259715525,
    0.169004726639267,
    0.190350578064785,
    0.204432940075298,
    0.209482141084728
  )
)

df_gk <- df_gk_ %>% 
  dplyr::mutate(locations = -locations) %>% 
  bind_rows(df_gk_) %>% 
  distinct() %>% 
  arrange(locations)

# Evaluate leave-one-out 
# @param llk piecewise log likelihood matrix
# @param ids Vector of multiple ed visits related to patients 
llkByPatient<-function(llk,ids){
  print(dim(llk))
  collapsed<-apply(llk,1L,function(row) tapply(row, ids, sum))
  print(dim(collapsed))
  # return transposed matrix
  return(t(collapsed))
}

# Evaluate the approximate leave-one-out mean posterior survival function using Pareto Smoothed Importance Sampling
#
# @param fit The fitted stan_surv model.
# @param times Vector of times for which the survival functions should be evaluated.
get_loo_mean_surv <- function(fit, times) {
  
  ntimes <- length(times)
  
  # log likelihood for the fitted model
  log_lik <- log_lik(fit)
  
  # extract posterior draws
  post <- as.array(fit) 
  
  # evaluate relative efficiencies (effective samples / total samples)
  r_eff <- relative_eff(
    exp(log_lik), 
    chain_id = rep(1:dim(post)[2], each = dim(post)[1])) 
  
  # evaluate loo
  loo_object <- loo(log_lik, 
                    r_eff     = r_eff,
                    cores     = 2, 
                    save_psis = TRUE)
  # empty list to store predicted survival
  ps <- list() 
  
  # evaluate predicted survival at each time
  for (q in 1:ntimes) {
    ps[[q]] <- posterior_survfit(object        = fit,
                                 times         = times[q], 
                                 extrapolate   = FALSE, 
                                 condition     = FALSE,
                                 return_matrix = TRUE,
                                 draws         = nrow(as.matrix(fit)))
  }
  
  do.call(cbind,
          map(1:ntimes, 
              ~ E_loo(x           = ps[[.]][[1]], 
                      psis_object = loo_object$psis_object, 
                      type        = "mean", 
                      log_ratios  = -log_lik)$value)
  ) 
}



# Evaluate the constrast (two groups) of RMST for the fitted model
#
# @param fit The fitted stan_surv model.
# @param tau The time horizon for evaluating RMST.
rmst_check <- function(fit, tau,newdata=NULL) {
  nd      <- newdata
  
  locs    <- df_gk$locations
  weights <- df_gk$weights
  qpts    <- tau/2 + tau/2 * locs
  # empty list to store predicted survival
  ps <- list() 
  
  # use all MCMC draws so that the sample draws
  # are used at each quadrature point
  ndraws <- nrow(as.matrix(fit))
  # evaluate predicted survival at each quadrature point
  for (q in 1:length(qpts)) {
    ps[[q]] <- posterior_survfit(object        = fit,
                                 newdata       = nd, 
                                 times         = qpts[q], 
                                 extrapolate   = FALSE, 
                                 condition     = FALSE,
                                 return_matrix = TRUE,
                                 draws         = ndraws)
  }
  
  # convert predicted survival to rmst using GK quadrature
  rmst <- 1:length(qpts) %>%
    map(~ ps[[.]][[1]] * weights[[.]] * tau/2) %>%
    Reduce('+', .)
  
  # print(head(rmst))
  
  if(nrow(nd)==1){
    # return only the rmst
    return (tibble(tau=tau,
                   rmstOne = rmst ))
  }else{
    # return the contrast as data frame
    tibble(tau=tau,
           rmstA = rmst[,1], 
           rmstB = rmst[,2],
           diffA_B = rmst[,1]-rmst[,2],
           ratioA_B=rmst[,1]/rmst[,2])
  }
  
}
# Plot RMST for the fitted model for one or two group
#
# @param fit The fitted stan_surv model.
# @param test_df Data frame with non-parametric RMST as calculated by `survRM2::rmst2`.
# @param tau The time horizon for evaluating RMST.
rmst_check_plot <- function(fit,test_df, tau = 12) {
  df_rmst <- rmst_check(fit, tau,newdata = test_df)
  
  if(nrow(test_df)==1){
    p1<-df_rmst %>% 
      bayesplot::mcmc_areas(pars = c("rmstOne"),prob = .95)
    return (list(df_rmst,p1=p1))
  }else{
    p1<-df_rmst %>% 
      bayesplot::mcmc_areas(pars = c("rmstA"),prob = .95) 
    # bayesplot::mcmc_hist(pars = c("rmstAA"), ) 
    # +vline_at(rmst_data$RMST.arm0$rmst["Est."], color='green') 
    
    p2<-df_rmst %>% 
      bayesplot::mcmc_areas(pars = c("rmstB"),prob = .95) 
    # +vline_at(rmst_data$RMST.arm1$rmst["Est."], color='green') 
    
    p3<-df_rmst %>% 
      bayesplot::mcmc_areas(pars = c("diffA_B"),prob = .95) 
    # +vline_at(rmst_data$RMST.arm1$rmst["Est."], color='green') 
    p4<-df_rmst %>% 
      bayesplot::mcmc_areas(pars = c("ratioA_B"),prob = .95) 
    # +vline_at(rmst_data$RMST.arm1$rmst["Est."], color='green') 
    return(list(df_rmst,p1=p1,p2=p2,p3=p3,p4=p4))
  }
}


add_knots <- function(x) {
  knots <- x$basehaz$knots
  if (is.null(knots))
    return(NULL)
  geom_vline(xintercept = knots, color = "green", alpha = 0.5)
}


############ Calibration
calibrate <- function(data, times, y, n_groups = 10,
                      tstart_col,
                      tstop_col,
                      status_col,
                      surv_col){
  n_groups = n_groups
  y2 <- Surv(y%>%pull(tstart_col),
             y%>%pull(tstop_col),
             y%>%pull(status_col))
  #survival object
  ymat <- as.matrix(y2)
  # non-survival object$tsart
  u=y%>%pull(tstart_col)
  
  survhat <- filter(data, !!sym(tstart_col) %in% times) %>% pull(surv_col)
  
  if (isTRUE(length(survhat) != nrow(ymat))){
    stop(paste0("The number of rows in 'y' must equal the number of columns ",
                "in 'x$surv'; that is, the number of observations in 'x' must ",
                "equal the number of observations in 'y'."))
  }
  
  surv <- data.table(u = u,
                     pred = c(survhat))
  surv <- cbind(surv,ymat)
  
  # (1) Cutpoints of predicted S(u|x)
  ntile <- function(x, n){
    cut(x, 
        breaks = quantile(x, 
                          probs = seq(0, 1, length = n + 1), 
                          na.rm = TRUE,
                          type = 2),
        include.lowest = TRUE,
        labels = FALSE)
  }
  if (n_groups > 1){
    #NOTE: ntile() and dplyr::ntile() should be equivalent but ntile() seems to cause more errors
    # when n + 1 quantiles cannot be recreated.
    surv[, interval := dplyr::ntile(pred, n = n_groups), by = "u"] 
  } else{
    surv[, interval := 1]
  }
  
  # (2) Average predicted S(u|x) in each interval
  surv_mean <- surv[, .(pred = mean(pred),
                        n_patients = .N),
                    by = c("u", "interval")]
  
  # (3) Compute KM for patients in each interval
  # y object before survival
  if(ncol(y) == 3){
    kmfit <- survival::survfit(survival::Surv(start, stop, status) ~ survival::strata(interval), 
                               data = surv)
  } else{
    kmfit <- survival::survfit(survival::Surv(time, status) ~ survival::strata(interval), 
                               data = surv)
  }
  kmfit_summary <- summary(kmfit, times = u, extend = TRUE)
  if (n_groups > 1){
    strata <- as.integer(gsub("survival::strata(interval)=interval=", "",
                              kmfit_summary$strata, fixed = TRUE))
  } else{
    strata <- 1
  }
  kmfit_df <-  data.frame(
    interval = strata,
    time = kmfit_summary$time,
    obs = kmfit_summary$surv
  )
  setnames(surv_mean, "u", "time")
  surv_mean <- merge(surv_mean, kmfit_df, 
                     by = c("time", "interval"), 
                     all.x = TRUE)
  # Return
  res <- data.frame(surv_mean)
  res$time<-round(res$time,1)
  class(res) <- c("calibrate", class(res))
  return(res)
}

autoplot.calibrate <- function(object, colour = NULL){
  object$f_time <- factor(object$time,
                          levels = object$time,
                          labels = paste0("Time = ", object$time))
  p <- ggplot(object)
  if (is.null(colour)){
    aes <- aes_string(x = "pred", y = "obs", label = "interval")
  } else{
    aes <- aes_string(x = "pred", y = "obs", col = colour, label = "interval")
  }
  p <- p +
    aes +
    geom_point(size = 1,
               # position=position_dodge(.3),
               position=position_jitter(w = 0.1, h = 0),
               # position = "jitter",
               alpha=0.5) +
    geom_abline(slope = 1) +
    facet_wrap(~f_time) +
    scale_shape_discrete(name = "Model") +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    xlab("Predicted survival probability") + 
    ylab("Observed survival probability") +theme_bw()
  return(p)
} 

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}


## Posterior uncertainty intervals
mcmc_post_ci<-function(model,prob=0.95,nvar=5){
  ## Posterior uncertainty intervals
  hr.coef<-(rstanarm::fixef(model))
  hr.ci<-posterior_interval(model,prob = prob)
  hr.table=data.frame(mean.exp=exp(hr.coef),
                      ci.lb=exp(hr.ci[,1]),
                      ci.ub=exp(hr.ci[,2]))
  #rmarkdown to excel
  hr.table<-hr.table %>% 
    slice(2:nvar) %>% 
    mutate(HR=paste0(round(mean.exp,2),
                     " (",round(ci.lb,2),
                     " to ", round(ci.ub,2),
                     ")" ))%>% 
    dplyr::select(HR)
  return(hr.table)
}

## cross-validation for MCMC models
cross_validation<-function(priormodel,newdata,formula,id,niter=10,seed=123,knots=10){
  patients<-newdata %>% 
    group_by(!!sym(id)) %>% summarize()
  ## set the seed to make your partition reproducible
  set.seed(seed)
  cox_test <- list()
  cox_train <- list()
  exp_test <- list()
  exp_train <- list()
  ms_test <- list()
  ms_train <- list()
  print(formula)
  i<-1
  niter=niter
  while(i<=niter){
    #shuffling
    df <- sample_n(patients, nrow(patients), replace = TRUE)
    smp_size <- floor(0.8 * nrow(df))
    
    #splitting indices
    train_ind <- sample(seq_len(nrow(df)), size = smp_size)
    train_p <- df[train_ind, ]
    test_p <- df[-train_ind, ]
    #filtering X and Targets selected
    train<-newdata %>% 
      dplyr::filter(PAT_MRN_ID %in% train_p$PAT_MRN_ID)
    y_trainsurv <- Surv(train$tstart,
                        train$tstop,
                        train$edvisit)
    test<-newdata %>% 
      dplyr::filter(PAT_MRN_ID %in% test_p$PAT_MRN_ID)
    
    y_testsurv <- Surv(test$tstart,
                       test$tstop,
                       test$edvisit)
    
    coxmd<-coxph(formula = formula,data =train)
    cox_train[[i]] <- concordance(coxmd)$concordance %>% round(3)
    
    cox_test[[i]] <- concordance(y_testsurv~predict(coxmd,newdata = test,"lp"),
                                 reverse=T)$concordance %>% round(3)
    expmd <- update(priormodel,
                    data = train ,
                    prior_PD = FALSE,
                    basehaz = "exp")
    train$exploghaz<-posterior_survfit(expmd,
                                       newdata = train,
                                       extrapolate = F,
                                       type="loghaz",
                                       draws = 500,
                                       return_matrix = F,
                                       times       = "tstart",
                                       last_time   = "tstop")$median
    test$exploghaz<-posterior_survfit(expmd,
                                      newdata = test,
                                      extrapolate = F,
                                      type="loghaz",
                                      draws = 500,
                                      return_matrix = F,
                                      times       = "tstart",
                                      last_time   = "tstop")$median
    exp_train[[i]] <- concordance(y_trainsurv~train$exploghaz,
                                  data=train,
                                  reverse=T)$concordance %>% round(3)
    exp_test[[i]] <- concordance(y_testsurv~test$exploghaz,
                                 data=test,
                                 reverse=T)$concordance %>% round(3)
    
    msmd <- update(priormodel,
                   data = train ,
                   prior_PD = FALSE,
                   basehaz = "ms",
                   basehaz_ops = list(df = knots))
    train$msloghaz<-posterior_survfit(msmd,
                                      newdata = train,
                                      extrapolate = F,
                                      type="loghaz",
                                      draws = 500,
                                      return_matrix = F,
                                      times       = "tstart",
                                      last_time   = "tstop")$median
    test$msloghaz<-posterior_survfit(msmd,
                                     newdata = test,
                                     extrapolate = F,
                                     type="loghaz",
                                     draws = 500,
                                     return_matrix = F,
                                     times       = "tstart",
                                     last_time   = "tstop")$median
    ms_train[[i]] <- concordance(y_trainsurv~train$msloghaz,
                                 data=train,
                                 reverse=T)$concordance %>% round(3)
    ms_test[[i]] <- concordance(y_testsurv~test$msloghaz,
                                data=test,
                                reverse=T)$concordance %>% round(3)
    
    # hr.table<-mcmc_post_ci(msmd) 
    # print(hr.table)
    i<-i+1
    print(paste0("iter:",i))
  }
  
  return(list("cox_train" = unlist(cox_train),
              "cox_test" = unlist(cox_test),
              "exp_train" = unlist(exp_train),
              "exp_test" = unlist(exp_test),
              "ms_train" = unlist(ms_train),
              "ms_test" = unlist(ms_test))
  )
}


#### matrix reshape for plotting disease progression
matrix2column<-function(mat){
  if(ncol(mat)>1){
    # Use the as.vector function to convert the matrix to a vector
    v <- as.vector(mat)
    # Use the matrix function to convert the vector to a one column matrix
    m1 <- matrix(v, ncol = 1)
    # Print the reshaped matrix
    # print(m1)
    colnames(m1)<-"fromto"
    return(m1)
  }else
    colnames(mat)<-"fromto"
  return(mat)
}

# proportion matrix as a function
get_prop_matrix <- function(m,margin=NULL) {
  mat.c <- matrix(0, ncol(m), ncol(m))
  mat.p <- matrix(0, ncol(m), ncol(m))
  mat.j <- matrix(0, ncol(m), ncol(m))
  for (i in 1:ncol(m)) {
    for (j in i:ncol(m)) {
      # print(paste("i:", i, colnames(m)[i], ", j:", j, colnames(m)[j]))
      my_table_0 <- table(m[, i], m[, j])
      # print.table(my_table_0)
      my_table_2 <- round(prop.table(my_table_0, margin = margin), 2)
      # # have a look at the table
      # print.table(my_table_2)
      if (mat.p[i, j] == 0.0) {
        # counts
        mat.c[i, j] <- my_table_0[2, 2]
        # proportion
        mat.p[i, j] <- my_table_2[2, 2]
        # jaccard index
        # mat.j[i, j] <- my_table_2[2, 2]/(my_table_2[2, 1]+my_table_2[1, 2]-my_table_2[2, 2])
        mat.j[i, j] <- round(jaccard::jaccard(m[, i], m[, j]),2)
        if (mat.p[j, i] == 0.0) {
          mat.c[j, i] <- my_table_0[2, 2]
          mat.p[j, i] <- my_table_2[2, 2]
          # mat.j[j, i] <- my_table_2[2, 2]/(my_table_2[2, 1]+my_table_2[1, 2]-my_table_2[2, 2])
          mat.j[j, i] <- round(jaccard::jaccard(m[, i], m[, j]),2)
          
        }
      }
      
    }
  }
  diag(mat.j)<-1.0
  #print.table(my_table_0)
  return(list(
    mat.counts = mat.c,
    mat.prop = mat.p,
    mat.jaccard=mat.j
  ))
}

#############3 Get upper triangle of the correlation/proportion matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

########### required for the ggplot risk table
RiskSetCount <- function(timeindex, patientsdf,strataoi) {
  atrisk <- NULL
  survivaltime<-patients$tstop[patients$Strata==strataoi]
  for (t in timeindex)
    atrisk <- c(atrisk, sum(survivaltime >= t))
  df<-data.frame(strata=strataoi,value=atrisk,time=timeindex)
  return(df)
}
