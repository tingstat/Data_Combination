Gen_nextStats <- function(data, nS){
  
  n_time <- dim(data)[1]
  num_time <- length(unique(data$time))
  
  S <-as.matrix(data[sapply(1:nS, function(i) paste("S", i, sep = ""))], n_time, nS)
  
  NextS <- rbind(S[-1,drop=F ], rep(0, nS) )
  
  last_time_position <- which( data$time == num_time)
  
  NextS[last_time_position, ] =0
  
  data <- cbind(data, NextS)
  
  colnames(data)[(dim(data)[2]-nS+1):dim(data)[2] ] <- paste("next.S", 1:nS, sep = "")
  
  return(data)
}


phi_basis <- function(X){
  nX <- dim(X)[2]
  
  phi_vector <- NULL
  
  for(i in 1:nX){
    phi_vector <- cbind( phi_vector , cbind(X[,i], X[, i]^2, X[, i]^3  ) )
  }
  
  return(phi_vector)
}

library(dplyr)


behavir_est <- function(data, treatment, method='All'){
  nS <- length(grep("S", colnames(data)))/2
  covariate_names<- paste0("S", 1:nS)
  formula <- as.formula(paste("action ~", paste( covariate_names, sep="", collapse=" + ")))
  
  if(method=='All'){
    model <- glm(formula, data = data, family = binomial)
    pre_prob <- predict(model, type = "response")
    data$pi_t <- data$action*pre_prob+(1-data$action)*(1-pre_prob)

    new_data <- data

    }else{

    data_split <- split(data, data$time)
  
    new_data <- do.call(rbind, lapply(data_split, function(subdata) {
    model <- glm(formula, data = subdata, family = binomial)
    pre_prob <- predict(model, type = "response")
    subdata$pi_t <- subdata$action*pre_prob+(1-subdata$action)*(1-pre_prob)
    subdata
  }))
  
  new_data <- new_data[order(new_data$user),]
  }
  

  return(new_data)
}



library(mgcv)
######
mu_t_est <- function(data, treatment){
  nS <- length(grep("S", colnames(data)))/2
  
  time_num <- length(unique(data$time))
  n_num <- length(unique(data$user))
  
  init_position <- which(data$time==1)
  
  data$omega_t <- 0
  data$omega_t[init_position]=1
  
  
  lambda_t <- aggregate(eta_t~user, data=data, FUN=cumprod)
  data$lambda_t <- as.numeric(t(lambda_t$eta_t))
  lambda_t <- data$lambda_t
  
  lambda_t_before <- c(0, lambda_t[-length(lambda_t)])
  
  lambda_t_response <- lambda_t_before[-which(data$time==1)]
  
  #lambda_t_position <- which(data$lambda_t>0)
  
  data_select <- data #[lambda_t_position,]
  
  n_time <- dim(data_select)[1]
  
  S <-as.matrix(data_select[sapply(1:nS, function(i) paste("S", i, sep = ""))], n_time, nS)
  next.S<- as.matrix(data_select[sapply(1:nS, function(i) paste("next.S", i, sep = ""))], n_time, nS)
  
  S <- S[-which(data$time==1),]
  
  covariate_names<- paste0("S", 1:nS)
  formula <- as.formula(paste("lambda_t_response ~", paste("s(", covariate_names, ")", sep="", collapse=" + ")))
  
  data_for_gam <- cbind(lambda_t_response, S)
  colnames(data_for_gam) <- c('lambda_t_response', covariate_names )
  data_for_gam <- as.data.frame(data_for_gam)
  
  lm_mu <- gam(  formula, data=data_for_gam, method='GCV.Cp' )
  
  omega_t_est <- lm_mu$fitted.value
  
  omega_t_est[which(omega_t_est<0)]=0
  
  data$omega_t[-init_position] <- omega_t_est 
  
  
  #data$mu_t <- 0
  mu_t_est <- data$eta_t*data$omega_t
  
  mu_t_est[which(mu_t_est<0)]=0
  #mu_t_est <- mu_t_est/mean(mu_t_est)
  
  #####
  mu_t_matrix <- matrix(mu_t_est, ncol=n_num)
  mu_t_mean <- mu_t_est/apply(mu_t_matrix, 1, mean)
  
  data$mu_t <- mu_t_mean
  
  mu_t_est <- mu_t_mean
  
  list(mu_t_est=mu_t_est, omega_t_est=data$omega_t,phi_S_eta_mu_all=NULL)
  
}


mu_t_est_backward <- function(data, treatment){
  nS <- length(grep("S", colnames(data)))/2
  
  time_num <- length(unique(data$time))
  n_num <- length(unique(data$user))
  n_time <- dim(data)[1]
  
  mu_t_est <- rep(0, n_time)
  
  
  S <-as.matrix(data[sapply(1:nS, function(i) paste("S", i, sep = ""))], n_time, nS)
  next.S<- as.matrix(data[sapply(1:nS, function(i) paste("next.S", i, sep = ""))], n_time, nS)
  
  if(Q_indicator=='Linear'){
    phi_S <- S #phi_basis(S)
    phi_next.S <- next.S  #phi_basis(next.S)
  }else{
    phi_S <- phi_basis(S)
    phi_next.S <- phi_basis(next.S)
  }
  
  phi_S_eta <- as.matrix(phi_S*data$eta_t)
  
  position <- which(data$action == treatment) 
  position_1 <- which(data$action == 1) 
  
  gamma_1 <- tcrossprod( solve( crossprod( phi_S [position, ] ) ), t(colSums(phi_S_eta)) )
  mu_k_est <- crossprod( t(phi_S), gamma_1)
  
  mu_t_est[intersect(position_1,position) ] <- mu_k_est[intersect(position_1,position)] 
  
  phi_S_eta_mu_all <- colSums(phi_S_eta) # for estimate omega_t_est
  
  for(k in c(2:time_num)){
    position_greater_than_k <- which(data$time >=k)
    position_temp <- intersect(position, position_greater_than_k)
    
    phi_S_eta_mu <-as.matrix(  phi_S[position_temp,]*data$eta_t[position_temp]*mu_k_est[position_temp-1] )
    
    gamma_k <- tcrossprod( solve( crossprod( phi_S [position_temp, ] ) ), t(colSums(phi_S_eta_mu) ) )
    
    mu_k_est <- crossprod( t(phi_S), gamma_k)
    
    position_k <- which(data$time ==k)
    mu_t_est[intersect(position_k,position) ] <- mu_k_est[intersect(position_k,position)] 
    
    phi_S_eta_mu_all <- cbind(phi_S_eta_mu_all, colSums(phi_S_eta_mu ) )
    
  }
  
  
  ##### check positivity and normalization ##########
  mu_t_est[which(mu_t_est<0)]=0

  mu_t_matrix <- matrix(mu_t_est, ncol=n_num)
  mu_t_est <- mu_t_est/apply(mu_t_matrix, 1, mean)
  

  
  list(mu_t_est=mu_t_est, phi_S_eta_mu_all=phi_S_eta_mu_all)
}

omega_t_est_backward <- function(data0, phi_S_eta_mu_all,n_D, n_H){
  
  nS <- length(grep("S", colnames(data0)))/2
  
  time_num <- length(unique(data0$time))
  n_num <- length(unique(data0$user))
  n_time <- dim(data0)[1]
  
  omega_t_est <- rep(0, n_time)
  
  
  S <-as.matrix(data0[sapply(1:nS, function(i) paste("S", i, sep = ""))], n_time, nS)
  next.S<- as.matrix(data0[sapply(1:nS, function(i) paste("next.S", i, sep = ""))], n_time, nS)
  
  if(Q_indicator=='Linear'){
    phi_S <- S #phi_basis(S)
    phi_next.S <- next.S  #phi_basis(next.S)
  }else{
    phi_S <- phi_basis(S)
    phi_next.S <- phi_basis(next.S)
  }
  
  for(k in c(1:time_num)){
    position_temp  < - which(data$time >=k)
    
    gamma_k <- tcrossprod( solve( crossprod( phi_S [position_temp, ] ) ), t(phi_S_eta_mu_all[,k] ) )*n_H/n_D
    
    omega_k_est <- crossprod( t(phi_S), gamma_k)
    
    position_k <- which(data$time ==k)
    omega_t_est[position_k ] <- omega_k_est[position_k] 

  }
  
  ##### check positivity and normalization ##########
  omega_t_est[which(omega_t_est<0)]=0
  
  omega_t_matrix <- matrix(omega_t_est, ncol=n_num)
  omega_t_est<- omega_t_est/apply(omega_t_matrix, 1, mean)
  
  return(omega_t_est)
  
}


Q_eta_est_backward <- function(data, treatment, prob_behavior){
  if(is.na(prob_behavior)){
    data <- behavir_est(data,treatment) 
  }else{
    data$pi_t <- as.numeric(data$action==1)*prob_behavior+ as.numeric(data$action==0)* (1 - prob_behavior)
  }
  
  eta_t <- as.numeric(data$action==treatment)/data$pi_t
  data$eta_t <- eta_t
  
  nS <- length(grep("S", colnames(data)))/2
  
  
  position <- which(data$action == treatment)
  nT_all <- dim(data)[1]
  
  n_num <- length(unique(data$user))
  time_num <- length(unique(data$time))
  
  data_select <- data[position, ]
  
  action <- data$action[position]
  reward <- data$reward[position]
  
  n_time <- dim(data)[1]
  
  S <-as.matrix(data[sapply(1:nS, function(i) paste("S", i, sep = ""))], n_time, nS)
  next.S<- as.matrix(data[sapply(1:nS, function(i) paste("next.S", i, sep = ""))], n_time, nS)
  
 

  if(Q_indicator=='Linear'){
    phi_S <- S #phi_basis(S)
    phi_next.S <- next.S  #phi_basis(next.S)
    }else{
      phi_S <- phi_basis(S)
      phi_next.S <- phi_basis(next.S)
    }
  
  ###k=1####################
  position_k <- which(  data$time ==(time_num)  ) # 
  position_tmp <- intersect( position, position_k)
  
  fitted_k  <- lm(data$reward[position_tmp]  ~ phi_S[position_tmp,] )
  fitted_value <- cbind(1, phi_S )%*%fitted_k$coe
  
  fitted_q_hat <- fitted_value
  fitted_q_next <- rep(0, length(fitted_value))
  
  
  for(k in c(1:(time_num-1) )){
    #print(k)
    if(time_num==1){break}
    
    position_k <- which(  data$time <=(time_num-k)  ) # let t*=t-k, the position of t*
    position_tmp <- intersect( position, position_k)
    
    V_hat_temp <- fitted_value[position_tmp + 1] 
    fitted_q_next[position_tmp] <- fitted_value[position_tmp+1]
    
    phi_S_temp <- phi_S[position_tmp,] 
    
    reward_temp <- data$reward[position_tmp] + V_hat_temp # 
    
    fitted_k  <- lm( reward_temp~ phi_S_temp )
    fitted_value <- cbind(1, phi_S )%*%fitted_k$coe
    #print(mean(fitted_k$res) )
    
    
    #position_t <- which(  data$time ==(time_num-k)  ) 
    
    fitted_q_hat[position_tmp ] <-  fitted_value[position_tmp]
    
  }
  

  if(time_num >1){
    #mu_t_estimation <- mu_t_est_backward(data, treatment)
    mu_t_estimation <- mu_t_est(data, treatment)
    data$mu_t <- mu_t_estimation$mu_t_est
    phi_S_eta_mu_all <- mu_t_estimation$phi_S_eta_mu_all
  }else{
    data$mu_t <- data$eta_t
    phi_S_eta_mu_all<- NULL
  }


  ###############################
  
  data$fitted_q_hat <- as.numeric(fitted_q_hat)
  data$fitted_q_next <- fitted_q_next
  
  data$eta <-  data$mu_t*(data$reward +data$fitted_q_next - data$fitted_q_hat)
  data$eta <- data$eta/time_num
  
  data$fitted_value <- fitted_value
  
  value0 <- mean(fitted_value[position_k] )/time_num
  eta_est <- value0 +  mean( aggregate( eta~user, data=data, FUN=sum)$eta) 
  
  TD_user <- fitted_value[position_k]/time_num + aggregate( eta~user, data=data, FUN=sum)$eta - eta_est 
  
  variance <- mean( TD_user^2 )
  TD_user_value0  <-fitted_value[position_k]/time_num -value0
  var_value0 <- mean(TD_user_value0^2 )
  
  ##### calculate variance of the value function #####
  cov_coe=vcov(fitted_k)
  variable_temp <- cbind(1, phi_S )[position_k,]/time_num
  E_var <- mean( diag( tcrossprod( crossprod(t(variable_temp), cov_coe),variable_temp )) )
  Var_E <- tcrossprod( crossprod(fitted_k$coe, cov(variable_temp)), t(fitted_k$coe) )
  
  var_v0S <- E_var+Var_E 
  
  
  list( eta_est = eta_est, value0=value0,  var_value0= var_value0, TD_user_value0=TD_user_value0, 
        TD_user=TD_user, variance=variance, data=data,cov_coe=vcov(fitted_k), last_coe= fitted_k$coe, 
        last_res=fitted_k$res, variable = cbind(1, phi_S ),var_v0S=var_v0S, phi_S_eta_mu_all=phi_S_eta_mu_all)
  
}

eta_est_historical  <- function(result_D, result_H, n_D,n_H, time_num,ratio_indicator='Given'){
  
  position_1_D <- which(result_D$data$time==1)
  fitted_value_for_D <- (result_D$variable%*%result_H$last_coe)[position_1_D]
  
    eta_est_for_D_part_1 <-mean(fitted_value_for_D )/time_num 
  eta_est_for_D_part_2 <- mean( aggregate( eta~user, data=result_H$data, FUN=sum)$eta)
  
  
  if(ratio_indicator=='Given'){
    omega_t_est <- rep(1, n_H*time_num)
  }else{
    omega_t_est <- omega_t_est_backward(result_H$data, result_D_0$phi_S_eta_mu_all,n_D, n_H)
    }
  
  #TD_user <- aggregate( eta~user, data=data, FUN=sum)$eta - mean( aggregate( eta~user, data=data, FUN=sum)$eta) 
  TD_user_for_D_part_1 <- fitted_value_for_D/time_num  - eta_est_for_D_part_1
  
  result_H$data$rho_t <-  omega_t_est*(result_H$data$reward +result_H$data$fitted_q_next - result_H$data$fitted_q_hat)/time_num

  
  TD_user_for_D_part_2 <- aggregate( eta~user, data=result_H$data, FUN=sum)$eta - eta_est_for_D_part_2
  
  variance_for_D_part_1 <- mean( TD_user_for_D_part_1^2 ) 
  variance_for_D_part_2 <- mean( TD_user_for_D_part_2^2 ) 
  cov_D_H <- mean( TD_user_for_D_part_1*result_D$TD_user )/n_D
  
  eta_H_for_D <- eta_est_for_D_part_1  + eta_est_for_D_part_2
  
  var_H_for_D <- variance_for_D_part_1/n_D +  variance_for_D_part_2/n_H
  
  list(eta_est=eta_H_for_D, var_H_for_D=var_H_for_D, TD_user_for_D_part_1=TD_user_for_D_part_1  )
}



############## pessimistic combine ########
pessimistic_combine <- function(result_D_1, result_D_0, result_H, n_D, n_H,time_num,ratio_indicator='Given'){
  eta_D_1 <- result_D_1$eta_est
  eta_D_0 <- result_D_0$eta_est
  eta_H <- result_H$eta_est
  eta_est_hist <- eta_est_historical(result_D_0, result_H, n_D,n_H, time_num,ratio_indicator)
  
  eta_H_for_D <- eta_est_hist$eta_est
  
  ATE_e <- eta_D_1 - eta_D_0
  ATE_h <- eta_D_1  - eta_H_for_D
  
  cov_D_1_0 <- mean(result_D_1$TD_user*result_D_0$TD_user )
  cov_D_1_H <- mean(result_D_1$TD_user*eta_est_hist$TD_user_for_D_part_1 )
  
  var_D <- (result_D_1$variance + result_D_0$variance - 2 *cov_D_1_0)/n_D
  var_H <- (result_D_1$variance + - 2 *cov_D_1_H)/n_D + eta_est_hist $var_H_for_D 
  cov_D_H <- mean( (result_D_1$TD_user-result_D_0$TD_user)*(result_D_1$TD_user-eta_est_hist$TD_user_for_D_part_1) )/n_D
  
  bias_square_UB_for_D <- ( abs(eta_D_0 - eta_H_for_D) + 1.64*sqrt(var_D+var_H-2*cov_D_H) )^2
  
  weight_for_D <- (var_H+ bias_square_UB_for_D - cov_D_H  )/(var_D  + var_H+ bias_square_UB_for_D - 2*cov_D_H )
  
  ATE_pessi_DR<- weight_for_D*ATE_e+ (1-weight_for_D)*ATE_h
  
  
  #### minimize MSE #####
  bias_square <- (eta_D_0 - eta_H_for_D)^2
  weight_1 <- (var_H + bias_square - cov_D_H )/(var_D  + var_H +bias_square- 2*cov_D_H)
  ATE_mse<- weight_1*ATE_e+ (1-weight_1)*ATE_h
  
  #### L1 penalty #####
  weight_2 <- (var_D - 0.4*bias_square)/(var_D  + var_H )
  weight_2 <- ifelse(weight_2 > 0 & weight_2 < 1, weight_2, 0)
  
  ATE_L1<- (1-weight_2)*ATE_e+ weight_2*ATE_h
  
  
  #### Based on Testing ########
  tstat <- (eta_D_0 - eta_H_for_D)/sqrt( var_D  + var_H-2*cov_D_H )
  
  if(abs(tstat) > qnorm(0.975) ){
    ATE_test005 <- ATE_e
  }else{
    ATE_test005 <-ATE_pessi_DR
  }
  
  
  ATEs <- c(ATE_e, ATE_h, ATE_pessi_DR, ATE_mse, ATE_L1)
  
  return(ATEs)
}



eta_est_historical_value  <- function(result_D, result_H, n_D,n_H, time_num,ratio_indicator='Given'){
  
  position_1_D <- which(result_D$data$time==1)
  fitted_value_for_D <- (result_D$variable%*%result_H$last_coe)[position_1_D]
  
  eta_est_for_D_part_1 <-mean(fitted_value_for_D )/time_num 

  TD_user_for_D_part_1 <- fitted_value_for_D/time_num  - eta_est_for_D_part_1
  
  
  variance_for_D_part_1 <- mean( TD_user_for_D_part_1^2 ) 

  
  eta_H_for_D <- eta_est_for_D_part_1  
  
  var_H_for_D <- variance_for_D_part_1/n_D 
  
  list(eta_est=eta_H_for_D, var_H_for_D=var_H_for_D, TD_user_for_D_part_1=TD_user_for_D_part_1,fitted_value_for_D= fitted_value_for_D )
}


############## pessimistic combine ########
pessimistic_combine_value <- function(result_D_1, result_D_0, result_H, n_D, n_H,time_num,ratio_indicator='Given'){
  eta_D_1 <- result_D_1$value0
  eta_D_0 <- result_D_0$value0
  eta_H <- result_H$value0
  eta_est_hist <- eta_est_historical_value(result_D_0, result_H, n_D,n_H, time_num,ratio_indicator)
  
  eta_H_for_D <- eta_est_hist$eta_est
  
  ATE_e <- eta_D_1 - eta_D_0
  ATE_h <- eta_D_1  - eta_H_for_D
  
  cov_D_1_0 <- mean(result_D_1$TD_user_value0 *result_D_0$TD_user_value0  )
  cov_D_1_H <- mean(result_D_1$TD_user_value0 *eta_est_hist$TD_user_for_D_part_1  )
  
  var_D <- (result_D_1$var_value0 + result_D_0$var_value0 - 2 *cov_D_1_0)/n_D
  var_H <- (result_D_1$var_value0 - 2 *cov_D_1_H)/n_D + eta_est_hist $var_H_for_D 
  cov_D_H <- mean( (result_D_1$TD_user_value0-result_D_0$TD_user_value0)*(result_D_1$TD_user_value0-eta_est_hist$TD_user_for_D_part_1) )/n_D
  
  bias_square_UB_for_D <- ( abs(eta_D_0 - eta_H_for_D) + 1.64*sqrt(var_D+var_H-2*cov_D_H) )^2
  
  weight_for_D <- (var_H+ bias_square_UB_for_D - cov_D_H  )/(var_D  + var_H+ bias_square_UB_for_D - 2*cov_D_H )
  
  ATE_pessi_DM<- weight_for_D*ATE_e+ (1-weight_for_D)*ATE_h
  
  
  #### minimize MSE #####
  bias_square <- (eta_D_0 - eta_H_for_D)^2
  weight_1 <- (var_H + bias_square - cov_D_H )/(var_D  + var_H +bias_square- 2*cov_D_H)
  ATE_mse<- weight_1*ATE_e+ (1-weight_1)*ATE_h
  
  #### L1 penalty #####
  weight_2 <- (var_D - 0.4*bias_square)/(var_D  + var_H )
  weight_2 <- ifelse(weight_2 > 0 & weight_2 < 1, weight_2, 0)
  
  ATE_L1<- (1-weight_2)*ATE_e+ weight_2*ATE_h
  
  
  #### Based on Testing ########
  tstat <- (eta_D_0 - eta_H_for_D)/sqrt( var_D  + var_H-2*cov_D_H )
  
  if(abs(tstat) > qnorm(0.95) ){
    ATE_test005 <- ATE_e
  }else{
    ATE_test005 <-ATE_pessi_DM
  }

    ##### depend on S #####
  position_1_D <- which(result_D_0$data$time==1) 

  fitted_value_D <- (result_D_0$variable%*%result_D_0$last_coe)[position_1_D]
  fitted_value_for_D <- eta_est_hist$fitted_value_for_D


  cov_D <-  (diag( tcrossprod( crossprod(t(result_D_0$variable), result_D_0$cov_coe), result_D_0$variable) ) )[position_1_D]/(time_num^2)
  
  cov_H <- ( diag( tcrossprod( crossprod(t(result_D_0$variable), result_H$cov_coe), result_D_0$variable) ) )[position_1_D]/(time_num^2)
  

  bias_square_UB_for_D_S <- ( abs(fitted_value_D - fitted_value_for_D)/time_num + 1.64*sqrt(cov_D+cov_H) )^2
  
  weight_for_D_S <- (cov_H + bias_square_UB_for_D_S )/(cov_D + cov_H + bias_square_UB_for_D_S )
  
  ATE_pessi_S<-eta_D_1- mean( weight_for_D_S*fitted_value_D+ (1-weight_for_D_S)*fitted_value_for_D )/time_num
  
  
  ATEs <- c(ATE_e, ATE_h, ATE_pessi_DM,ATE_pessi_S, ATE_mse, ATE_L1,ATE_test005)
  
  return(ATEs)
}



################################## nS=1#################
#################
######## reward ###########
rwrd.gen = function(S, a, S.next,error_tmp,t,d) {
  R <- 10 +d*mu_diff+ b*a + b * S  + (2 + 1*delta*d )*error_tmp 
  return(R)
}


### next state #####
next.S.gen = function(S, a,t,d) {
  
  S.next <- 0.25*S + 0.2*a + (2+d*p_s)*rnorm(1) #0.08*t+
  #S.next <- phi_0[t]+phi_1[t]*S +0.08*t+ a*alpha[t]+ (1+a*p_s)*rnorm(1,0, sqrt(0.85))
  
  return(S.next)
}

######### reward ###########
#rwrd.gen = function(S, a, S.next,error_tmp) {
#10  + 2*a + 0.25 * S * a * (0.04 + 0.02*S) + (0.4 + delta*a )*rnorm(1) # (0.16 + 0.8*a )*rnorm(1)  #0.5*a
# 10  + 2*a + 2 * S  + (0.4 + delta*a )*error_tmp # (0.16 + 0.8*a )*rnorm(1)  #0.5*a

#}
### next state #####
#next.S.gen = function(S, a) {

# S.next <- 0.25*S + 2*a + (2+a*p_s)*rnorm(1)

#return(S.next)
#}

Data_gen <- function(n, nT, prob,d,TI=1){
  
  data <- NULL;
  error_sigma <- 0.5^( abs( outer(1:nT, 1:nT, "-") ) )
  #error_sigma <- exp( -abs( outer(1:nT, 1:nT, "-") )/(2*nT) )
  eps_cors1=diag(x=rep(1,nT),nrow = nT)
  
  errors <- mvrnorm(n, rep(0, nT), error_sigma ) #+ mvrnorm(n, rep(0, nT), eps_cors1 )
  
  if(TI==nT){
    A1 <- rep(1,TI)
    A2<- rep(0, TI)
  }else{
    A1 <- rep( rep(c(0, 1), each=TI), nT/TI/2)
    A2 <- rep( rep(c(1, 0), each=TI), nT/TI/2)
  }
  
  
  
  for(i in 1:n){
    
    next.S <- rnorm(1)
    
    for(t in 1:nT){
      
      S <- next.S;
      
      ###
      if(is.na(prob)){
        if(i%%2==0){
          A= as.numeric(t%%2==0)  ## alternation design
        }else{
          A= as.numeric((t+1)%%2==0)
        }
        
      }else if(prob==3){
        if(i%%2==0){
          A =A1[t] #as.numeric(t%%TI==0)  #as.numeric( t <= (nT /2) )
        }else{
          A=A2[t] #as.numeric(t > (nT /2))
        }
        #print(A)
      }else if(prob==4){
        A = as.numeric( i <= (n /2) )
      }else{
        A <-rbinom(1,1,prob) #(runif(1) < prob);
      }
      
      ###
      next.S <- next.S.gen(S, A,t,d);
      R <- rwrd.gen(S, A,  next.S, errors[i,t],t ,d);
      
      
      data <- rbind(data, c(i, t, S, A, prob, R, next.S))
    }
    
  }
  
  colnames(data) <- c("user", "time",'S1', "action", "prob", "reward", 'next.S1')
  
  data <- data.frame(data)
  
  #data$next.S[which(data$time==nT)]=0
  
  return(data)
  
}



# #################
# ######### nS =3 ###########
# rwrd.gen = function(S, a, S.next, error_tmp,d) {
#   10  + b*a - 0*S[3] +b*S[2]+ b * S[1] * a * (0.04 + 0.02*S[1]) +  (1 + delta*d )*error_tmp #(0.16 + 0.8*a )*rnorm(1)(for T 1000)   #0.5*a
# }
# 
# ##### error term with AR structure ######
# library(MASS)
# 
# 
# 
# ### next state #####
# next.S.gen = function(S, a,t,d) {
# 
#   S1 <- 0.5*S[1] +  2*rnorm(1)  #0.08*t
#   S2 <- 0.25*S[2] + 0.2*a + (2+ p_s*d )*rnorm(1)  #0.05*t
#   S3 <- 0.6*S[3] +  0.05*a*S[3] + 0.5*a + rnorm(1) #0.06*t
# 
#   S.next <- c(S1, S2, S3)
# 
#   return(S.next)
# }
# 
# Data_gen <- function(n, nT, prob, d,nS=3){
# 
#   error_sigma <- 0.5^( abs( outer(1:nT, 1:nT, "-") ) )
# 
#   #error_sigma <- exp( -abs( outer(1:nT, 1:nT, "-") )/(2*nT) )
# 
#   eps_cors1=diag(x=rep(1,nT),nrow = nT)
# 
#   errors <- mvrnorm(n, rep(0, nT),eps_cors1 )
# 
#   data <- NULL;
#   for(i in 1:n){
# 
#     next.S <- rnorm(nS)
# 
#     for(t in 1:nT){
# 
#       S <- next.S;
# 
#       ###
#       if(is.na(prob)){
#         if(i%%2==0){
#           A= as.numeric(t%%2==0)  ## alternation design
#         }else{
#           A= as.numeric((t+1)%%2==0)
#         }
# 
#       }else if(prob==3){
#         A = as.numeric( i <= (n /2) )
#         #print(A)
#       }else{
#         A <- (runif(1) < prob);
#       }
# 
#       ###
#       next.S <- next.S.gen(S, A,t,d);
#       error_tmp <- errors[i, t]
#       R <- rwrd.gen(S,A,  next.S,error_tmp,d);
# 
# 
#       data <- rbind(data, c(i, t, S, A, prob, R, next.S))
#     }
# 
#   }
# 
#   colnames(data) <- c("user", "time",paste("S", 1:nS, sep = ""), "action", "prob", "reward", paste("next.S", 1:nS, sep = ""))
# 
#   data <- data.frame(data)
#   #data$next.S[which(data$time==nT)]=0
# 
#   return(data.frame(data))
# 
# }





