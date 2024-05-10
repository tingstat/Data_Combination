library(MASS)
source('Est_CombiningData.R')


multi <- 3
mu_diff <- 0
prob_behavior <- NA
delta <- 1
print(c(multi ,mu_diff, prob_behavior, delta))


n <- 48
n0 <- multi*n

nn <- 1
p_s <-0
b=0 # 1

t0 <- Sys.time()

######## check the distribution of direct method ##########
Q_indicator <- 'Linear'
method <- 'Simple'
ratio_indicator='Given'

################################################
simu_times <- 200

data_100_1000_1 <- Data_gen(5000, nn, 1,d=1)
data_100_1000_0 <- Data_gen(5000, nn, 0,d=1)

ATE_emprical_true <-  mean(data_100_1000_1$reward) - mean(data_100_1000_0$reward) #10.777, 9.26
ATE_true_1 <- mean(data_100_1000_1$reward)
ATE_true_0 <- mean(data_100_1000_0$reward)

#simu <- function(simu_times){
 
  ATEs_all <- NULL
  
  for(i in 1:simu_times){
    print(i)
    
    data <- Data_gen(n, nn, prob_behavior, d=1)  ## 1/2
    data0 <- Data_gen(n0, nn, 0,d=0) 
    
    ATE_1_temp <- try( Q_eta_est_backward(data, 1, prob_behavior),silent=T) #Q_eta_est_poly
    ATE_0_temp_D <- try( Q_eta_est_backward(data, 0, prob_behavior), silent=T)
    ATE_0_temp_H <- try( Q_eta_est_backward(data0, 0, prob_behavior), silent=T)
    
    ATEs <- try( pessimistic_combine(ATE_1_temp, ATE_0_temp_D, ATE_0_temp_H,n,n0,nn,ratio_indicator='Given'), silent=T)

     if('try-error'%in%class(ATE_1_temp) |  'try-error'%in%class(ATE_0_temp_D) |  'try-error'%in%class(ATE_0_temp_H) |  'try-error'%in%class(ATEs)  ){next}

    
   
    ATEs_all <- rbind(ATEs_all,c(ATEs, ATE_emprical_true) )
  }
  
colnames(ATEs_all) <- c('ATE_e', 'ATE_h', 'ATE_pessi_DR', 'ATE_mse', 'ATE_L1',  'True')
  
  MSEs_ATE <- colMeans( (ATEs_all  - ATE_emprical_true)^2 )

  
  names(MSEs_ATE) <- c('ATE_e', 'ATE_h', 'ATE_pessi_DR', 'ATE_mse', 'ATE_L1',
                        'True')

  
  #list(MSEs_ATE=MSEs_ATE,  ATEs_all= ATEs_all )
  
#}


MSEs_results<- list(MSEs_ATE=MSEs_ATE,  ATEs_all= ATEs_all )
print(Sys.time()-t0)

print(MSEs_results$MSEs)

dir.create(paste0('./', 'Results_n',n, '_nn_', nn))

save(MSEs_results, file=paste0('./', 'Results_n',n, '_nn_', nn,'/Mses_different_method_multi_',multi, '_muDiff_', mu_diff,  'prob', prob_behavior, '_delta_', delta,'.Rdata') )
