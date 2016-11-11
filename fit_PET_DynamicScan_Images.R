# this script fits the PET images and calculates the permeability of the imaging tracer.
# it can also calculate the target concentration.

library(plotrix)
library(deSolve)
library(ggplot2)
library(optimx)
library(reshape2)
require(gridExtra)
library(RColorBrewer)

#######################
##    User  INPUT    ##
#######################
filepath = 'xxx.csv'
time_shift = 0 # min

ignore_dyna_points = 14
animal_ind = 211

# log10_fit = 'no'
log10_fit = 'yes'

MW = 100 # kDa
dose = 3 # mg/kg
mouse_weight = 20 # g
K_D = 30000000 # nM
v_e = 0.2
dose_amount_mg = dose * mouse_weight*1e-3 # mg in mouse
nM_converstion = (1/100)*dose_amount_mg/(MW*1000)*(1e9) # convert %ID/g to nM. MW is in kDa.

## calculate P and epsilon
source("/xxx/calculate_P_and_epsilon.R")
P_and_epsilon = calculate_P_and_epsilon(MW=MW,K_D=K_D) # P = cm/hr

## calculate blood PK
source("/xxx/calculate_blood_PK.R")
calculate_initial_blood_PK_parameters(MW=MW,dose=100,hematocrit=0.45) # CL, V_ss, V_c, CL_d, L_1, L_2, lambda_1, lambda_2
initial_blood_PK_parameters[c('L_1','L_2','lambda_1','lambda_2')] = c(15.025, 15.0001, 0.0045, 0.0041)

## set up initial parms for the fit
initial_parms=c(
  k_trans=6.69E-4, # cm/min
  v_p=0.02,
  epsilon=as.numeric(P_and_epsilon["epsilon"]),
  k_ns=1E-6,
  k_elim=1E-6,
  k_endo=1E-6,
  Ag = 0.01
)

# L_1, L_2, lambda_1, lambda_2
lower_bounds_blood <- c(0.00001, 0.00001, 0.00001,0.00001) 
upper_bounds_blood <- c(100,100,10,1)
# k_trans, v_p, epsilon, k_ns, k_elim, k_endo, Ag
lower_bounds_tumor <- c( 1E-7, 0.001, 0.001)
upper_bounds_tumor <- c( 0.9E-2, 0.5, 0.5)

#######################
## End of User INPUT ##
#######################

## transform the data
data_raw = read.csv(name_of_dataset)
colnames(data_raw)    <- c('animal','day','tumor_time','tumor','tumor_upper','tumor_lower','tumor_std','heart_time','heart','heart_upper','hear_lower','heart_std')
# shift the time points
data_raw['tumor_time']=data_raw['tumor_time'] + signif(time_shift, 3)
data_raw['heart_time']=data_raw['heart_time'] + signif(time_shift, 3)

data = data_raw[,c('animal','day','tumor_time','tumor','heart')]
data_sd = data_raw[,c('animal','day','tumor_time','tumor_std','heart_std')]
colnames(data_sd) = c('animal','day','tumor_time','tumor','heart')

data_ana = data[data['animal']==animal_ind,]
num_of_day = length(unique(data_ana$day))
num_of_animal = length(unique(data$animal))

# setting the colour palette
myColors_tumor <- brewer.pal(num_of_animal,"RdYlBu")
names(myColors_tumor) <- levels(factor(data_raw$animal))
colScale_tumor <- scale_colour_manual(name = "animal",values = myColors_tumor)

p_tumor=ggplot(data_ana, aes(x=tumor_time, y=tumor, colour = factor(day))) +
  theme_bw() + # shape = factor(animal)
  # geom_errorbar(aes(ymin=tumor-tumor_std, ymax=tumor+tumor_std),colour='skyblue2', width=0.5) +
  geom_point(size=2, alpha=0.7) +
  scale_x_continuous(limits=c(0, 60), breaks = seq(0, 65, 10), minor_breaks=NULL) +
  scale_y_continuous(limits=c(0.2, 4), breaks = seq(0, 10, 0.5), minor_breaks=NULL) +
  labs(x='minutes',y='Tumor %ID/g') + 
  theme(
    #axis.line=element_blank(),
    axis.line=element_line(size=0.5, colour="grey16"),
    axis.text=element_text(size=16, vjust=0.5),
    #axis.text.y=element_text(size=20, vjust=0.5),
    axis.title.x=element_text(size=20, vjust=0),
    axis.title.y=element_text(size=20, vjust=1),
    axis.title.y=element_blank(),
    panel.grid.major=element_line(colour="#D0D0D0",size=.5),
    legend.key = element_blank(),
    legend.title = element_text(colour="grey12", size=14),
    legend.text = element_text(colour="grey12", size = 14))
    #legend.position="none") 
    #legend.justification=c(1,0), legend.position=c(1,0)) 
  #colScale_tumor

p_heart=ggplot(data_ana, aes(x=tumor_time, y=heart, colour=factor(day))) + 
  theme_bw() + # shape = factor(animal)
  # geom_errorbar(aes(ymin=heart-heart_std, ymax=heart+heart_std), colour='firebrick2', width=0.5) +
  geom_point(size=2.5, alpha=0.7) +
  scale_x_continuous(limits=c(0, 60), breaks = seq(0, 65, 10), minor_breaks=NULL) +
  scale_y_continuous(limits=c(10, 38), breaks = seq(0, 62, 2), minor_breaks=NULL) +
  #scale_y_continuous(breaks = seq(0, 60, 2), minor_breaks=NULL) +
  labs(x='minutes',y='Heart %ID/g') +
  theme(
    #axis.line=element_blank(),
    axis.line=element_line(size=0.5, colour="grey16"),
    axis.text=element_text(size=16, vjust=0.5),
    #axis.text.y=element_text(size=20, vjust=0.5),
    axis.title.x=element_text(size=20, vjust=0),
    axis.title.y=element_text(size=20, vjust=1),
    panel.grid.major=element_line(colour="#D0D0D0",size=.5),
    legend.key = element_blank(),
    legend.title = element_text(colour="grey12", size=14),
    legend.text = element_text(colour="grey12", size = 14))
    #legend.justification=c(0,1), legend.position=c(0.5,1),
    #legend.justification=c(0,1), legend.position="none")
  #colScale_tumor

data_heart=dcast(data_ana, tumor_time ~ day, value.var='heart') # take the heart data from long to wide
data_heart=data_heart[ignore_dyna_points:length(data_heart$'tumor_time'),] # remove the first 2 data point
data_heart=melt(data_heart, id='tumor_time') # melt the data to long form
H_blood=list()
value_blood=list()
aa = lapply(1:num_of_day, function(x) {
  animal_day_tick = unique(data_ana$'day')[x]
  time_blood = data_heart[data_heart['variable']==animal_day_tick,]$tumor_time
  DV = data_heart[data_heart['variable']==animal_day_tick,]$value
  
  if (log10_fit == 'yes'){
    blood_PK_bo2_fit <- optim(par=log10(initial_blood_PK_parameters[c('L_1','L_2','lambda_1','lambda_2')]),fn=log10_objective_function_for_biexponential_blood_PK_model, lower=log10(lower_bounds_blood), time=time_blood,DV=DV,method='L-BFGS-B',hessian=T,control=list(maxit=1000))
    H_blood[[x]] <<- blood_PK_bo2_fit$hessian
    value_blood[[x]] <<- blood_PK_bo2_fit$value
    10^(blood_PK_bo2_fit$par)
  } else {
    blood_PK_bo2_fit <- optim(par=initial_blood_PK_parameters[c('L_1','L_2','lambda_1','lambda_2')],fn=objective_function_for_biexponential_blood_PK_model, lower=lower_bounds_blood, time=time_blood,DV=DV,method='L-BFGS-B',hessian=T,control=list(maxit=1000))
    H_blood[[x]] <<- blood_PK_bo2_fit$hessian
    blood_PK_bo2_fit$par
    value_blood[[x]] <<- blood_PK_bo2_fit$value
  }
})

## generate plots for the blood curve
prediction_times   <- seq(0,max(data_heart$tumor_time),length=1000)
day_list = sapply(1:num_of_day, function(x) {
  rep(unique(data_ana$'day')[x], length(prediction_times)) # generate the animal list which will be combined with the simulated dataframe
})
day_list = melt(day_list)['value']

xx <- seq(0,max(data_heart$tumor_time),length=length(prediction_times))
yy <- lapply(1:num_of_day, function(x) {
  aa[[x]]['L_1']*exp(-aa[[x]]['lambda_1']*xx) + aa[[x]]['L_2']*exp(-aa[[x]]['lambda_2']*xx)
})
yy = do.call(cbind,yy)  

heart_fit_data=data.frame(xx,yy)
heart_fit_data = melt(heart_fit_data, id='xx')
heart_fit_data['variable']=day_list
colnames(heart_fit_data)=c('xx','day','value')

## differential equation of the model that will be solved and fitted
model <- function(t,y,parms){

  # parms['epsilon']=0.25 # if epsilon is fixed
  hematocrit=0.45
  #epsilon = P_and_epsilon["epsilon"]
  Rl_blood = aa[[parms['ticker']]]['L_1']*(nM_converstion)*exp(-aa[[parms['ticker']]]['lambda_1']/60*t) + aa[[parms['ticker']]]['L_2']*(nM_converstion)*exp(-aa[[parms['ticker']]]['lambda_2']/60*t)
  Rl_p = Rl_blood/(1-hematocrit) 
   
  AAgKD <- (parms['Ag']/v_e + K_D)
  f_u_temp <<- # ran into problem when divided by y[1], so use f_u_temp instead
  
  Rl_i.st <<- y[1] # C_extravascular
  yd1 <-  parms['k_trans']*Rl_p - parms['k_trans']*f_u_temp - parms['k_endo']*Rl_i.st/parms['epsilon'] + parms['k_endo']*f_u_temp - parms['k_ns']*f_u_temp
  # Intracellular radiolabel, concentration in total tumor volume
  Rl_i.c <<- y[2]
  yd2 <- parms['k_endo']*Rl_i.st/parms['epsilon'] - parms['k_endo']*f_u_temp + parms['k_ns']*f_u_temp - parms['k_elim']*Rl_i.c/(1-v_e)
  
  # [Rl] in plasma, concentration in total tumor volume:
  RRl_p <- Rl_p*parms['v_p']    
  
  list(c(yd1,yd2),Rl_p, y[1] + y[2] + RRl_p)
}

## Fit the model
FitIt <- function(initial_parms,ticker,time,data_fit){
   
  # for lsoda 
  time_rep = unique(time)
  initial_conditions <- c(0,0)  # initial values for Rl_s
  
  if (log10_fit == 'yes'){
    parms <- 10^(initial_parms) 
  } else {
    parms <- initial_parms
  }
   
  # for fitting
  DV = data_fit*(nM_converstion)
  parms['ticker']=ticker
  # First time point for lsoda must be the time for which the initial conditions hold: 
  #RRl <<- lsoda(y=initial_conditions, times=c(0,time), func=model, parms=parms, rtol=1e-6, atol=1e-6)
  RRl <<- lsoda(y=initial_conditions, times=c(0,time_rep), func=model, parms=parms, rtol=1e-6, atol=1e-6)
  colnames(RRl) <<- c('mins','Rl_i.st','RL_i.c','Rl_plasma','Rl_total')
  # Drop prediction for zero time point because there is no DV at that time
  # (this time point was only added for  lsoda  as the time point for which 
  #  the initial conditions are given):
  
  temp_sig=RRl[is.element(RRl[,'mins'],time_rep),'Rl_total']
  temp_time = RRl[is.element(RRl[,'mins'],time_rep),'mins']
  
  m=match(time, temp_time)
  PRED      <- temp_sig[m]
  log_PRED  <- log(PRED)
  log_DV    <- log(DV)
  # Sum of Squares:
  
  #SSQ <- sum((log_DV - log_PRED)^2)
  SSQ <- sum((DV - PRED)^2)
  return(SSQ)
}

data_tumor=dcast(data_ana, tumor_time ~ day, value.var='tumor') # take the heart data from long to wide
data_tumor=data_tumor[ignore_dyna_points:length(data_tumor$'tumor_time'),] # remove the first 3 data point
data_tumor=melt(data_tumor, id='tumor_time') # melt the data to long form
H_tumor=list()
value_tumor=list()
convergence = list()

bb = lapply(1:num_of_day, function(x) {
  animal_day_tick = unique(data_ana$'day')[x]
  tumor_time = data_tumor[data_tumor['variable']==animal_day_tick,]$tumor_time
  tumor_data = data_tumor[data_tumor['variable']==animal_day_tick,]$value
  
  if (log10_fit == 'yes'){
    FitIt_result <- optim(par=log10(initial_parms), fn=FitIt, ticker=x, time=tumor_time, data_fit=tumor_data, method='L-BFGS-B', lower=log10(lower_bounds_tumor), hessian=T, control=list(maxit=1000, ndeps=c(10^(-6), 10^(-3), 10^(-3), 10^(-6), 10^(-6), 10^(-6), 10^(-3))))
    # FitIt_result <- optimx(par=log10(initial_parms), fn=FitIt, time=tumor_time, data_fit=tumor_data, method='Nelder-Mead',  control=list(maxit=1000,trace=4))
    H_tumor[[x]] <- FitIt_result$hessian
    H_tumor[[x]] <<- H_tumor[[x]]
    value_tumor[[x]] <<- FitIt_result$value
    convergence[[x]] <<- FitIt_result$convergence
    fitted_parms=10^(FitIt_result$par)
  } else {
    FitIt_result <- optim(par=initial_parms, fn=FitIt, time=tumor_time, data_fit=tumor_data, method='L-BFGS-B', lower=lower_bounds_tumor, hessian=T, control=list(maxit=10000, ndeps=c(10^(-6), 10^(-3), 10^(-3), 10^(-6), 10^(-6), 10^(-6), 10^(-3))))
    H_tumor[[x]] <- FitIt_result$hessian
    H_tumor[[x]] <<- H_tumor[[x]]
    value_tumor[[x]] <<- FitIt_result$value
    fitted_parms=FitIt_result$par
    # FitIt_result <- optimx(par=initial_parms, fn=FitIt, time=tumor_time, data_fit=tumor_data, method='Nelder-Mead',  control=list(maxit=1000,trace=4))
  }
})

#fitted_parms['epsilon']=0.25  # if epsilon is fixed
# generate plot for the tumor fit
initial_conditions <- c(0,0)

Rl_temp = lapply(1:num_of_day, function(x) {
  fitted_parms = bb[[x]]
  fitted_parms['ticker'] = x
  lsoda(y=initial_conditions, times=prediction_times,func=model,parms=fitted_parms, rtol=1e-6,atol=1e-6)
})
Rl <<- do.call(rbind,Rl_temp)
colnames(Rl)=c('time','Rl_i.st','RL_i.c','Rl_plasma','Rl_tot')
# making the plot for fitted tumor curve
tumor_fit_data=data.frame(day_list, Rl[,'time'],Rl[,'Rl_tot']/(nM_converstion))
colnames(tumor_fit_data)=c('day','time','Rl_tot')

# plot the fitted lines with a shared legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend = get_legend(p_heart)

p_heart_fit=p_heart + 
  geom_line(data=heart_fit_data, aes(x=xx, y=value), size=1, alpha=0.7) +
  theme(legend.position="none")
p_tumor_fit=p_tumor + 
  geom_line(data=tumor_fit_data, aes(x=time, y=Rl_tot), size = 1, alpha=0.7) +
  theme(legend.position="none")
grid.arrange(p_heart_fit, p_tumor_fit, legend, ncol=3, widths=c(2.3, 2.3, 0.8))

## Calculating Errors from Hessian
total_data_num_blood = dim(data_heart)[1]/num_of_animal
total_parms_num_blood = 4
total_data_num_tumor = dim(data_tumor)[1]/num_of_animal
total_parms_num_tumor = 3

std_lsqr_blood=list()
std_lsqr_tumor=list()
std_test=list()
std_lsqr = lapply(1:num_of_day, function(x) {
  if (log10_fit == 'yes'){
    std_lsqr_10_blood = 10^sqrt(diag(2*value_blood[[x]]/(total_data_num_blood - total_parms_num_blood) * solve(H_blood[[x]]))) # this is for least-square fit
    # std_lsqr_10_blood = 10^sqrt(diag(2*blood_PK_bo2_fit$value/(total_data_num_blood - total_parms_num_blood) * solve(H_blood[1:3,1:3]))) # this is for least-square fit
    std_lsqr_blood[[x]] <<- aa[[x]]*(std_lsqr_10_blood - 1)
    #std_lsqr_blood = aa[1:3]*(std_lsqr_10_blood - 1)
    std_lsqr_10_tumor = 10^sqrt(diag(2*value_tumor[[x]]/(total_data_num_tumor - total_parms_num_tumor) * solve(H_tumor[[x]]))) # this is for least-square fit
    #std_lsqr_tumor = fitted_parms[1:7]*(std_lsqr_10_tumor - 1)
    std_lsqr_tumor[[x]] <<- bb[[x]][1:3]*(std_lsqr_10_tumor - 1)
  } else {
    std_lsqr_blood <<- sqrt(diag(2*value_blood[[x]]/(total_data_num_blood - total_parms_num_blood) * solve(H_blood[[x]])))
    std_lsqr_tumor <<- sqrt(diag(2*value_tumor/(total_data_num_tumor - total_parms_num_tumor) * solve(H_tumor[[x]])))
  }
})

aa
std_lsqr_blood

bb
std_lsqr_tumor

convergence

####################################################
## Calculate the confidence interval the hard way ##
####################################################
# get the SSQ value for each day
# cc = lapply(1:num_of_day, function(x) {
#   animal_day_tick = unique(data_ana$'day')[x]
#   tumor_time = data_tumor[data_tumor['variable']==animal_day_tick,]$tumor_time
#   tumor_data = data_tumor[data_tumor['variable']==animal_day_tick,]$value
#   fitted_parms = bb[[x]]
#   
#   LOF_0 <<- FitIt(log10(fitted_parms), x, tumor_time, tumor_data)
#   
#   k_trans_eval=c(seq(3.1e-4,bb[[x]]['k_trans'],0.1e-4),seq(bb[[x]]['k_trans'],16e-4,0.1e-4))
#   v_p_eval=c(seq(0.028,bb[[x]]['v_p'],0.0005),seq(bb[[x]]['v_p'],0.032,0.0005))
#   epsilon_eval=c(seq(0.025,bb[[x]]['epsilon'],0.0005),seq(bb[[x]]['epsilon'],0.055,0.0005))
#   
#   grid = expand.grid(k_trans=k_trans_eval, v_p=v_p_eval, epsilon=epsilon_eval)
#   Grid_LOF<<-cbind(grid, result=apply(grid, 1, function(y) FitIt(log10(y), x, tumor_time, tumor_data)))
#   LOF<<-Grid_LOF['result']
#   
#   chisq <<- abs((LOF_0 - LOF)/LOF*(total_data_num_tumor-3)) # <= pchisq(0.95,1)
#   Grid_LOF[chisq<=pchisq(0.95,3),]
# })
# 
# cc
# feels like this confidence interval is much more narrower than the standard error predicted by hessian. 
# the values of standard error predicted by hessian look more reasonable. 

# P = P_and_epsilon[['P']] # 1/hr
# P = P/60 # 1/min
# 
# S_over_V_ini = initial_parms['k_trans']/P
# k_pip_fit = fitted_parms['k_pip']
# S_over_V = k_pip_fit/P
# names(S_over_V)='S/V'
# 
# initial_blood_PK_parameters[7:8]=initial_blood_PK_parameters[7:8]*24 # convert lambda unit from 1/hr to 1/day
# aa[3:4]=aa[3:4]*24 # convert lambda unit from 1/hr to 1/day
# std_lsqr_blood[3:4]=std_lsqr_blood[3:4]*24
# 
# # Result Output
# initial_tumor_parms = c(initial_blood_PK_parameters[c('L_1','L_2','lambda_1','lambda_2')],'S_over_V'=S_over_V_ini,initial_parms, 'value'=0)
# fit_tumor_parms = c(aa[1:4],'S_over_V'=S_over_V,fitted_parms[1:7],log10(fitted_parms['value']))
# 
# #rbind(initial_blood_parms,fit_blood_parms)
# final_table = rbind(initial_tumor_parms,fit_tumor_parms)
# 
# std_lsqr_blood
# std_lsqr_tumor
# 
# total_animals = dim(unique(data['animal']))[1]
# std_error_mean_blood = std_lsqr_blood/sqrt(total_animals)
# std_error_mean_tumor = std_lsqr_tumor/sqrt(total_animals)
# 
# final_table
# std_error_mean_blood
# std_error_mean_tumor
# 
# FitIt_result$convcode
