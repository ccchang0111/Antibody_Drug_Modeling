library(ggplot2) #library for plotting
library(ggthemes)
library(reshape)
library(deSolve)
require(gridExtra)

# Simulation of Antibody Uptake using Saturable vs Non-Saturable Model
# Dose Dependent
doses = c(15.25e-5)
mouse_weight = 20 # g
# mouse_volume_distribution = 3 # mL
KD = 2 # nM
MW = 15 # kDa
hematocrit = 0.45

# calculate Ktrans from estimated P_day
source("/xxx/Initial_parameters.R")
iniP = Initial_parameters(MW=MW) # can be commented out if you know the permeability of the molecule
P_day = iniP[["P"]] # Permeability unit: cm/day
S_over_V = 23 # 1/cm
k_trans = S_over_V*P_day/24/60 # 1/min
L_1 = 29 # iniP[['L_1']]
L_2 = 1 # iniP[['L_2']]
lbd1 = 54.3 # iniP[['lambda_1']], 54.3
lbd2 = 2.37 # iniP[['lambda_2']], 2.37

parms = c(L1=L_1, L2=L_2, lambda1=lbd1, lambda2=lbd2, epsilon = 0.165, v_e=0.165, v_p = 0.020, k_trans = k_trans, k_ns = 0.21, k_elim = 0, Ag = 25, K_D=KD, k_endo = 0) 

###################
## Binding Model ##
###################
tumor_residualizing_predictions_saturation <- function(t,y,parms=parms){
  # This model is for 'nonsaturable target'.
  p        <- parms
  # mouse_PK_parms = mouse_PK(dose)
  L1 = p[['L1']]
  L2 = p[['L2']]
  lambda1 = p[['lambda1']]
  lambda2 = p[['lambda2']]
  # Extract from  parms  the first row where  parms[,'group']==group :  
  
  Ag=p['Ag']
  v_p=p['v_p']
  k_ns=p['k_ns']
  k_elim=p['k_elim']
  k_endocytosis=p['k_endo']
  K_D=p['K_D']
  ve=p['v_e']
  
  Rl_p <- (L1*exp(-lambda1*t) +  L2*exp(-lambda2*t))/(1-hematocrit)
  
  # Interstitial radiolabel, concentration in total tumor volume
  f_u <<- K_D/(Ag/ve + K_D)
  
  Rl_i.st_unsat = y[1]
  yd1 = (p['k_trans']*(60*24)*Rl_p - p['k_trans']*(60*24)*f_u*Rl_i.st_unsat/p['epsilon'] - (k_endocytosis*(1 - f_u) + k_ns*f_u)*Rl_i.st_unsat/p['epsilon'])
  Rl_i.c_unsat = y[2]
  yd2 = (k_endocytosis*(1 - f_u) + k_ns*f_u)*Rl_i.st_unsat/p['epsilon'] - k_elim*Rl_i.c_unsat/(1-ve)
  
  # [Rl] in plasma, concentration in total tumor volume:
  RRl_p <- Rl_p*v_p   
  list(c(yd1,yd2),Rl_p*(1-hematocrit),y[1] + y[2] + RRl_p,f_u)
} # end of Function

##############################
## Generate Simulated Curve ##
##############################

initial_conditions = c(0,0)
t_hr = 24 # time course of uptake curve
t_day = t_hr/24
prediction_times = seq(0,t_day,length=1000)

# create animal list that will be combine with the data
dose_list = sapply(1:length(doses), function(x) {
  rep(doses[x], length(prediction_times)) # generate the animal list which will be combined with the simulated dataframe
})
dose_list = melt(dose_list)['value']

## Non-Binding control ##
datalist_sat_con = list()
for (i in 1:length(doses)){
  dose = doses[i]
  parms['K_D'] = 100000 # 100 uM, like a nonbinding control
  #mouse_PK_parms = mouse_PK(dose)
  #new_parms = c(parms,mouse_PK_parms)
  
  Rl_sat_con <- lsoda(y=initial_conditions,times=prediction_times,func=tumor_residualizing_predictions_saturation,parms=parms, rtol=1e-6, atol=1e-6)
  colnames(Rl_sat_con) <- c('days','Rl_i.st','Rl_i.c','Rl_blood','Rl_total_sat','f_u')
  
  # convert the nM unit to %ID/g
  # nM_pIDpg_conversion = 100*(1e-9)/((dose * mouse_weight*1e-3)/(MW*1000))
  time_hr = Rl_sat_con[,'days']*24
  datalist_sat_con[[i]] = data.frame(time_hr,Rl_sat_con[,'Rl_blood'], Rl_sat_con[,'Rl_total_sat'])
  colnames(datalist_sat_con[[i]]) = c('hrs','heart','tumor')
}
datalist_sat_con = do.call(rbind, datalist_sat_con)
datalist_sat_con['dose'] = dose_list

## Binding Peptide ##
datalist_sat = list()
for (i in 1:length(doses)){
  dose = doses[i]
  parms['K_D'] = KD # change the KD back to the binding peptide's KD
  #mouse_PK_parms = mouse_PK(dose)
  #new_parms = c(parms,mouse_PK_parms)
  
  Rl_sat <- lsoda(y=initial_conditions,times=prediction_times,func=tumor_residualizing_predictions_saturation,parms=parms, rtol=1e-6, atol=1e-6)
  colnames(Rl_sat) <- c('days','Rl_i.st','Rl_i.c','Rl_blood','Rl_total_sat','f_u')

  # convert the nM unit to %ID/g
  #nM_pIDpg_conversion = 100*(1e-9)/((dose * mouse_weight*1e-3)/(MW*1000))
  time_hr = Rl_sat[,'days']*24
  datalist_sat[[i]] = data.frame(time_hr,Rl_sat[,'Rl_blood'], Rl_sat[,'Rl_total_sat'])
  colnames(datalist_sat[[i]]) = c('hrs','heart','tumor')
}
datalist_sat = do.call(rbind, datalist_sat)
datalist_sat['dose'] = dose_list

# color for the plot
cols <- c("#666666", "#999999")

p_heart_all = ggplot(data = datalist_sat, aes(x=hrs, y=heart, colour = factor(dose))) +
  geom_line(color="#FF0000", size=1, alpha=0.7) +
  geom_line(data = datalist_sat_con, aes(x=hrs, y=heart, colour = factor(dose)), color = "#003399", size=1, alpha=0.7, linetype=2) +
  scale_colour_manual(values = cols) +
  # scale_x_log10() +
  labs(x = 'hrs', y = 'Heart (%ID/g)') +
  theme(legend.position="none")

p_tumor_all = ggplot(data = datalist_sat, aes(x=hrs, y=tumor, colour = factor(dose))) +
  geom_line(color="#FF0000", size=1, alpha=0.7) +
  geom_line(data = datalist_sat_con, aes(x=hrs, y=tumor, colour = factor(dose)), color = "#003399", size=1, alpha=0.7, linetype=2) +
  scale_colour_manual(values = cols) +
  labs(x = 'hrs', y = 'Tumor (%ID/g)') +
  theme(legend.position="none")

grid.arrange(p_heart_all, p_tumor_all, ncol=2) 