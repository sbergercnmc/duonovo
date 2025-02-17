### PPV assessment when varying parameters
setwd("C:/Users/lboukas/duoNovo_results/Parameter_sensitivity")
all_outputs <- list.files()

########### Father-proband duos
###########
all_outputs_PF <- all_outputs[grep("PF", all_outputs)]
all_outputs_PF_pass_QC <- all_outputs_PF[pass_QC]


load(file = all_outputs_PF_pass_QC[1]) #loading the first data frame in order to use colnames and dimensions
ppv_parameter_sensitivity <- matrix(NA, nrow = 40, ncol = dim(ppv_df_pf)[2])
colnames(ppv_parameter_sensitivity) <- colnames(ppv_df_pf)
rownames(ppv_parameter_sensitivity) <- all_outputs_PF_pass_QC

for (i in 1:length(all_outputs_PF_pass_QC)){
  load(file = all_outputs_PF_pass_QC[i])
  ppv_df <- ppv_df_pf
  ppv_parameter_sensitivity[i, ] <- apply(ppv_df, 2, function(xx) xx[1]/xx[2])
}
ppv_parameter_sensitivity_pf <- ppv_parameter_sensitivity

number_called_parameter_sensitivity_pf <- matrix(NA, nrow = 40, 
                                                 ncol = dim(ppv_parameter_sensitivity_pf)[2])
for (i in 1:length(all_outputs_PF_pass_QC)){
  load(file = all_outputs_PF_pass_QC[i])
  ppv_df <- ppv_df_pf
  number_called_parameter_sensitivity_pf[i, ] <- apply(ppv_df, 2, function(xx) xx[2])
}
rownames(number_called_parameter_sensitivity_pf) <- rownames(ppv_parameter_sensitivity_pf)
colnames(number_called_parameter_sensitivity_pf) <- colnames(ppv_parameter_sensitivity_pf)


dn_called_parameter_sensitivity_pf <- matrix(NA, nrow = 40, 
                                                 ncol = dim(ppv_parameter_sensitivity_pf)[2])
for (i in 1:length(all_outputs_PF_pass_QC)){
  load(file = all_outputs_PF_pass_QC[i])
  ppv_df <- ppv_df_pf
  dn_called_parameter_sensitivity_pf[i, ] <- apply(ppv_df, 2, function(xx) xx[1])
}
rownames(dn_called_parameter_sensitivity_pf) <- rownames(ppv_parameter_sensitivity_pf)
colnames(dn_called_parameter_sensitivity_pf) <- colnames(ppv_parameter_sensitivity_pf)

########### Mother-proband duos
###########
all_outputs_PM <- all_outputs[grep("PM", all_outputs)]
all_outputs_PM_pass_QC <- all_outputs_PM[pass_QC]


load(file = all_outputs_PM_pass_QC[1]) #loading the first data frame in order to use colnames and dimensions
ppv_parameter_sensitivity <- matrix(NA, nrow = 40, ncol = dim(ppv_df_pm)[2])
colnames(ppv_parameter_sensitivity) <- colnames(ppv_df_pm)
rownames(ppv_parameter_sensitivity) <- all_outputs_PF_pass_QC

for (i in 1:length(all_outputs_PM_pass_QC)){
  load(file = all_outputs_PM_pass_QC[i])
  ppv_df <- ppv_df_pm
  ppv_parameter_sensitivity[i, ] <- apply(ppv_df, 2, function(xx) xx[1]/xx[2])
}
ppv_parameter_sensitivity_pm <- ppv_parameter_sensitivity


number_called_parameter_sensitivity_pm <- matrix(NA, nrow = 40, 
                                                 ncol = dim(ppv_parameter_sensitivity_pm)[2])
for (i in 1:length(all_outputs_PM_pass_QC)){
  load(file = all_outputs_PM_pass_QC[i])
  ppv_df <- ppv_df_pm
  number_called_parameter_sensitivity_pm[i, ] <- apply(ppv_df, 2, function(xx) xx[2])
}
rownames(number_called_parameter_sensitivity_pm) <- rownames(ppv_parameter_sensitivity_pm)
colnames(number_called_parameter_sensitivity_pm) <- colnames(ppv_parameter_sensitivity_pm)

dn_called_parameter_sensitivity_pm <- matrix(NA, nrow = 40, 
                                             ncol = dim(ppv_parameter_sensitivity_pm)[2])
for (i in 1:length(all_outputs_PM_pass_QC)){
  load(file = all_outputs_PM_pass_QC[i])
  ppv_df <- ppv_df_pm
  dn_called_parameter_sensitivity_pm[i, ] <- apply(ppv_df, 2, function(xx) xx[1])
}
rownames(dn_called_parameter_sensitivity_pm) <- rownames(ppv_parameter_sensitivity_pm)
colnames(dn_called_parameter_sensitivity_pm) <- colnames(ppv_parameter_sensitivity_pm)


###
###
collective_ppv_pf_parameter_sensitivity <- colSums(dn_called_parameter_sensitivity_pf)/colSums(number_called_parameter_sensitivity_pf)
  
save(dn_called_parameter_sensitivity_pf, dn_called_parameter_sensitivity_pm, 
     number_called_parameter_sensitivity_pf, number_called_parameter_sensitivity_pm,
     file = "C:/Users/lboukas/duoNovo_results/parameter_sensitivity.rda")


