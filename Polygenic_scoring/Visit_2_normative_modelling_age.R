#Maria Di Biase, Feb 6th 2023

#This code will optimise a gamlss model and compute deviations based on the optimised model
#Two models are optimised: the first incorporates batch and the second incorporates batch and bmi
#The code can be split into 4 steps:
#1) GAMLSS optimisation, which uses all the samples
#2) Running the optimised model in half of the sample
#3) Predicting parameters and quantiles for hypothetical data (e.g., for batch 1)
#4) Computing deviations and bands in the other half of the sample (computed based on an individual's specific batch number)

#Inputs:
#Dataframe with age, phenotype (retinol), bmi and batch 
#Can amend how the half sample is defined at line 50
#Please search for $AMEND$ to find code requiring amendments to run on real data

#Outputs:
#1) fit parameters for each model evaluated (BIC, AIC and deviance)
#2) Batch-specific quantiles (can use an average for manuscript plotting purposes) and KS statistic 
#3) individual deviations. 
#4) individual bands labelled, band_infra, band_supra and band_normal, respectively

rm(list = ls())
setwd("~/Desktop/Retinol_GWAS/TwinsUK_app/E1205_20.01.2023/Maria_normative_modelling/") #$AMEND$
source('helpers.R')
library(ggpointdensity)
library(cowplot)
library(patchwork)
library(qqplotr)
library(R.matlab)
library(Hmisc)
library(rstatix)
library(ggpubr)
library(scoring)
library(gamlss)
library(Hmisc)
library(pracma)
library(ComplexHeatmap)
library(haven)

set.seed(5353672)                                                                                                                                                          
#set out path $AMEND$
PATH_OUT <- paste0(getwd(), "/gamlss_output_visit_2/")

# read data currently modeling baseline 
df <- read.table("../Long_format_full_retinol_TwinsUK.txt", sep = "\t", header = TRUE)

## Remove any missing values - retinol visit 2 as a test
df <- df %>% filter(!is.na(Retinol_visit2)) %>% filter(!is.na(BMI_visit2))

## Remove males

df <- df %>% filter(sex == "F")

mydata_tmp <- df %>% select(Age_visit2, Batch_visit2, BMI_visit2, Retinol_visit2, Last_digit, PublicID)
names(mydata_tmp)<-c("age","batch","bmi","retinol", "Last_digit", "ID")

mydata_tmp$batch <- as.factor(mydata_tmp$batch)
levels(mydata_tmp$batch) <- c(1, 2, 3, 4, 5)
mydata_tmp$batch <- as.integer(mydata_tmp$batch)

batch_length <- length(unique(mydata_tmp$batch))

#select half of the subjects for model training (STEP 2) - split based on twin pairs (Subset 2 as training to match genetics)
idx1 <- which(mydata_tmp$Last_digit == 2, arr.ind = T)
idx2 <- which(mydata_tmp$Last_digit == 1, arr.ind = T)


########################################################################################
print(paste("STEP 1 - Optimise Model"))
########################################################################################

#evaluate the following family distributions. Add more if you wish 
fams<-c("NO","NOF","LNO","LOGNO","GG","WEI","WEI2","WEI3","GB1","GB2","GA","BCT","GAF",
        "ST1","ST2","ST3","SEP","SEP1","SEP2","SEP3","SEP4","JSU","exGAUS","SHASH","GT","TF","IG")


#initialise matrices to store model 1 fit-related outputs
BIC_1<-matrix(NaN, nrow = length(fams), ncol = 1) #Bayesian information criterion BIC)
gvd_1<-matrix(NaN, nrow = length(fams), ncol = 1) #deviance
AIC_1<-matrix(NaN, nrow = length(fams), ncol = 1) #Akaike information criterion (AIC)
iterations_1<-matrix(NaN, nrow = length(fams), ncol = 1) #number of iterations for model to converge
converged_bin_1<-matrix(NaN, nrow = length(fams), ncol = 1) #did model converge (y/n)
models_list_1 <- list()

#initialise matrices to store model 2 fit-related outputs
BIC_2<-matrix(NaN, nrow = length(fams), ncol = 1) #Bayesian information criterion BIC)
gvd_2<-matrix(NaN, nrow = length(fams), ncol = 1) #deviance
AIC_2<-matrix(NaN, nrow = length(fams), ncol = 1) #Akaike information criterion (AIC)
iterations_2<-matrix(NaN, nrow = length(fams), ncol = 1) #number of iterations for model to converge
converged_bin_2<-matrix(NaN, nrow = length(fams), ncol = 1) #did model converge (y/n)
models_list_2 <- list()

for (mm in 1:length(fams)){ #loop over family distributions
  try({ #we need to add this because it is unlikely that all models will converge 
    
    #MODEL TRAINING
    control <- gamlss.control(n.cyc = 100, pars.start = list(mu = c(-10, 10), sigma = c(0.1, 10),
                                                             maxit = 100, nu.fo = c(1,5), nu.tau = c(0.1,1)))      
    
    #MODEL 1 - model batch only
    mdl1 <- gamlss(retinol ~ fp(age,npoly=1) + batch, # formula for mu
                   sigma.fo=~fp(age,npoly=1),
                   family=fams[mm],
                   data=mydata_tmp, #optimise on all data 
                   control = control,
                   robust=T)
    
    #MODEL 1 OUTPUT
    models_list_1[[mm]] <- mdl1
    names(models_list_1)[mm] <- paste0(fams[mm],mm)
    
    if (mdl1$converged) {
      print("Model 1 Converged : YES")
      converged_bin_1[mm]<-1
    } else {
      print("Model 1 Cnverged: NO")
      converged_bin_1[mm]<-0
    }
    
    
    #fit indices - model 1
    gvd_1[mm]<-deviance(mdl1,what="G") # validation global deviance  
    BIC_1[mm]<-mdl1$sbc
    AIC_1[mm]<-mdl1$aic
    iterations_1[mm]<-mdl1$iter 
    
    
    #MODEL 2 - model batch and BMI
    mdl2 <- gamlss(retinol ~ fp(age,npoly=1) + batch + bmi, # formula for mu
                   sigma.fo=~fp(age,npoly=1),
                   family=fams[mm],
                   data=mydata_tmp,  #optimise on all data 
                   control = control,
                   robust=T)
    
    models_list_2[[mm]] <- mdl2
    names(models_list_2)[mm] <- paste0(fams[mm],mm)
    
    if (mdl2$converged) {
      print("Model 2 Converged : YES")
      converged_bin_2[mm]<-1
    } else {
      print("Model 2 Cnverged: NO")
      converged_bin_2[mm]<-0
    }
    
    #fit indices - model 2
    gvd_2[mm]<-deviance(mdl2,what="G") # validation global deviance  
    BIC_2[mm]<-mdl2$sbc
    AIC_2[mm]<-mdl2$aic
    iterations_2[mm]<-mdl2$iter 
    
  }
  , silent = T) #END try
  
} #END loop over mm


#find best fitted family distribution 
best_fit_bic_idx_1 <- which.min(BIC_1)
best_fit_aic_idx_1 <- which.min(AIC_1)
deltaAIC_1 <- AIC_1[best_fit_bic_idx_1] - AIC_1[best_fit_aic_idx_1]

best_fit_bic_idx_2 <- which.min(BIC_2)
best_fit_aic_idx_2 <- which.min(AIC_2)
deltaAIC_2 <- AIC_2[best_fit_bic_idx_2] - AIC_2[best_fit_aic_idx_2]

#Use family distribution with lowest BIC (fit based on BIC and AIC typically agree but worth checking) 
fam_dist1=fams[best_fit_bic_idx_1] 
print(paste("Best family distribution for Model 1: ",fam_dist1,sep=""))

## [1] "Best family distribution for Model 1: exGAUS"
if(best_fit_bic_idx_1==best_fit_aic_idx_1 | abs(deltaAIC_1)<2) {
  print("MODEL 1 - AIC and BIC agree: YES")
} else {
  print("MODEL 1 - AIC and BIC agree: NO")
}
##"MODEL 1 - AIC and BIC agree: NO"
fam_dist2=fams[best_fit_bic_idx_2] 
print(paste("Best family distribution for Model 2: ",fam_dist2,sep=""))

## [1] "Best family distribution for Model 2: exGAUS

if(best_fit_bic_idx_2==best_fit_aic_idx_2 | abs(deltaAIC_2)<2) {
  print("MODEL 2 - AIC and BIC agree: YES")
} else {
  print("MODEL 2 - AIC and BIC agree: NO")
}
## "MODEL 2 - AIC and BIC agree: YES"
## Test BIC between with and without BMI - very similar < 2 between them

## Wit BMI fits marginally better delta AIC = 0.03

AIC_2[best_fit_bic_idx_2] - AIC_1[best_fit_bic_idx_1]

outname<-paste(PATH_OUT,"GAMLSSout_FIT_eval_famDist.txt",sep="")

## Outcome df 

Out_df <- as.data.frame(cbind(AIC_1, gvd_1, BIC_1, AIC_2, gvd_2, BIC_2, fams))
names(Out_df)<-c("AIC_no_BMI","gvd_no_BMI","BIC_no_BMI","AIC_w_BMI","gvd_w_BMI","BIC_w_BMI", "GAMLSS_family")

## Convert to numeric
outname<-paste(PATH_OUT,"GAMLSSout_FIT_eval_famDist.mat",sep="")
writeMat(outname,  AIC1=AIC_1, gvd1=gvd_1,  BIC1=BIC_1, 
         AIC2=AIC_2, gvd2=gvd_2,  BIC2=BIC_2, fams=fams)

## Matrix of fit metrics

Input_mat <- as.matrix(Out_df)

rownames(Input_mat) <- Input_mat[, 7]

Input_mat <- Input_mat[,c(1,3,4,6)]

Input_mat <- `dimnames<-`(`dim<-`(as.numeric(Input_mat), dim(Input_mat)), dimnames(Input_mat))

colnames(Input_mat) <- c("AIC (no BMI)", "BIC (no BMI)", "AIC (with BMI)", 
                         "BIC (with BMI)")

write.table(Input_mat, file = outname, sep="\t", quote = F)

Heatmap(na.omit(Input_mat), rect_gp = gpar(col = "white", lwd = 1),
        column_names_rot = 50, column_title_gp = grid::gpar(fontsize = 11.5, fontface="bold"),
        name = "Model fit value", column_title = "GAMLSS model fit", show_row_dend = F, show_column_dend = F,
        row_title = "GAMLSS family distribution")

## BCT seems to perform the best across both conditions considering both BIC and AIC in general across all visits

## Rename fam dist 1 var

fam_dist1 = "BCT"
fam_dist2 = "BCT"

########################################################################################
print(paste("STEP 2 - model half sample using optimised family distribution"))
########################################################################################

# initialise matrix for subject z-scores (deviations)
z_scores_unseen_subs<-matrix(NaN, nrow = length(mydata_tmp$age), ncol = 1) 
band_infra<-matrix(NaN, nrow = length(mydata_tmp$age), ncol = 1) 
band_supra<-matrix(NaN, nrow = length(mydata_tmp$age), ncol = 1) 
band_normal<-matrix(NaN, nrow = length(mydata_tmp$age), ncol = 1) 

## Just do non-BMI for now

for (idx_mdl in 1:2) { #loop over models 
  
  print(paste("running STEP 2 for model ",idx_mdl,sep=""))
  
  for (idx_sample in 1:2) { #loop over 2 half samples: model 1 half sample and predict deviations in the other
    
    if (idx_sample==1) {
      mydata <- mydata_tmp[idx1,] #index one half sample for modeling
      predict_sample <- idx2 #index second half sample for prediction
      mydata_predict <- mydata_tmp[idx2,]
      print(paste("Modeling first half sample..."))
    } else if (idx_sample==2) {
      mydata <- mydata_tmp[idx2,] #index one half sample for modeling
      predict_sample <- idx1 #index second half sample for prediction
      mydata_predict <- mydata_tmp[idx1,]
      print(paste("Modeling second half sample..."))
    }
    
    
    #Re-run best fit model on half sample
    if (idx_mdl==1) {
      fam_dist=fam_dist1
      mdl <- gamlss(retinol ~ fp(age,npoly=1) + batch, # formula for mu
                    sigma.fo=~fp(age,npoly=1),
                    family=fam_dist1,
                    data=mydata, #only use half sample
                    control = control,
                    robust=T)
    } else if (idx_mdl==2) {
      fam_dist=fam_dist2
      mdl <- gamlss(retinol ~ fp(age,npoly=1) + batch + bmi, # formula for mu
                    sigma.fo=~fp(age,npoly=1),
                    family=fam_dist2,
                    data=mydata, #only use half sample
                    control = control,
                    robust=T)
    }
    
    
    #Predict parameters on half sample
    if(fam_dist=="GB1" | fam_dist=="BCT" | fam_dist=="ST1" | fam_dist=="ST2" | fam_dist=="ST3" |
       fam_dist=="SEP" | fam_dist=="SEP1" | fam_dist=="SEP2" | fam_dist=="SEP3" | fam_dist=="SEP4"|
       fam_dist=="GB2" | fam_dist=="JSU" | fam_dist=="SHASH" |
       fam_dist=="GT" | fam_dist=="IG") {
      
      params<-predictAll(mdl,data=mydata, newdata=mydata,
                         output='matrix',type="response",
                         y.value="median",what=c("mu", "sigma", "nu", "tau"))
      
    } else if(fam_dist=="NOF" | fam_dist=="LNO" | fam_dist=="GG" | fam_dist=="exGAUS" | fam_dist=="TF"  ) {
      
      params<-predictAll(mdl,data=mydata, newdata=mydata,
                         output='matrix',type="response",
                         y.value="median",what=c("mu", "sigma", "nu"))
      
    } else if(fam_dist=="NO" | fam_dist=="NO2" | fam_dist=="LOGNO" | fam_dist=="WEI" |
              fam_dist=="WEI2"| fam_dist=="WEI3" | fam_dist=="GA"| fam_dist=="GAF") {
      params<-predictAll(mdl,data=mydata, newdata=mydata,
                         output='matrix',type="response",
                         y.value="median",what=c("mu", "sigma"))
      
    }
    
    
    #generate quantiles using best fit
    if(fam_dist=="NO") {
      quantiles <- pNO(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="NO2") {
      quantiles <- pNO2(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="NOF") {
      quantiles <- pNOF(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4]);
    } else if(fam_dist=="LNO") {
      quantiles <- pLNO(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4]);
    } else if(fam_dist=="LOGNO") {
      quantiles <- pLOGNO(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="GG") {
      quantiles <- pGG(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4]);
    } else if(fam_dist=="WEI") {
      quantiles <- pWEI(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="WEI2") {
      quantiles <- pWEI2(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="WEI3") {
      quantiles <- pWEI3(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="GB1") {
      quantiles <- pGB1(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="GB2") {
      quantiles <- pGB2(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="GA") {
      quantiles <- pGA(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="GAF") {
      quantiles <- pGAF(params[,1],mu=params[,2],sigma=params[,3]);
    } else if(fam_dist=="BCT") {
      quantiles <- pBCT(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="ST1") {
      quantiles <- pST1(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="ST2") {
      quantiles <- pST2(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="ST3") {
      quantiles <- pST3(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="ST4") {
      quantiles <- pST4(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="SEP") {
      quantiles <- pSEP(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="SEP1") {
      quantiles <- pSEP1(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="SEP2") {
      quantiles <- pSEP2(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="SEP3") {
      quantiles <- pSEP3(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="SEP4") {
      quantiles <- pSEP4(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="JSU") {
      quantiles <- pJSU(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="exGAUS") {
      quantiles <- pexGAUS(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4]);
    } else if(fam_dist=="SHASH") {
      quantiles <- pSHASH(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="GT") {
      quantiles <- pGT(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    } else if(fam_dist=="TF") {
      quantiles <- pTF(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4]);
    } else if(fam_dist=="IG") {
      quantiles <- pIG(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]);
    }
    
    
    #can check the following outputs to ensure model is well fit (although KS is conservative)
    z_randomized<-qnorm(quantiles, mean = 0, sd =1);
    ks_out<-ks.test(z_randomized,"pnorm") #Kolmogorov–Smirnov test to confirm appropriate fit
    if(ks_out$p.value<.05) {
      print(paste("model not well fit according to KS test, p=",ks_out$p.value,sep=""))
    }
    
    ########################################################################################
    print(paste("STEP 3 - predict hypothetical data"))
    ########################################################################################
    
    par(mfrow = c(2, 3))
    
    for (idx_batch in 1:batch_length) { #generate quantiles for each batch
      
      print(paste("run for batch",idx_batch,sep=" "))
      min_age<-min(mydata_tmp$age)
      max_age<-max(mydata_tmp$age)
      age_test<-linspace(min_age,max_age,1000)
      
      plot(retinol~age, data=mydata,  col="lightgray", pch=10, ylim=c(0,3.2),
           main=paste("Batch", idx_batch, sep=" "), xlab="Age", ylab="Retinol")
      predictions_quantiles<-matrix(data=0,ncol=5,nrow=1000)
      centile_vec=c(.05, 0.25, 0.5, 0.75, 0.95)
      for (i in 1:length(centile_vec)){
        Qua <- getQuantile(mdl, quantile=centile_vec[i],term="age",fixed.at=list(batch=idx_batch))
        out<-curve(Qua, min_age, max_age,  lwd=2, lty=1, add=T,col="red",n = 1000)
        predictions_quantiles[,i]=as.vector(out$y)
      }
      
      
      ########################################################################################
      print(paste("STEP 4 - compute deviations on independent sample half"))
      ########################################################################################
      
      # Interpolate the 50th percentile curve at the ages of the new subjects
      ind_subs<-which(mydata_predict$batch==idx_batch)
      mydata_new <- mydata_predict[ind_subs,] #index only subjects in batch idx_batch
      names(mydata_new)<-c("age","batch","bmi","retinol");
      
      # Compute the standard deviation
      sd_phenotype <- sd(mydata_new$retinol, na.rm = TRUE)
      
      for (idx_sub in 1:length(ind_subs)) {
        sub_age <- mydata_new$age[idx_sub] # ages of independent subjects
        sub_phenotype <- mydata_new$retinol[idx_sub]
        interpolated_percentile <- approx(age_test, predictions_quantiles[,3], xout = sub_age)$y 
        interpolated_percentile95 <- approx(age_test, predictions_quantiles[,5], xout = sub_age)$y 
        interpolated_percentile05 <- approx(age_test, predictions_quantiles[,1], xout = sub_age)$y 
        
        # Subtract the interpolated 50th percentile values from the responses of the new subjects
        difference <- sub_phenotype - interpolated_percentile
        
        # Divide the differences by the standard deviation of the response and store
        z_scores_unseen_subs[predict_sample[ind_subs[idx_sub]]] <- difference / sd_phenotype 
        
        # Determine individual bands (infra/supra normal)
        band_infra[predict_sample[ind_subs[idx_sub]]]<-ifelse(sub_phenotype<interpolated_percentile05, 1, 0)
        band_supra[predict_sample[ind_subs[idx_sub]]]<-ifelse(sub_phenotype>interpolated_percentile95, 1, 0)
        band_normal[predict_sample[ind_subs[idx_sub]]]<-ifelse(sub_phenotype >= interpolated_percentile05 & sub_phenotype <= interpolated_percentile95, 1, 0)
        
      } #END loop over idx_sub
      
      # index and plot infranormal subjects. Note that bands are computed  is based on quantiles predicted for a specific batch, whereas 
      # the below is plotting subjects that are considered infra/supra normal across all batches. This is why there appears to be 
      # miss-categorized data points 
      ind_subs_infra <- which(band_infra==1) 
      points(mydata_tmp$age[ind_subs_infra], mydata_tmp$retinol[ind_subs_infra], col = "red", pch = 19)
      
      #index and plot supranormal subjects 
      ind_subs_supra <- which(band_supra==1) 
      points(mydata_tmp$age[ind_subs_supra], mydata_tmp$retinol[ind_subs_supra], col = "blue", pch = 19)
      
      #spit out quantiles in preferred format (MATLAB)
      outname<-paste(PATH_OUT,"GAMLSSout_QUANTILES_BATCH",idx_batch,"_MODEL",idx_mdl,".mat",sep="")
      writeMat(outname, predictions_quantiles=as.matrix(predictions_quantiles), age=age_test,zscores=z_randomized,ks_out=ks_out)
      
    } #END loop over idx_batches
    
    
  } #end loop over idx_sample
  
  #spit out deviations in preferred format (MATLAB)
  outname<-paste(PATH_OUT,"GAMLSSout_INDV_DEVIATIONS","_MODEL",idx_mdl,".mat",sep="")
  writeMat(outname, z_scores_unseen_subs=z_scores_unseen_subs, band_infra=band_infra, band_supra=band_supra,band_normal=band_normal)
  
}  #END loop over idx_mdl


print(paste("SCRIPT complete"))

## Test identifying deviations

TE <- readMat("gamlss_output_visit_2/GAMLSSout_INDV_DEVIATIONS_MODEL1.mat")

Ids <- mydata_tmp %>% select(ID)

Dev_1 <- as.data.frame(cbind(Ids, band_normal, band_infra, band_supra))

Dev_1 <- rename(Dev_1, "publicid"="ID")

write.table(Dev_1, file="Visit_2_deviations_SHASH.txt",
            sep = "\t", row.names = F)
