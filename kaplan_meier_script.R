######################################################################################################
##
##
##                KAPLAN MEIER CURVES COMPARING SURVIVAL BETWEEN IMMUNE CLUSTERS
##
##
#####################################################################################################

library(survminer)
library(survival)
library(dplyr)


#read in deconvoluted data
immun <- as.data.frame(read.csv("/Users/jameskelly/Documents/IMMUNEESTIMATION/IMMUNE_PROPORTION.csv"))
#load file containing associated clinical info
load("/Users/jameskelly/Downloads/ucec.raw_counts.RData")


CLINICAL_INFO <- ucec_cdr_immune_clin
CLINICAL_INFO <- CLINICAL_INFO[,-1]

rownames(immun) <- immun[,1]
immun <- immun[,-1]
rownames(immun) == CLINICAL_INFO$sample ##TRUE
clinical_n_immun <- cbind(immun, CLINICAL_INFO)
diseas_spec_surv <- clinical_n_immun[,c(23,49,50,26)] 
diseas_spec_surv <-  dplyr::filter(diseas_spec_surv, P.value<.05)

immun <- dplyr::filter(immun, P.value<.05)
rownames(immun) == diseas_spec_surv$sample ## TRUE


immun_cells_only <- immun[,1:22] # remove everything except immune cell estimations
clstr <- kmeans(immun_cells_only, 2) # Create a column of associated clusters
clstr <- clstr$cluster

diseas_spec_surv$CLUSTER <- clstr

surv_object <- Surv(time = diseas_spec_surv$DSS.time, event = diseas_spec_surv$DSS)
CLUSTER <- diseas_spec_surv$CLUSTER
fit1_. <- survfit(surv_object ~ CLUSTER, data = diseas_spec_surv)
ggsurvplot(fit1_., data = diseas_spec_surv, pval = T,
           title = "Dis Spec Surv between clstrs",
           legend = "bottom",
           xlab = "Time(Days)",
           ylab = "Survival Probability",
           legend.title = "Cluster"
)
