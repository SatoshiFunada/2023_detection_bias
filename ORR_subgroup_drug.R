library(tidyverse)
library(metafor)
data <- read.csv("data/Include_114.csv")

data_main <- data %>%
  select(ID, Blinding_Central, Trial_phase, Tumor, Funding, Trial_arm, Treatment_line,Treatment_line_re,
         Drug_classification_intervention, Primary_outcome, Response_criteria,
         Randomized_patient_number_intervention,
         Randomized_patient_number_control,
         RHR, lnRHR, lnRHR_SE_p_0, lnRHR_SE_p_0.25, lnRHR_SE_p_0.5,
         lnRHR_SE_p_0.75, lnRHR_SE_p_0.95,
         ROR, lnROR, lnROR_SE_p_0, lnROR_SE_p_0.25, lnROR_SE_p_0.5,
         lnROR_SE_p_0.75, lnROR_SE_p_0.95,
         authors, year, pubmed_id) %>%
  mutate(Tumor=as.factor(Tumor)) %>% 
  mutate(Drug_class=case_when(Drug_classification_intervention=="Immune checkpoint inhibitor"~"ICI",
                              Drug_classification_intervention=="Immune checkpoint inhibitor and chemotherapy"~"ICI",
                              Drug_classification_intervention=="Immune checkpoint inhibitor and target therapy"~"ICI",
                              Drug_classification_intervention=="Immune checkpoint inhibitor and target therapy and chemotherapy"~"ICI",
                              Drug_classification_intervention=="Target therapy"~"Target_therapy",
                              Drug_classification_intervention=="Target therapy and chemotherapy"~"Target_therapy",
                              Drug_classification_intervention=="Target therapy and Hormone therapy"~"Target_therapy",
                              Drug_classification_intervention=="Target therapy plus best supportive care"~"Target_therapy",
                              Drug_classification_intervention=="Immunotherapy"~"Immunotherapy",
                              Drug_classification_intervention=="Immunetherapy and Chemotherapy"~"Immunotherapy",
                              Drug_classification_intervention=="Chemotherapy"~"Chemotherapy"))%>% 
  mutate(Trial_phase_re=case_when(Trial_phase=="Phase 3"~"Phase_3",
                                  Trial_phase=="Phase 2B"~"Phase_2",
                                  Trial_phase=="Phase 2"~"Phase_2"))%>%
  mutate(RHR=as.numeric(RHR),
         lnRHR=as.numeric(lnRHR),
         lnRHR_SE_p_0=as.numeric(lnRHR_SE_p_0),
         lnRHR_SE_p_0.25=as.numeric(lnRHR_SE_p_0.25),
         lnRHR_SE_p_0.5=as.numeric(lnRHR_SE_p_0.5),
         lnRHR_SE_p_0.75=as.numeric(lnRHR_SE_p_0.75),
         lnRHR_SE_p_0.95=as.numeric(lnRHR_SE_p_0.95),
         ROR=as.numeric(ROR),
         lnROR=as.numeric(lnROR),
         lnROR_SE_p_0=as.numeric(lnROR_SE_p_0),
         lnROR_SE_p_0.25=as.numeric(lnROR_SE_p_0.25),
         lnROR_SE_p_0.5=as.numeric(lnROR_SE_p_0.5),
         lnROR_SE_p_0.75=as.numeric(lnROR_SE_p_0.75),
         lnROR_SE_p_0.95=as.numeric(lnROR_SE_p_0.95))

data_main_ORR <- data_main%>%
  drop_na(ROR)


#yi=effect size, vi=variance, sei=standard error
#ROR
res_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML")

### a little helper function to add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, res_ORR) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res_ORR$QE, digits=2, format="f")),
                    ", df = ", .(res_ORR$k - res_ORR$p),
                    ", p ", .(metafor:::.pval(res_ORR$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    I^2, " = ", .(formatC(res_ORR$I2, digits=1, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res_ORR$tau2, digits=2, format="f")), ")")))}

#Subgroup=Drug classification
summary(as.factor(data_main_ORR$Drug_class))

pdf("ORR_sub_drug.pdf", width = 12, height = 11)
forest(res_ORR, atransf=exp, order=Drug_class, 
       showweights=TRUE,textpos=c(-5,4),
       slab=paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), ylim=c(-1, 92), cex=.4,
       rows=c(3:22,27:36,41:43,48:88),
       mlab=mlabfun("Overall", res_ORR))
op <- par(cex=0.4, font=2)
text(0, 91, "Ratio of odds ratios (95% CI)")
text(2.7, 91, "Weight")
### add text for the subgroups
text(-5, c(89,44,37,23), pos=4, c("Target therapy", 
                                   "Immunotherapy",
                                 "Immune checkpoint inhibitor",
                                 "Chemotherapy"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.t <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Drug_class=="Target_therapy"))
res.i <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Drug_class=="Immunotherapy"))
res.ici <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Drug_class=="ICI"))
res.c <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Drug_class=="Chemotherapy"))

### add summary polygons for the three subgroups
addpoly(res.t, row=46, mlab=mlabfun("Subtotal", res.t))
addpoly(res.i, row= 39, mlab=mlabfun("Subtotal", res.i))
addpoly(res.ici, row= 25, mlab=mlabfun("Subtotal", res.ici))
addpoly(res.c, row= 1, mlab=mlabfun("Subtotal", res.c))

### fit meta-regression model to test for subgroup differences
res_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main, method="REML", mods = ~ Drug_class)

### add text for the test of subgroup differences
text(-5, -2, pos=4, cex=.4, bquote(paste("Test for Subgroup Differences: ",
                                             Q[M], " = ", .(formatC(res_ORR$QM, digits=2, format="f")), ", df = ", .(res_ORR$p - 1),
                                             ", p = ", .(formatC(res_ORR$QMp, digits=2, format="f")))))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

