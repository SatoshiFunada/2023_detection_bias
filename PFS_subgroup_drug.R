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

data_main_PFS <- data_main%>%
  drop_na(RHR)


#yi=effect size, vi=variance, sei=standard error
#RHR
res_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML")

### a little helper function to add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, res_PFS) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res_PFS$QE, digits=2, format="f")),
                    ", df = ", .(res_PFS$k - res_PFS$p),
                    ", p ", .(metafor:::.pval(res_PFS$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    I^2, " = ", .(formatC(res_PFS$I2, digits=1, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res_PFS$tau2, digits=2, format="f")), ")")))}

#Subgroup=Drug classification
summary(as.factor(data_main_PFS$Drug_class))

pdf("PDF_sub_drug.pdf", width = 12, height = 11)
forest(res_PFS, atransf=exp, order=Drug_class, 
       showweights=TRUE,textpos=c(-5,4),
       slab=paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), ylim=c(-1, 106), cex=.4,
       rows=c(3:20,25:36,41:102),
       mlab=mlabfun("Overall", res_PFS))
op <- par(cex=.4, font=2)
text(0, 105, "Ratio of hazard ratios (95% CI)")
text(2.7, 105, "Weight")
### add text for the subgroups
text(-5, c(103,37,21), pos=4, c("Target therapy",
                               "Immune checkpoint inhibitor",
                               "Chemotherapy"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.t <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Drug_class=="Target_therapy"))
res.ici <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Drug_class=="ICI"))
res.c <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Drug_class=="Chemotherapy"))

### add summary polygons for the three subgroups
addpoly(res.t, row=39, mlab=mlabfun("Subtotal", res.t))
addpoly(res.ici, row= 23, mlab=mlabfun("Subtotal", res.ici))
addpoly(res.c, row= 1, mlab=mlabfun("Subtotal", res.c))

### fit meta-regression model to test for subgroup differences
res_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main, method="REML", mods = ~ Drug_class)

### add text for the test of subgroup differences
text(-5, -2, pos=4, cex=.4, bquote(paste("Test for Subgroup Differences: ",
                                              Q[M], " = ", .(formatC(res_PFS$QM, digits=2, format="f")), ", df = ", .(res_PFS$p - 1),
                                              ", p = ", .(formatC(res_PFS$QMp, digits=2, format="f")))))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()
