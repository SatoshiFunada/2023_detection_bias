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

#Subgroup=Tumor site
summary(as.factor(data_main_ORR$Tumor))

pdf("ORR_sub_type.pdf", width = 12, height = 11)
forest(res_ORR, atransf=exp, order=Tumor, 
       showweights=TRUE, textpos=c(-5,4),
       slab=paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), ylim=c(-1, 132), cex=.4,
       rows=c(3:16,21:23,26,31:34,37,42:46,51:52,55,60:63,68:85,90:94,99:100,103,108:116,121:122,125,128),
       mlab=mlabfun("Overall", res_ORR))
op <- par(cex=.4, font=2)
text(0, 131, "Ratio of odds ratios (95% CI)")
text(2.7, 131, "Weight")
### set font expansion factor (as in forest() above) and use a bold font

### add text for the subgroups
text(-5, c(129,126,123,117,104,101,95,86,64,56,53,47,38,35,27,24,17), pos=4, c("Squamous cell carcinoma of the lung", "Small-cell lung cancer", "Sarcoma", "Renal cell cancer", 
                                                                                   "Pleural mesohelioma","Pancretic cancer", "Ovarian cancer", "Non-small cell lung cancer",
                                                                                   "Melanoma", "Malignant mesothelioma", "Hepatocellular cancer",
                                                                                   "Head and Neck Cancer", "Gastrointestinal stromal tumours", 
                                                                                   "Gastric and gastroesophageal junction cancer", "Endometrial carcinoma", "Colorectal cancer", "Breast cancer"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.sarcoma <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Sarcoma"))
res.rcc <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Renal cell cancer"))
res.panc <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Pancretic cancer"))
res.ovarian <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Ovarian cancer"))
res.nsclc <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Non-small cell lung cancer"))
res.melanoma <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Melanoma"))
res.hepa <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Hepatocellular cancer"))
res.head <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Head and Neck Cancer"))
res.gast <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Gastric and gastroesophageal junction cancer"))
res.colo <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Colorectal cancer"))
res.breast <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", subset = (Tumor=="Breast cancer"))

### add summary polygons for the three subgroups
addpoly(res.sarcoma, row=119, mlab=mlabfun("Subtotal", res.sarcoma))
addpoly(res.rcc, row= 106, mlab=mlabfun("Subtotal", res.rcc))
addpoly(res.panc, row= 97, mlab=mlabfun("Subtotal", res.panc))
addpoly(res.ovarian, row=88, mlab=mlabfun("Subtotal", res.ovarian))
addpoly(res.nsclc, row= 66, mlab=mlabfun("Subtotal", res.nsclc))
addpoly(res.melanoma, row=58, mlab=mlabfun("Subtotal", res.hepa))
addpoly(res.hepa, row= 49, mlab=mlabfun("Subtotal", res.head))
addpoly(res.head, row= 40, mlab=mlabfun("Subtotal", res.head))
addpoly(res.gast, row=29, mlab=mlabfun("Subtotal", res.gast))
addpoly(res.colo, row=19, mlab=mlabfun("Subtotal", res.colo))
addpoly(res.breast, row= 1, mlab=mlabfun("Subtotal", res.breast))

### fit meta-regression model to test for subgroup differences
res_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main_ORR, method="REML", mods = ~ Tumor)

### add text for the test of subgroup differences
text(-5, -2, pos=4, cex=.4, bquote(paste("Test for Subgroup Differences: ",
                                             Q[M], " = ", .(formatC(res_ORR$QM, digits=2, format="f")), ", df = ", .(res_ORR$p - 1),
                                             ", p = ", .(formatC(res_ORR$QMp, digits=2, format="f")))))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

