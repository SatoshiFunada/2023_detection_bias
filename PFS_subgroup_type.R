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

#Subgroup=Tumor site
summary(as.factor(data_main_PFS$Tumor))

pdf("PFS_sub_type.pdf", width = 12, height = 11)
forest(res_PFS, atransf=exp, order=Tumor, 
       showweights=TRUE, textpos=c(-5,4),
       slab=paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), ylim=c(-1, 148), cex=.4,
       rows=c(3:22,25,28,33:35,38,43:45,48,51,54,59:66,69,74:96,101:108,113:114,117,122:134,139:141,144),
       mlab=mlabfun("Oveall", res_PFS))
op <- par(cex=.4, font=2)
text(0, 147, "Ratio of hazard ratios (95% CI)")
text(2.7, 147, "Weight")
### set font expansion factor (as in forest() above) and use a bold font
### add text for the subgroups
text(-5, c(145,142,135,118,115,109,97,70,67,55,52,49,46,39,36,29,26,23), pos=4, c("Squamous cell carcinoma of the lung", "Sarcoma", "Renal cell cancer", "Prostate cancer",
                               "Pancretic cancer", "Ovarian cancer", "Non-small cell lung cancer", "Neuroendocrine tumors",
                               "Melanoma", "Malignant mesothelioma", "Liposarcoma or Leiomyosarcoma", "Hepatocellular cancer",
                               "Head and Neck Cancer", "Glioblastoma", "Gastrointestinal stromal tumours", 
                               "Gastric and gastroesophageal junction cancer", "Endometrial carcinoma", "Breast cancer"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.sarcoma <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Sarcoma"))
res.rcc <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Renal cell cancer"))
res.panc <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Pancretic cancer"))
res.ovarian <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Ovarian cancer"))
res.nsclc <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Non-small cell lung cancer"))
res.melanoma <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Melanoma"))
res.head <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Head and Neck Cancer"))
res.gist <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Gastrointestinal stromal tumours"))
res.breast <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", subset = (Tumor=="Breast cancer"))

### add summary polygons for the three subgroups
addpoly(res.sarcoma, row=137, mlab=mlabfun("Subtotal", res.sarcoma))
addpoly(res.rcc, row= 120, mlab=mlabfun("Subtotal", res.rcc))
addpoly(res.panc, row= 111, mlab=mlabfun("Subtotal", res.panc))
addpoly(res.ovarian, row=99, mlab=mlabfun("Subtotal", res.ovarian))
addpoly(res.nsclc, row= 72, mlab=mlabfun("Subtotal", res.nsclc))
addpoly(res.melanoma, row=57, mlab=mlabfun("Subtotal", res.melanoma))
addpoly(res.head, row= 41, mlab=mlabfun("Subtotal", res.head))
addpoly(res.gist, row=31, mlab=mlabfun("Subtotal", res.gist))
addpoly(res.breast, row= 1, mlab=mlabfun("Subtotal", res.breast))

### fit meta-regression model to test for subgroup differences
res_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main_PFS, method="REML", mods = ~ Tumor)

### add text for the test of subgroup differences
text(-5, -2, pos=4, cex=.4, bquote(paste("Test for Subgroup Differences: ",
                                             Q[M], " = ", .(formatC(res_PFS$QM, digits=2, format="f")), ", df = ", .(res_PFS$p - 1),
                                             ", p = ", .(formatC(res_PFS$QMp, digits=2, format="f")))))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

