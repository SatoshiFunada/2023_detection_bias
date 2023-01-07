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

#RHR for PFS
#yi=effect size, vi=variance, sei=standard error
res_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0, data=data_main, method="REML")
res0.25_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0.25, data=data_main, method="REML")
res0.5_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0.5, data=data_main, method="REML")
res0.75_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0.75, data=data_main, method="REML")
res0.95_PFS <- rma(yi=lnRHR, sei=lnRHR_SE_p_0.95, data=data_main, method="REML")

#ρ=0
pdf("PFS_0.pdf", width = 11, height = 11)
forest(res_PFS, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
text(0, 94, "Ratio of hazard ratios (95% CI)")
text(1.7, 94, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res_PFS$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res_PFS$QE, digits=2, format="f")), ", df = ", .(res_PFS$k - res_PFS$p),
                                         " (P = ", .(formatC(res_PFS$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res_PFS$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res_PFS$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res_PFS$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.25
pdf("PFS_0.25.pdf", width = 11, height = 11)
forest(res0.25_PFS, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
text(0, 94, "Ratio of hazard ratios (95% CI)")
text(1.7, 94, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.25_PFS$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.25_PFS$QE, digits=2, format="f")), ", df = ", .(res0.25_PFS$k - res0.25_PFS$p),
                                         " (P = ", .(formatC(res0.25_PFS$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.25_PFS$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.25_PFS$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.25_PFS$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.5
pdf("PFS_0.50.pdf", width = 11, height = 11)
forest(res0.5_PFS, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
text(0, 94, "Ratio of hazard ratios (95% CI)")
text(1.7, 94, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.5_PFS$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.5_PFS$QE, digits=2, format="f")), ", df = ", .(res0.5_PFS$k - res0.5_PFS$p),
                                         " (P = ", .(formatC(res0.5_PFS$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.5_PFS$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.5_PFS$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.5_PFS$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.75
pdf("PFS_0.75.pdf", width = 11, height = 11)
forest(res0.75_PFS, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
text(0, 94, "Ratio of hazard ratios (95% CI)")
text(1.7, 94, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.75_PFS$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.75_PFS$QE, digits=2, format="f")), ", df = ", .(res0.75_PFS$k - res0.75_PFS$p),
                                         " (P = ", .(formatC(res0.75_PFS$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.75_PFS$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.75_PFS$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.75_PFS$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.95
pdf("PFS_0.95.pdf", width = 11, height = 11)
forest(res0.95_PFS, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","RHR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
text(0, 94, "Ratio of hazard ratios (95% CI)")
text(1.7, 94, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.95_PFS$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.95_PFS$QE, digits=2, format="f")), ", df = ", .(res0.95_PFS$k - res0.95_PFS$p),
                                         " (P = ", .(formatC(res0.95_PFS$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.95_PFS$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.95_PFS$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.95_PFS$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Local exaggerate","Central exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ROR for ORR
#yi=effect size, vi=variance, sei=standard error
res_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0, data=data_main, method="REML")
res0.25_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0.25, data=data_main, method="REML")
res0.5_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0.5, data=data_main, method="REML")
res0.75_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0.75, data=data_main, method="REML")
res0.95_ORR <- rma(yi=lnROR, sei=lnROR_SE_p_0.95, data=data_main, method="REML")

#ρ=0
pdf("ORR_0.pdf", width = 11, height = 11)
forest(res_ORR, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
### add text with heterogeneity statistics
text(0, 76, "Ratio of odds ratios (95% CI)")
text(1.7, 76, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res_ORR$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res_ORR$QE, digits=2, format="f")), ", df = ", .(res_ORR$k - res_ORR$p),
                                         " (P = ", .(formatC(res_ORR$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res_ORR$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res_ORR$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res_ORR$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.25
pdf("ORR_0.25.pdf", width = 11, height = 11)
forest(res0.25_ORR, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
### add text with heterogeneity statistics
text(0, 76, "Ratio of odds ratios (95% CI)")
text(1.7, 76, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.25_ORR$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.25_ORR$QE, digits=2, format="f")), ", df = ", .(res0.25_ORR$k - res0.25_ORR$p),
                                         " (P = ", .(formatC(res0.25_ORR$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.25_ORR$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.25_ORR$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.25_ORR$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.5
pdf("ORR_0.50.pdf", width = 11, height = 11)
forest(res0.5_ORR, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
### add text with heterogeneity statistics
text(0, 76, "Ratio of odds ratios (95% CI)")
text(1.7, 76, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.5_ORR$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.5_ORR$QE, digits=2, format="f")), ", df = ", .(res0.5_ORR$k - res0.5_ORR$p),
                                         " (P = ", .(formatC(res0.5_ORR$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.5_ORR$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.5_ORR$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.5_ORR$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.75
pdf("ORR_0.75.pdf", width = 11, height = 11)
forest(res0.75_ORR, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
### add text with heterogeneity statistics
text(0, 76, "Ratio of odds ratios (95% CI)")
text(1.7, 76, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.75_ORR$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.75_ORR$QE, digits=2, format="f")), ", df = ", .(res0.75_ORR$k - res0.75_ORR$p),
                                         " (P = ", .(formatC(res0.75_ORR$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.75_ORR$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.75_ORR$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.75_ORR$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()

#ρ=0.95
pdf("ORR_0.95.pdf", width = 11, height = 11)
forest(res0.95_ORR, atransf=exp, order="obs", 
       showweights=TRUE,textpos=c(-4,3),
       slab = paste(authors, year, sep=", "),
       header=c("Author(s) and Year","ROR, 95% CI"),
       at=log(c(.25,.5,1,2,4)), xlim=c(-10,8), cex=.4, mlab="")

### add text with heterogeneity statistics
op <- par(cex=.4, font=2)
### add text with heterogeneity statistics
text(0, 76, "Ratio of odds ratios (95% CI)")
text(1.7, 76, "Weight")
op <- par(op)
text(-4, -1, pos=4, cex=.4, bquote(paste("Heterogeneity: ", "Tau"^2, " = ",
                                         .(formatC(res0.95_ORR$tau2, digits=2, format="f")), "; ", "Chi"^2, " = ",
                                         .(formatC(res0.95_ORR$QE, digits=2, format="f")), ", df = ", .(res0.95_ORR$k - res0.95_ORR$p),
                                         " (P = ", .(formatC(res0.95_ORR$QEp, digits=2, format="f")), "); ", I^2, " = ",
                                         .(formatC(res0.95_ORR$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(-4, -2, pos=4, cex=.4, bquote(paste("Test for overall effect: Z = ",
                                         .(formatC(res0.95_ORR$zval, digits=2, format="f")),
                                         " (P = ", .(formatC(res0.95_ORR$pval, digits=2, format="f")), ")")))
text(log(c(.25,4)), font=2, -4, cex=.4, c("Central exaggerate","Local exaggerate"), pos=c(4,2), offset=-0.5)
dev.off()