##################### Add conceptual model

addConceptual.Model <- function(data) {
  data$row.orig <- 1:(dim(data)[1])
  # Use DR.Type as a starting point - treating all dichotomous as quantal-deterministic
  data$Conceptual.Model <- data$DR.Type
  if (sum(is.dichot<-data$DR.Type == "Dichotomous")>0) {
    data$Conceptual.Model[is.dichot] <- "Quantal-Deterministic"
  }
  # Modify based on effect type
  effectvec <- data$Effect
  if (sum(is.quant.deter <- effectvec == "nonneoplastic histopathology" |
          effectvec == "clinical signs" |
          effectvec == "gross pathology") > 0) {
    data$Conceptual.Model[is.quant.deter] <- "Quantal-Deterministic"
  }
  if (sum(is.cont <- effectvec == "body weight" |
          effectvec == "organ weight" |
          effectvec == "clinical chemistry" |
          effectvec == "urinalysis" |
          effectvec == "enzyme activity" |
          effectvec == "hematology" |
          effectvec == "food and/or water consumption" |
          effectvec == "neurotransmitter") > 0) {
    data$Conceptual.Model[is.cont] <- "Continuous"
  }
  # For reproductive and developmental effects, treat dichotomous effects as quantal-stochastic
  if (sum(is.quant.stoch <- data$DR.Type == "Dichotomous" & 
          (effectvec == "reproduction" | effectvec == "development"|
           effectvec == "mortality/survival" )) > 0) {
    data$Conceptual.Model[is.quant.stoch] <- "Quantal-Stochastic"
  }
  # Mortality/survival -- Quantal-stochastic
  data$Conceptual.Model[effectvec == "mortality/survival"] <- "Quantal-Stochastic"
  # For unspecified developmental and reproductive, need to do both Continuous and Quantal-Stochastic
  if (sum(multirows.1 <- data$DR.Type == "" & 
          (effectvec == "development" | effectvec == "reproduction"))>0) {
    data.newrows1 <- data[multirows.1,]
    data$Conceptual.Model[multirows.1] <- "Continuous"
    data.newrows1$Conceptual.Model <- "Quantal-Stochastic"
  }
  # For unspecified none, multiple, or other, need to do both Continuous and Quantal-Deterministic
  if (sum(multirows.2 <- data$DR.Type == "" & 
          (effectvec == "none" | effectvec == "multiple" |
                                       effectvec == "other"))>0) {
    data.newrows2 <- data[multirows.2,]
    data$Conceptual.Model[multirows.2] <- "Continuous"
    data.newrows2$Conceptual.Model <- "Quantal-Deterministic"
  }
  # 
  if (sum(multirows.1) > 0) {
    data <- rbind(data,data.newrows1)
  }
  if (sum(multirows.2) > 0) {
    data <- rbind(data,data.newrows2)
  }
  data <- data[order(data$row.orig),]
  data
}


##################### Add Magnitude of Effect

addEffect.Magnitude <- function(data) {
  # Begin with BMR Type 
  data$Effect.Magnitude <- data$BMR.Type
  # For 1SD change, convert to % RD
  is1SD <- data$Effect.Magnitude == "1SD"
  is1SD[is.na(is1SD)] <- FALSE
  isSE <- data$CV.Type == "SE"
  isSE[is.na(isSE)] <- FALSE
  # Convert to % Change, using SD  / mean response
  pctchange <- data$numSDResponse/data$numMeanResponse*100
  # Convert SE to SD if necessary
  if (sum(isSE) > 0) {
    pctchange[isSE] <- pctchange[isSE] * sqrt(data$intAnimalsControl[isSE])
  }
  data$Effect.Magnitude[is1SD] <- paste(round(pctchange[is1SD]),"% RD",sep="")
  # For LOAELS or NOAELS, fixed by Assessment Factor method 
  isnotbmdl <- data$POD.type != "BMDL"
  iscontinuous <- data$Conceptual.Model == "Continuous"
  data$Effect.Magnitude[isnotbmdl & iscontinuous] <- "5% RD"
  isquantstoch <- data$Conceptual.Model == "Quantal-Stochastic"
  data$Effect.Magnitude[isnotbmdl & isquantstoch] <- "10% ER"
  # For Quant-Deterministic, always ED50
  isquantdeter <- data$Conceptual.Model == "Quantal-Deterministic"
  data$Effect.Magnitude[isquantdeter] <- "ED50"
  data
}


##################### Apply probablistic Adjustment factors - lognormal approximation

addPrAF.BMDL <- function(data, p95p50.default=3) { 
  data$PrAF.BMDL.p50 <- 1
  data$PrAF.BMDL.p95p50 <- 1
  isbmdl <- data$POD.type == "BMDL"
  # First use default value of BMD/BMDL ratio
  data$PrAF.BMDL.p50[isbmdl] <- 1/p95p50.default
  data$PrAF.BMDL.p95p50[isbmdl] <- p95p50.default
  # Adjust to ED50 if quantal-deterministic and effect magnitude not 50%
  isquantaldeterm <- data$Conceptual.Model == "Quantal-Deterministic"
  ised50 <- data$BMR.Type == "ED50"
  data$PrAF.BMDL.p50[isbmdl & isquantaldeterm & !ised50] <- 
    data$PrAF.BMDL.p50[isbmdl & isquantaldeterm & !ised50]/3
  data$PrAF.BMDL.p95p50[isbmdl & isquantaldeterm & !ised50] <- 
    exp(sqrt(log(data$PrAF.BMDL.p95p50[isbmdl & isquantaldeterm & !ised50])^2 + log(1.5)^2))
  # Use actual value of BMD/BMDL ratio if present
  isbmd <- !is.na(data$numBMD)
  data$PrAF.BMDL.p50[isbmdl & isbmd] <- data$numPOD[isbmdl & isbmd]/data$numBMD[isbmdl & isbmd]
  data$PrAF.BMDL.p95p50[isbmdl & isbmd] <- data$numBMD[isbmdl & isbmd]/data$numPOD[isbmdl & isbmd]
  # Adjust to ED50 if quantal-deterministic and effect magnitude not 50%
  data$PrAF.BMDL.p50[isbmdl & isbmd & isquantaldeterm & !ised50] <- 
    data$PrAF.BMDL.p50[isbmdl & isbmd & isquantaldeterm & !ised50]/3
  data$PrAF.BMDL.p95p50[isbmdl & isbmd & isquantaldeterm & !ised50] <- 
    exp(sqrt(log(data$PrAF.BMDL.p95p50[isbmdl & isbmd & isquantaldeterm & !ised50])^2 + log(1.5)^2))
  # Return data
  data
}

addPrAF.LOAEL <- function(data, p95p50.default=3) {
  data$PrAF.LOAEL.p50 <- 1
  data$PrAF.LOAEL.p95p50 <- 1
  isloael <- data$POD.type == "LOAEL" | data$POD.type == "LEL" | data$POD.type == "LOEL"
  data$PrAF.LOAEL.p50[isloael] <- data$numUFl[isloael]
  data$PrAF.LOAEL.p95p50[isloael] <- p95p50.default
  data
}

addPrAF.NOAEL <- function(data) { 
  data$PrAF.NOAEL.p50 <- 1
  data$PrAF.NOAEL.p95p50 <- 1
  isnoael <- data$POD.type == "NOAEL" | data$POD.type == "NOEL" | 
    data$POD.type == "LOAEL" | data$POD.type == "LEL" | data$POD.type == "LOEL"
  iscontinuous <- data$Conceptual.Model == "Continuous"
  isquantdeterm <- data$Conceptual.Model == "Quantal-Deterministic"
  isquantstoch <- data$Conceptual.Model == "Quantal-Stochastic"
  isreprodev <- (data$Effect == "reproduction" | data$Effect == "development")
  # Continuous - not repro/dev
  data$PrAF.NOAEL.p50[isnoael & iscontinuous & !isreprodev] <- 1/3
  data$PrAF.NOAEL.p95p50[isnoael & iscontinuous & !isreprodev] <- 4.7
  # Continuous - repro/dev
  data$PrAF.NOAEL.p50[isnoael & iscontinuous & isreprodev] <- 1/3
  data$PrAF.NOAEL.p95p50[isnoael & iscontinuous & isreprodev] <- 7
  # Quantal-deterministic
  data$PrAF.NOAEL.p50[isnoael & isquantdeterm] <- 2/9
  data$PrAF.NOAEL.p95p50[isnoael & isquantdeterm] <- 5
  # Quantal-stochastic
  data$PrAF.NOAEL.p50[isnoael & isquantstoch] <- 2/3
  data$PrAF.NOAEL.p95p50[isnoael & isquantstoch] <- 4.7
  data
}

addPrAF.Subchronic <- function(data) {
  issubchronic <- data$numUFs > 1
  data$PrAF.Subchronic.p50 <- 1
  data$PrAF.Subchronic.p95p50 <- 1
  data$PrAF.Subchronic.p50[issubchronic] <- 2
  data$PrAF.Subchronic.p95p50[issubchronic] <- 4
  data
}

addPrAF.BW <- function(data, BW.human=70, 
                       # All EPA 1988 values are for Chronic
                       BW.mouse=0.0363, # Mean of M and F B6C3F1 from EPA 1988
                       BW.rat=0.3045, # Mean of M and F F344 from EPA 1988
                       BW.rabbit=3.845, # Mean of M and F NZ
                       BW.dog=10.45, # Mean of M and F beagles
                       BW.primate=9.45 # Mean of Mand F Rhesus
) {
  data$PrAF.BW.p50 <- 1
  data$PrAF.BW.p95p50 <- 1
  ismouse <- data$Species == "Mouse"
  data$PrAF.BW.p50[ismouse] <- (BW.human/BW.mouse)^0.3
  data$PrAF.BW.p95p50[ismouse] <- (BW.human/BW.mouse)^0.04
  israt <- data$Species == "Rat"
  data$PrAF.BW.p50[israt] <- (BW.human/BW.rat)^0.3
  data$PrAF.BW.p95p50[israt] <- (BW.human/BW.rat)^0.04
  israbbit <- data$Species == "Rabbit"
  data$PrAF.BW.p50[israbbit] <- (BW.human/BW.rabbit)^0.3
  data$PrAF.BW.p95p50[israbbit] <- (BW.human/BW.rabbit)^0.04
  isdog <- data$Species == "Dog"
  data$PrAF.BW.p50[isdog] <- (BW.human/BW.dog)^0.3
  data$PrAF.BW.p95p50[isdog] <- (BW.human/BW.dog)^0.04
  isprimate <- data$Species == "Non-human primate"
  data$PrAF.BW.p50[isprimate] <- (BW.human/BW.primate)^0.3
  data$PrAF.BW.p95p50[isprimate] <- (BW.human/BW.primate)^0.04
  isbw <- !is.na(data$numBW)
  data$PrAF.BW.p50[isbw] <- (BW.human/data$numBW[isbw])^0.3
  data$PrAF.BW.p95p50[isbw] <- (BW.human/data$numBW[isbw])^0.04
  data
}

addPrAF.TKTD <- function(data) {
  ishuman <- data$Species == "Human"
  data$PrAF.TKTD.p50 <- 1
  data$PrAF.TKTD.p95p50 <- 1
  data$PrAF.TKTD.p95p50[!ishuman] <- 3
  data
}

addHDM50.approx <- function(data) {
  data$HDM50.approx.p50 <- data$numPOD / (data$PrAF.BMDL.p50 * 
                                           data$PrAF.LOAEL.p50 * 
                                           data$PrAF.NOAEL.p50 * 
                                           data$PrAF.Subchronic.p50 * 
                                           data$PrAF.BW.p50 * 
                                           data$PrAF.TKTD.p50)
  data$HDM50.approx.p95p50 <- exp(sqrt(log(data$PrAF.BMDL.p95p50)^2+
                                        log(data$PrAF.LOAEL.p95p50)^2+
                                        log(data$PrAF.NOAEL.p95p50)^2+
                                        log(data$PrAF.Subchronic.p95p50)^2+
                                        log(data$PrAF.BW.p95p50)^2+
                                        log(data$PrAF.TKTD.p95p50)^2
  ))
  data
}

addPrAF.Human.approx <- function(data, incidence=0.01, log10gsdh=0.42390, log10gsdu=0.16621) {
  # uses normal approximation for sigmaH that has same p05 and p95 as lognormal sigmaH
  # lognormal sigmaH 
  z <- qnorm(1-incidence)
  data$Human.log10gsdh.approx <- log10gsdh
  data$Human.log10gsdu.approx <- log10gsdu
  log10gsdh.p05 <- log10(exp(log(10^(log10gsdh+log10gsdu*qnorm(0.05)))))
  log10gsdh.p95 <- log10(exp(log(10^(log10gsdh+log10gsdu*qnorm(0.95)))))
  data$Human.log10gsdh.p50 <- sqrt(log10gsdh.p05*log10gsdh.p95)
  data$Human.log10gsdh.p95p50 <- sqrt(log10gsdh.p95/log10gsdh.p05)
  data$Human.approx.incidence.target <- incidence
  data$PrAF.Human.approx.p50 <- exp(log10gsdh*z*log(10))
  data$PrAF.Human.approx.p95p50 <- exp(log10gsdu*z*log(10)*qnorm(0.95))
  data
}

addHDMI.approx <- function(data) {
  data$HDMI.approx.p50 <- data$numPOD / (data$PrAF.BMDL.p50 * 
                                           data$PrAF.LOAEL.p50 * 
                                           data$PrAF.NOAEL.p50 * 
                                           data$PrAF.Subchronic.p50 * 
                                           data$PrAF.BW.p50 * 
                                           data$PrAF.TKTD.p50 * 
                                           data$PrAF.Human.approx.p50)
  data$HDMI.approx.p95p50 <- exp(sqrt(log(data$PrAF.BMDL.p95p50)^2+
                                        log(data$PrAF.LOAEL.p95p50)^2+
                                        log(data$PrAF.NOAEL.p95p50)^2+
                                        log(data$PrAF.Subchronic.p95p50)^2+
                                        log(data$PrAF.BW.p95p50)^2+
                                        log(data$PrAF.TKTD.p95p50)^2+
                                        log(data$PrAF.Human.approx.p95p50)^2
  ))
  data
}

addPrRfD.approx <- function(data) {
  data$PrRfD.approx <- data$HDMI.approx.p50/data$HDMI.approx.p95p50
  data$PrRfDText <- ""
  isquantdet<-data$Conceptual.Model == "Quantal-Deterministic"
  effectandeffectPODsame <- data$strEffectPOD == data$Effect
  
  data$PrRfDText[isquantdet] <- 
    paste("Probabilistic RfD of ",signif(data$PrRfD.approx[isquantdet],3),
          " mg/kg-d = estimate, based on the study data, of the dose at which, with 95% confidence, at most ",
          data$Human.approx.incidence.target[isquantdet]*100,
          "% of the population will have \'",
          data$strEffectPOD[isquantdet]," (",data$Organ[isquantdet],
          " effect: ",data$Effect[isquantdet],")\'.",sep="")

  data$PrRfDText[isquantdet & effectandeffectPODsame] <- 
    paste("Probabilistic RfD of ",signif(data$PrRfD.approx[isquantdet &
                                                           effectandeffectPODsame],3),
          " mg/kg-d = estimate, based on the study data, of the dose at which, with 95% confidence, at most ",
          data$Human.approx.incidence.target[isquantdet & effectandeffectPODsame]*100,
          "% of the population will have \'",data$Organ[isquantdet & effectandeffectPODsame],
          " effect: ",data$Effect[isquantdet & effectandeffectPODsame],"\'.",sep="")

  data$PrRfDText[!isquantdet]<- 
    paste("Probabilistic RfD of ",signif(data$PrRfD.approx[!isquantdet],3),
          " mg/kg-d = estimate, based on the study data, of the dose at which, with 95% confidence, at most ",
          data$Human.approx.incidence.target[!isquantdet]*100,
          "% of the population will have \'",
          data$strEffectPOD[!isquantdet]," (",data$Organ[!isquantdet],
          " effect: ",data$Effect[!isquantdet],")\' of magnitude ",
          data$Effect.Magnitude[!isquantdet]," or greater.",sep="")

  data$PrRfDText[!isquantdet & effectandeffectPODsame]<- 
    paste("Probabilistic RfD of ",signif(data$PrRfD.approx[!isquantdet & effectandeffectPODsame],3),
          " mg/kg-d = estimate, based on the study data, of the dose at which, with 95% confidence, at most ",
          data$Human.approx.incidence.target[!isquantdet & effectandeffectPODsame]*100,
          "% of the population will have \'",data$Organ[!isquantdet & effectandeffectPODsame],
          " effect: ",data$Effect[!isquantdet & effectandeffectPODsame],"\' of magnitude ",
          data$Effect.Magnitude[!isquantdet & effectandeffectPODsame]," or greater.",sep="")
  
  data
}

addHDMI.mc <- function(data, n = 1e7, seed=3.14159, showprogress=FALSE) {
 set.seed(seed)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  for (i in 1:nrow(data)) {
    if (showprogress) cat(i,"...",sep="")
    hdm50.mc <- data$HDM50.approx.p50[i] * 
      (data$HDM50.approx.p95p50[i])^(z1/qnorm(0.95))
    sigmah.mc <- log(10) * data$Human.log10gsdh.p50[i] * 
      (data$Human.log10gsdh.p95p50[i])^(z2/qnorm(0.95))
    afh.mc <- exp(sigmah.mc*qnorm(1-data$Human.approx.incidence.target[i]))
    hdmi.mc <- hdm50.mc / afh.mc
    quants <- quantile(hdmi.mc,prob=c(0.05,0.5,0.95))
    data$HDMI.mc.p50[i] <- quants[2]
    data$HDMI.mc.p05[i] <- quants[1]
    data$HDMI.mc.p95[i] <- quants[3]
  }
  if (showprogress) cat("\n")
  data
}

# adddr.slope.E.bgrd.mc <- function(data, n = 1e7, seed=3.14159) {
#   
# # Background exposure: B50.p50	B50.p95p50	
# # Variability in background exposure: B95B50.p50	B95B50.95p50
# #  z = (log B50 - log HD50) / sigmaTot
# #  sigmaTot = sqrt(sigmaH^2 + sigmaB^2)
# #  slope = dnorm(z) / (B50 * sigmaTot)
# #  EF = slope * 1e6 / (t in days * BW in kg) - do separately
#   
#   set.seed(seed)
#   z1 <- rnorm(n)
#   z2 <- rnorm(n)
#   z3 <- rnorm(n)
#   z4 <- rnorm(n)
#   for (i in 1:nrow(data)) {
#     hdm50.mc <- data$HDM50.approx.p50[i] * 
#       (data$HDM50.approx.p95p50[i])^(z1/qnorm(0.95))
#     sigmah.mc <- log(10) * data$Human.log10gsdh.p50[i] * 
#       (data$Human.log10gsdh.p95p50[i])^(z2/qnorm(0.95))
#     B50.mc <- data$B50.p50[i] * data$B50.p95p50[i]^(z3/qnorm(0.95))
#     sigmaB.mc <- log(data$B95B50.p50[i] * data$B95B50.p95p50[i]^(z4/qnorm(0.95)))/qnorm(0.95)
#     sigmatot.mc <- sqrt(sigmah.mc^2+sigmaB.mc^2)
#     z.mc <- (log(B50.mc) - log(hdm50.mc))/sigmatot.mc
#     dr.slope.E.bgrd.mc <- dnorm(z.mc)/(B50.mc * sigmatot.mc)
#     quants <- quantile(dr.slope.E.bgrd.mc,prob=c(0.05,0.5,0.95))
#     data$dr.slope.E.bgrd.mc.p50[i] <- quants[2]
#     data$dr.slope.E.bgrd.mc.p05[i] <- quants[1]
#     data$dr.slope.E.bgrd.mc.p95[i] <- quants[3]
#     data$dr.slope.E.bgrd.mc.mean[i] <- mean(dr.slope.E.bgrd.mc)
#   }
#   data
# }
# 
# adddr.slope.I.bgrd.mc <- function(data, n = 1e7, seed=3.14159, addrisk=0.01) {
#   
# #  zeff defined by qnorm(background incidence)
# #  Beff = HD50 * exp(zeff*sigmaH)
# #  slope = dnorm(zeff) / (HD50 * sigmaH * exp(zeff*sigmaH))
# 
#   set.seed(seed)
#   z1 <- rnorm(n)
#   z2 <- rnorm(n)
#   z3 <- rnorm(n)
#   for (i in 1:nrow(data)) {
#     hdm50.mc <- data$HDM50.approx.p50[i] * 
#       (data$HDM50.approx.p95p50[i])^(z1/qnorm(0.95))
#     sigmah.mc <- log(10) * data$Human.log10gsdh.p50[i] * 
#       (data$Human.log10gsdh.p95p50[i])^(z2/qnorm(0.95))
#     z.eff.mc <- qnorm(data$I.bgrd.p50[i] * data$I.bgrd.p95p50[i]^(z3/qnorm(0.95)))
#     Beff.mc <- hdm50.mc*exp(z.eff.mc*sigmah.mc)
#     dr.slope.I.bgrd.mc <- dnorm(z.eff.mc)/(Beff.mc * sigmah.mc)
#     
#     prob.rfd.additive.fromslope.mc <- addrisk/dr.slope.I.bgrd.mc
#     zstar.mc <- qnorm(pnorm(z.eff.mc) + addrisk)
#     prob.rfd.additive.exact.mc <- hdm50.mc*exp(sigmah.mc*zstar.mc) - Beff.mc
# 
#     quants <- quantile(prob.rfd.additive.fromslope.mc,prob=c(0.05,0.5,0.95))
#     data$prob.rfd.additive.fromslope.mc.p50[i] <- quants[2]
#     data$prob.rfd.additive.fromslope.mc.p05[i] <- quants[1]
#     data$prob.rfd.additive.fromslope.mc.p95[i] <- quants[3]
#     
#     quants <- quantile(prob.rfd.additive.exact.mc,prob=c(0.05,0.5,0.95))
#     data$prob.rfd.additive.exact.mc.p50[i] <- quants[2]
#     data$prob.rfd.additive.exact.mc.p05[i] <- quants[1]
#     data$prob.rfd.additive.exact.mc.p95[i] <- quants[3]
#     
#     quants <- quantile(dr.slope.I.bgrd.mc,prob=c(0.05,0.5,0.95))
#     data$dr.slope.I.bgrd.mc.p50[i] <- quants[2]
#     data$dr.slope.I.bgrd.mc.p05[i] <- quants[1]
#     data$dr.slope.I.bgrd.mc.p95[i] <- quants[3]
#     data$dr.slope.I.bgrd.mc.mean[i] <- mean(dr.slope.I.bgrd.mc)
#   }
#   data
# }

#####

HDMI.approx.fun <- function(incidence,prob=0.05,datarow=datarow) {
  z <- qnorm(1-incidence)
  PrAF.Human.approx.p50.z <- exp(0.42390*z*log(10))
  PrAF.Human.approx.p95p50.z <- exp(0.16621*z*log(10)*qnorm(0.95))
  HDMI.approx.p50 <- datarow$numPOD / (datarow$PrAF.BMDL.p50 * 
                                         datarow$PrAF.LOAEL.p50 * 
                                         datarow$PrAF.NOAEL.p50 * 
                                         datarow$PrAF.Subchronic.p50 * 
                                         datarow$PrAF.BW.p50 * 
                                         datarow$PrAF.TKTD.p50 * 
                                         PrAF.Human.approx.p50.z)
  HDMI.approx.p95p50 <- exp(sqrt(log(datarow$PrAF.BMDL.p95p50)^2+
                                   log(datarow$PrAF.LOAEL.p95p50)^2+
                                   log(datarow$PrAF.NOAEL.p95p50)^2+
                                   log(datarow$PrAF.Subchronic.p95p50)^2+
                                   log(datarow$PrAF.BW.p95p50)^2+
                                   log(datarow$PrAF.TKTD.p95p50)^2+
                                   log(PrAF.Human.approx.p95p50.z)^2
  ))
  HDMI.approx.p <- HDMI.approx.p50*(HDMI.approx.p95p50)^(qnorm(prob)/qnorm(0.95))
  HDMI.approx.p
}

add.zRfD <- function(data) {
  data$zRfD <-(log(data$RfDEndpoint)-log(data$HDMI.approx.p50))/log((data$HDMI.approx.p95p50)^(1/qnorm(0.95)))
  data
}

add.I.approx.at.RfD <- function(data,HQ=1) {
  numPrRfD <- nrow(data)
  I95<-numeric(numPrRfD)
  I50<-numeric(numPrRfD)
  I05<-numeric(numPrRfD)
  for (i in 1:numPrRfD) {
    if (HDMI.approx.fun(1e-15,datarow=data[i,]) < (HQ*data$RfDEndpoint[i])) {
      if (HDMI.approx.fun(10^(-1e-15),datarow=data[i,]) < (HQ*data$RfDEndpoint[i])) {
        I95[i]<-1
      } else {
        I95[i]<-10^uniroot(function(x) {HDMI.approx.fun(10^x,datarow=data[i,]) - 
            HQ*data$RfDEndpoint[i]},lower=-15,upper=-1e-15)$root
      }
    }
    if (HDMI.approx.fun(1e-15,prob=0.5,datarow=data[i,]) < (HQ*data$RfDEndpoint[i])) {
      if (HDMI.approx.fun(10^(-1e-15),prob=0.5,datarow=data[i,]) < (HQ*data$RfDEndpoint[i])) {
        I50[i]<-1
      } else {
        I50[i]<-10^uniroot(function(x) {HDMI.approx.fun(10^x,prob=0.5,datarow=data[i,]) - 
            HQ*data$RfDEndpoint[i]},lower=-15,upper=-1e-15)$root
      }
    }
    if (HDMI.approx.fun(1e-15,prob=0.95,datarow=data[i,]) < (HQ*data$RfDEndpoint[i])) {
      if (HDMI.approx.fun(10^(-1e-15),prob=0.95,datarow=data[i,]) < (HQ*data$RfDEndpoint[i])) {
        I05[i]<-1
      } else {
        I05[i]<-10^uniroot(function(x) {HDMI.approx.fun(10^x,prob=0.95,datarow=data[i,]) - 
            HQ*data$RfDEndpoint[i]},lower=-15,upper=-1e-15)$root
      }
    }
  }
  if (HQ==1) {
    data$I95.approx.RfD <- I95
    data$I50.approx.RfD <- I50
    data$I05.approx.RfD <- I05
  } else {
    data[[paste("I95.approx.HQ",HQ,sep="")]]<-I95
    data[[paste("I50.approx.HQ",HQ,sep="")]]<-I50
    data[[paste("I05.approx.HQ",HQ,sep="")]]<-I05
  }
  data
}


zDRfunc.bound <- function(dose, prob=0.95,datarow=datarow, sigmah=0.42390, sigmau=0.16621) {
  zp <- qnorm(prob)
  z95 <- qnorm(0.95)
  y <- log(dose * datarow$PrAF.BMDL.p50 *
             datarow$PrAF.LOAEL.p50 *
             datarow$PrAF.NOAEL.p50 *
             datarow$PrAF.Subchronic.p50 *
             datarow$PrAF.BW.p50 *
             datarow$PrAF.TKTD.p50 / datarow$numPOD)
  v2 <- (log(datarow$PrAF.BMDL.p95p50)^2+
           log(datarow$PrAF.LOAEL.p95p50)^2+
           log(datarow$PrAF.NOAEL.p95p50)^2+
           log(datarow$PrAF.Subchronic.p95p50)^2+
           log(datarow$PrAF.BW.p95p50)^2+
           log(datarow$PrAF.TKTD.p95p50)^2)/z95^2
  zi <- (y+zp*sqrt(sigmau^2/sigmah^2*y^2 + v2*(1-zp^2*sigmau^2/sigmah^2)))/
    ((1-zp^2*sigmau^2/sigmah^2)*log(10)*sigmah)
  zi
}

addVarContrib <- function(data) {
  data$PrAF.BMDL.Var <- log(data$PrAF.BMDL.p95p50)^2 
  data$PrAF.LOAEL.Var <- log(data$PrAF.LOAEL.p95p50)^2 
  data$PrAF.NOAEL.Var <- log(data$PrAF.NOAEL.p95p50)^2 
  data$PrAF.Subchronic.Var <- log(data$PrAF.Subchronic.p95p50)^2 
  data$PrAF.BW.Var <- log(data$PrAF.BW.p95p50)^2 
  data$PrAF.TKTD.Var <- log(data$PrAF.TKTD.p95p50)^2 
  data$PrAF.Human.approx.Var <- log(data$PrAF.Human.approx.p95p50)^2 
  data
}

addPctContrib <- function(data) {
  data$PrAF.BMDL.pct <- log(data$PrAF.BMDL.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data$PrAF.LOAEL.pct <- log(data$PrAF.LOAEL.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data$PrAF.NOAEL.pct <- log(data$PrAF.NOAEL.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data$PrAF.Subchronic.pct <- log(data$PrAF.Subchronic.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data$PrAF.BW.pct <- log(data$PrAF.BW.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data$PrAF.TKTD.pct <- log(data$PrAF.TKTD.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data$PrAF.Human.approx.pct <- log(data$PrAF.Human.approx.p95p50)^2 / log(data$HDMI.approx.p95p50)^2
  data
}
