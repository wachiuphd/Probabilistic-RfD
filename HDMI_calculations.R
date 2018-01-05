source("HDMI_functions.R")

### Load database
# Read input file
RfD_data_in<-read.csv("RfD_InputData.csv",as.is=TRUE)
# Number of RfDs & Endpoints
nRfdEndpoint <- nrow(RfD_data_in)
# Number of unique chemicals
nChem <- length(unique(RfD_data_in$strName))
# Number of chemicals with multiple RfDs or endpoints
nChemMulti <- sum(table(RfD_data_in$strName) > 1)
# Distribution of endpoints
tabEndpoints <- table(RfD_data_in$Effect)
write.csv(tabEndpoints,"RfDEffectTypes.csv")
# Cross-tabulation of PODs by Source
tabPODs <- table(RfD_data_in[,c("Source","POD.type")])
write.csv(tabPODs,"RfDPODs.csv")
# Cross-tabulation of composite UF by Source
tabUFcomps <- table(RfD_data_in[,c("numUF","Source")])
write.csv(tabUFcomps,"RfDUFcomps.csv")

### STEP 1: Convert to endpoint-specific RfDs by removing database factor
RfD_data_Endpoint <- RfD_data_in
# Set missing UF values to 1
RfD_data_Endpoint$numUFa[is.na(RfD_data_Endpoint$numUFa)] <- 1
RfD_data_Endpoint$numUFh[is.na(RfD_data_Endpoint$numUFh)] <- 1
RfD_data_Endpoint$numUFs[is.na(RfD_data_Endpoint$numUFs)] <- 1
RfD_data_Endpoint$numUFl[is.na(RfD_data_Endpoint$numUFl)] <- 1
RfD_data_Endpoint$numUFd[is.na(RfD_data_Endpoint$numUFd)] <- 1
# Identify which RfDs have database factor
hasUFd <- RfD_data_Endpoint$numUFd > 1
numhasUFd <- sum(hasUFd)
# Calculate endpoints-specific RfDs and new composite UF
RfD_data_Endpoint$RfDEndpoint <- RfD_data_Endpoint$numValue * RfD_data_Endpoint$numUFd
RfD_data_Endpoint$UFEndpoint <- RfD_data_Endpoint$numUF / RfD_data_Endpoint$numUFd
# Save output as csv file
write.csv(RfD_data_Endpoint,file="RfD_Endpoint.csv",row.names=FALSE)

### STEP 2: Assign conceptual model(s) and magnitudes of 
### effect to each endpoint-specific RfD
RfD_data_out <- addConceptual.Model(RfD_data_Endpoint)
RfD_data_out <- addEffect.Magnitude(RfD_data_out)

### STEP 3: Assign uncertainty distributions for each POD and UF
# Adjustments to POD
RfD_data_out <- addPrAF.BMDL(RfD_data_out)
RfD_data_out <- addPrAF.LOAEL(RfD_data_out)
RfD_data_out <- addPrAF.NOAEL(RfD_data_out)
# Adjustment for subchronic study
RfD_data_out <- addPrAF.Subchronic(RfD_data_out)
# Adjustment for interspecies
RfD_data_out <- addPrAF.BW(RfD_data_out, BW.human = 60) # Modified to WHO default
RfD_data_out <- addPrAF.TKTD(RfD_data_out)
# Adjustment for human variability - Approximate calculation (analytic solution)
RfD_data_out <- addPrAF.Human.approx(RfD_data_out)

### STEP 4: Combine POD and UF uncertainties probabilistically
# Median and P95/P50 of HDM50 (dose for median individual)
RfD_data_out <- addHDM50.approx(RfD_data_out)
# Median and P95/P50 of HDMI (dose for incidence I)
RfD_data_out <- addHDMI.approx(RfD_data_out)
nHDMI <- nrow(RfD_data_out)
# Probabilistic RfD = P05 of HDMI
RfD_data_out <- addPrRfD.approx(RfD_data_out)
# z-score at RfD
RfD_data_out <- add.zRfD(RfD_data_out)
# P05,median,P95 of incidence at the endpoint-specific RfD
RfD_data_out <- add.I.approx.at.RfD(RfD_data_out)
# P05,median,P95 of incidence at the traditional HQ=3
RfD_data_out <- add.I.approx.at.RfD(RfD_data_out,HQ=3)
# P05,median,P95 of incidence at the traditional HQ=10
RfD_data_out <- add.I.approx.at.RfD(RfD_data_out,HQ=10)
# P05,median,P95 of incidence at the traditional HQ=0.1
RfD_data_out <- add.I.approx.at.RfD(RfD_data_out,HQ=0.1)
# P05,median,P95 of incidence at the traditional HQ=0.3
RfD_data_out <- add.I.approx.at.RfD(RfD_data_out,HQ=0.3)


# Contributions to variance
RfD_data_out <- addVarContrib(RfD_data_out)
RfD_data_out <- addPctContrib(RfD_data_out)

### Save output as csv file
write.csv(RfD_data_out,file="RfD_HDMI_results.csv",row.names=FALSE)

###
RfD_data_out_mc <- addHDMI.mc(RfD_data_out,showprogress=TRUE)
write.csv(RfD_data_out_mc,"RfD_HDMI_mc_results.csv",row.names=FALSE)
