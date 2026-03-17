##############################################################################################
#' @title CREATE A 3-PANEL FIGURE SHOWING THE RATING CURVE CONTROLS, POSTERIOR
#' RATING CURVE, AND CONTINUOUS DISCHARGE TIMESERIES FOR A NEON SITE

#' @author
#' Zachary Nickerson \email{nickerson@battelleecology.org} \cr

#' @description First, the hydraulic controls are plotted that define the prior
#' rating curve model at NEON sites. The controls are plotted overlaying the 
#' discharge cross-section data. Second, the posterior rating curve is modelled
#' and plotted with uncertainty and the empirical stage/discharge gaugings that
#' define the posterior model. Third, the resulting continuous discharge is
#' plotted with uncertainty and empirical discharge bouts. The three plots are 
#' combined into a single multi-panel plot.

#' @details
#' Users should run this script within the .Rproj file in the root directory.
#' 
#' It is recommended that downloaders of NEON data have a NEON Data Portal
#' personal access token (PAT) set as an environment variable - 
#' Sys.setenv(NEON_PAT="YOUR PAT"). See instructions on obtaining a PAT at 
#' https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-tokens-tutorial 
#' 
#' Plotting the hydraulic controls cannot be done solely with data available
#' from the portal. Some knowledge is required of survey point IDs that were 
#' used to delineate the left and right extent of each control, which is not
#' included in the published tables. To determine the survey point IDs needed
#' for the D07 WALK cross-section survey completed on 2022-09-15, the following
#' script was referenced: 
#' https://github.com/NEONScience/NEON-stream-discharge/blob/main/hydrologicControls/D07_WALK/WALK_controls_20220915.R

#' @return NEON Portal downloads are saved as .rds to the /data subdirectory
#' @return Output plot (.png) is saved to the /out subdirectory

# changelog and author contributions / copyrights
#   Zachary Nickerson (2026-03-16)
#     original creation
##############################################################################################

# Load libraries ####
library(neonUtilities)
library(dplyr)
library(ggplot2)
library(cowplot)
require(stageQCurve)

# Set options ####
options(stringsAsFactors=F,
        scipen = 999)

# Set constants ####
queryStartDate <- "2023-10"
queryEndDate <- "2024-09"
queryRelease <- "RELEASE-2026"

# Set filepath to BaM executable (see header) ####
DIRPATH <- "C:/Users/nickerson/Documents/GitHub/NEON-hydro-clean-fill-proc/shiny-cleanFlow/BaM_beta/"
BAMWS <- "BaM_BaRatin/"

# Create subdirectories if they do not exist ####
if(!dir.exists("./data")){dir.create("./data")}
if(!dir.exists("./out")){dir.create("./out")}

# Set output filepath for final plot
plotSavePath <- paste0(getwd(),"/out/walk_fig_",queryRelease,".png")

# Panel 1: Discharge cross section with controls ####

## Download RELEASED stage-discharge rating curve data ####
sdrc_WALK <- neonUtilities::loadByProduct(
  dpID = "DP4.00133.001",
  site = "WALK",
  release = queryRelease,
  package = "expanded",
  check.size = F,
  token = Sys.getenv("NEON_PAT")
)
# Save downloaded data as RDS for quicker access in future runs
saveRDS(sdrc_WALK,paste0("./data/NEON.DP4.00133.001_WALK_",queryRelease,".rds"))

# If data is saved locally, read it in
# sdrc_WALK <- readRDS(paste0("./data/NEON.DP4.00133.001_WALK_",queryRelease,".rds"))

## Get the date of the cross-section survey from curveIndentification ####
sdrc_curveIdentification <- sdrc_WALK$sdrc_curveIdentification
sdrc_curveIdentification <- sdrc_curveIdentification[
  sdrc_curveIdentification$curveID=="WALK.2024",
]
surveyDate <- format(sdrc_curveIdentification$controlSurveyEndDateTime,
                     "%Y-%m")

## Download RLEASED stream morphology data to retrieve survey points ####
geo_WALK <- neonUtilities::loadByProduct(
  dpID = "DP4.00131.001",
  site = "WALK",
  startdate = surveyDate,
  enddate = surveyDate,
  release = queryRelease,
  package = "expanded",
  check.size = F,
  token = Sys.getenv("NEON_PAT")
)
# Save downloaded data as RDS for quicker access in future runs
saveRDS(geo_WALK,paste0("./data/NEON.DP4.00131.001_WALK_",queryRelease,".rds"))

# If data is saved locally, read it in
# geo_WALK <- readRDS(paste0("./data/NEON.DP4.00131.001_WALK_",queryRelease,".rds"))

## Derive discharge cross section from relative easting, northing, height ####
geo_surveyPoints <- geo_WALK$geo_surveyPoints
geo_surveyPoints_dsc <- geo_surveyPoints[
  geo_surveyPoints$mapCode=="Transect_DSC",
]
geo_surveyPoints_dsc<-geo_surveyPoints_dsc[
  order(geo_surveyPoints_dsc$relativeNorthing),
]
geo_surveyPoints_gag <- geo_surveyPoints[
  geo_surveyPoints$mapCode=="Gauge",
]
rownames(geo_surveyPoints_dsc)<-seq(length=nrow(geo_surveyPoints_dsc))
leftPoint<-"DSC_LFH" # Manual identification of left-most survey point

### Assign raw distance to each point relative to left-most survey point ####
geo_surveyPoints_dsc$distanceRaw <- NA
for(i in 1:nrow(geo_surveyPoints_dsc)){
  pointNorth<-geo_surveyPoints_dsc$relativeNorthing[i]
  pointEast<-geo_surveyPoints_dsc$relativeEasting[i]
  geo_surveyPoints_dsc$distanceRaw[i]<-sqrt(
    ((pointNorth-geo_surveyPoints_dsc$relativeNorthing[
      geo_surveyPoints_dsc$surveyPointID==leftPoint])^2)
    +((pointEast-geo_surveyPoints_dsc$relativeEasting[
      geo_surveyPoints_dsc$surveyPointID==leftPoint])^2))
}
lbPIN <- geo_surveyPoints_dsc$distanceRaw[
  geo_surveyPoints_dsc$surveyPointID=="DSC_LB_PIN"
]

### Transform raw distance to adjusted distance based on reference distance ####
geo_surveyPoints_dsc$distanceAdj <- NA
for(i in 1:(length(geo_surveyPoints_dsc$surveyPointID))){
  geo_surveyPoints_dsc$distanceAdj[i]<-geo_surveyPoints_dsc$distanceRaw[i]-lbPIN
}
geo_surveyPoints_dsc <- geo_surveyPoints_dsc[
  order(geo_surveyPoints_dsc$distanceAdj),
]

### Convert elevations of survey points in DSC transect to gauge height ####
gagElev <- geo_surveyPoints_gag$relativeHeight
gagMark <- as.numeric(regmatches(geo_surveyPoints_gag$surveyPointID,
                                 gregexpr("[-+]?\\d*\\.?\\d+", 
                                          geo_surveyPoints_gag$surveyPointID, 
                                          perl = TRUE))[[1]])
geo_surveyPoints_dsc$gaugeHeight <- NA
geo_surveyPoints_dsc$gaugeHeight<-round(
  geo_surveyPoints_dsc$relativeHeight - (gagElev - gagMark),
  digits=2
)

## Derive activation stages from rating curve priors ####
sdrc_priorParameters_WALK <- sdrc_WALK$sdrc_priorParameters[
  grepl(surveyDate,sdrc_WALK$sdrc_priorParameters$endDate),
]

## Plot control 1 rectangle ####

### Define boundaries for control 1 rectangle ####
c1Activate <- sdrc_priorParameters_WALK$priorActivationStage[
  sdrc_priorParameters_WALK$controlNumber==1
]
c1Left <- geo_surveyPoints_dsc$distanceAdj[
  geo_surveyPoints_dsc$surveyPointID == "DSC_XS10"
]# See header
c1Right <- geo_surveyPoints_dsc$distanceAdj[
  geo_surveyPoints_dsc$surveyPointID == "DSC_XS20"
]# See header
c1Max <- sdrc_priorParameters_WALK$priorActivationStage[
  sdrc_priorParameters_WALK$controlNumber==2
]# control 1 is not active in higher controls; see header

### Get x and y coordinates for control 1 rectangle ####
c1x <- c(c1Left,c1Right,c1Right,c1Left,c1Left)
c1y <- c(c1Activate,c1Activate,c1Max,c1Max,c1Activate)

## Plot control 2 rectangle ####

### Define boundaries for control 2 rectangle ####
c2Activate <- sdrc_priorParameters_WALK$priorActivationStage[
  sdrc_priorParameters_WALK$controlNumber==2
]
c2Left <- geo_surveyPoints_dsc$distanceAdj[
  geo_surveyPoints_dsc$surveyPointID == "DSC_LEW"
]# See header
c2Right <- geo_surveyPoints_dsc$distanceAdj[
  geo_surveyPoints_dsc$surveyPointID == "DSC_REW"
]# See header
c2Max <- max(
  geo_surveyPoints_dsc$gaugeHeight
)# control 2 is active in higher controls; see header

### Get x and y coordinates for control 2 rectangle ####
c2x <- c(c2Left,c2Right,c2Right,c2Left,c2Left)
c2y <- c(c2Activate,c2Activate,c2Max,c2Max,c2Activate)

## Plot control 3 rectangle ####

### Define boundaries for control 3 rectangle ####
c3Activate <- sdrc_priorParameters_WALK$priorActivationStage[
  sdrc_priorParameters_WALK$controlNumber==3
]
c3Left <- geo_surveyPoints_dsc$distanceAdj[
  geo_surveyPoints_dsc$surveyPointID == "DSC_LFH"
]# See header
c3Right <- geo_surveyPoints_dsc$distanceAdj[
  geo_surveyPoints_dsc$surveyPointID == "DSC_LEW"
]# See header
c3Max <- max(
  geo_surveyPoints_dsc$gaugeHeight
)# control 3 is highest control; see header

### Get x and y coordinates for control 3 rectangle ####
c3x <- c(c3Left,c3Right,c3Right,c3Left,c3Left)
c3y <- c(c3Activate,c3Activate,c3Max,c3Max,c3Activate)


## Compile final plot ####
xsplot <- geo_surveyPoints_dsc%>%
  ggplot2::ggplot(aes(x=distanceAdj,y=gaugeHeight))+
  ggplot2::geom_line(color="black")+
  ggplot2::geom_polygon(
    data=data.frame(x=c1x,y=c1y),
    aes(x=x,y=y,fill="C1"),
    linewidth=0.5,
    alpha=0.5)+
  ggplot2::geom_polygon(
    data=data.frame(x=c2x,y=c2y),
    aes(x=x,y=y,fill="C2"),
    linewidth=0.5,
    alpha=0.5)+
  ggplot2::geom_polygon(
    data=data.frame(x=c3x,y=c3y),
    aes(x=x,y=y,fill="C3"),
    linewidth=0.5,
    alpha=0.5)+
  ggplot2::geom_point(aes(fill="XS"), size=0, alpha=0)+
  ggplot2::scale_fill_manual(name=NULL, 
    values=c("XS"="black", 
             "C1"="#E69F00", 
             "C2"="#56B4E9", 
             "C3"="#009E73"),
    guide = guide_legend(
      override.aes = list(
        linetype = c(1, 1, 1, 1),
        shape = c(22, 22, 22, 95),
        size = c(4, 4, 4, 6),
        alpha = c(0.5, 0.5, 0.5, 1)),
      order = 1))+
  ggplot2::theme_bw(base_size = 12)+
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.margin = margin(4, 6, 4, 6),
    legend.box.spacing = unit(3, "pt"),
    legend.text = element_text(size = 8))+
  ggplot2::labs(x="Distance (m)",
                y="Stage (m)")+
  ggplot2::scale_y_continuous(n.breaks=6)+
  ggplot2::scale_x_continuous(n.breaks=6)

# Panel 2: Rating curve for WY 2024 ####

## Download RELEASED continuous discharge data ####
csd_WALK <- neonUtilities::loadByProduct(
  dpID = "DP4.00130.001",
  site = "WALK",
  startdate = queryStartDate,
  enddate = queryEndDate,
  release = queryRelease,
  package = "basic",
  check.size = F,
  token = Sys.getenv("NEON_PAT")
)
# Save downloaded data as RDS for quicker access in future runs
saveRDS(csd_WALK,paste0("./data/NEON.DP4.00130.001_WALK_",queryRelease,".rds"))

# If data is saved locally, read it in
# csd_WALK <- readRDS(paste0("./data/NEON.DP4.00130.001_WALK_",queryRelease,".rds"))

## Configure inputs and run posterior rating curve prediction model ####

### Configure Model Run ####
RunOptionsName <- "Config_RunOptions.txt"
stageQCurve::txt.out.run.opts(runType = "pred", 
                              RunOptionsPath = paste0(DIRPATH, 
                                                      BAMWS, 
                                                      RunOptionsName))
predMasterName <- "Config_Pred_Master.txt"
Config_Pred_Master <- readLines(paste0(DIRPATH, BAMWS, predMasterName))
Config_Pred_Master[1] <- gsub("[0-9] ","4 ",Config_Pred_Master[1])
Config_Pred_Master[2] <- gsub("'.*'","'Config_Pred_Prior.txt'",
                              Config_Pred_Master[2])
Config_Pred_Master[3] <- gsub("'.*'","'Config_Pred_RCMaxpost.txt'",
                              Config_Pred_Master[3])
Config_Pred_Master[4] <- gsub("'.*'","'Config_Pred_RCParamU.txt'",
                              Config_Pred_Master[4])
Config_Pred_Master[5] <- gsub("'.*'","'Config_Pred_RCTotalU.txt'",
                              Config_Pred_Master[4])
writeLines(Config_Pred_Master,paste0(DIRPATH, BAMWS, predMasterName))

### Configure and write out the model matrix from the rating curve priors ####
sdrc_controlInfo_WALK <- sdrc_WALK$sdrc_controlInfo[
  grepl(surveyDate,sdrc_WALK$sdrc_controlInfo$endDate),
]
numCtrls <- as.numeric(max(sdrc_controlInfo_WALK$controlNumber))
Config_ControlMatrix <- matrix(data=NA, nrow = numCtrls, ncol = numCtrls)
for(rw in 1:numCtrls){
  for(cl in 1:numCtrls){
    Config_ControlMatrix[rw,cl] <- as.numeric(
      sdrc_controlInfo_WALK$controlActivationState[
        sdrc_controlInfo_WALK$controlNumber == cl
        &sdrc_controlInfo_WALK$segmentNumber == rw
      ]
    )
  }
}
write.table(Config_ControlMatrix, 
            paste0(DIRPATH, BAMWS, "Config_ControlMatrix.txt"),
            row.names = F, col.names = F)

### Configure and write out the hydraulic control configurations ####
Config_Model <- matrix(data = NA, nrow = (4 + 12*numCtrls))
Config_Model[1] <- '"BaRatin"'
Config_Model[2:3] <- 1
Config_Model[4] <- 3 * numCtrls
for(j in 1:numCtrls){
  offset <- (j-1)*12
  
  #Divide by two and round to three places after the decimal
  kUnc <- format((as.numeric(sdrc_priorParameters_WALK$priorActivationStageUnc[
    sdrc_priorParameters_WALK$controlNumber == j
  ])/1.96), digits = 3)
  aUnc <- format((as.numeric(sdrc_priorParameters_WALK$priorCoefficientUnc[
    sdrc_priorParameters_WALK$controlNumber == j
  ])/1.96), digits = 3)
  cUnc <- format((as.numeric(sdrc_priorParameters_WALK$priorExponentUnc[
    sdrc_priorParameters_WALK$controlNumber == j
  ])/1.96), digits = 3)
  
  Config_Model[offset+5] <- paste0('"k', j, '"')
  Config_Model[offset+6] <- sdrc_priorParameters_WALK$priorActivationStage[
    sdrc_priorParameters_WALK$controlNumber == j
  ]
  Config_Model[offset+7] <- "'Gaussian'"
  Config_Model[offset+8] <- paste(
    sdrc_priorParameters_WALK$priorActivationStage[
      sdrc_priorParameters_WALK$controlNumber == j
    ],
    as.character(kUnc),
    sep = ",")
  Config_Model[offset+9] <- paste0('"a', j, '"')
  Config_Model[offset+10] <- sdrc_priorParameters_WALK$priorCoefficient[
    sdrc_priorParameters_WALK$controlNumber == j
  ]
  Config_Model[offset+11] <- "'Gaussian'"
  Config_Model[offset+12] <- paste(
    sdrc_priorParameters_WALK$priorCoefficient[
      sdrc_priorParameters_WALK$controlNumber == j
    ],
    as.character(aUnc),
    sep = ",")
  Config_Model[offset+13] <- paste0('"c', j, '"')
  Config_Model[offset+14] <- sdrc_priorParameters_WALK$priorExponent[
    sdrc_priorParameters_WALK$controlNumber == j
  ]
  Config_Model[offset+15] <- "'Gaussian'"
  Config_Model[offset+16] <- paste(
    sdrc_priorParameters_WALK$priorExponent[
      sdrc_priorParameters_WALK$controlNumber == j
    ],
    as.character(cUnc),
    sep = ",")
}
write.table(Config_Model, 
            paste0(DIRPATH, BAMWS, "Config_Model.txt"), 
            row.names = F, col.names = F, quote = F)

### Write out sequence of stage values between the min and max for WY 2024 ####
Hgrid <- seq(from = -0.1, 
             to = max(csd_WALK$csd_15_min$stageContinuous,na.rm=T),
             length.out = 181)
write.table(Hgrid,
            paste0(DIRPATH,BAMWS,"data/Hgrid.txt"),
            row.names = F,col.names = F)

### Write out the gaugings that make up the WY 2024 rating curve ####
sdrc_gaugeDischargeMeas_WALK <- sdrc_WALK$sdrc_gaugeDischargeMeas[
  sdrc_WALK$sdrc_gaugeDischargeMeas$curveID=="WALK.2024",
]
gagNam <- c('H','uH','bH','bHindx','Q','uQ','bQ','bQindx')
gaugings <- data.frame(
  matrix(data=NA, 
         ncol=length(gagNam), 
         nrow=length(sdrc_gaugeDischargeMeas_WALK$gaugeHeight))
)
names(gaugings) <- gagNam
gaugings$H <- sdrc_gaugeDischargeMeas_WALK$gaugeHeight
gaugings$uH <- 0.00
gaugings$bH <- 0.00
gaugings$bHindx <- 0.00
gaugings$Q <- sdrc_gaugeDischargeMeas_WALK$streamDischarge
gaugings$uQ <- sdrc_gaugeDischargeMeas_WALK$streamDischargeUnc
gaugings$bQ <- 0.00
gaugings$bQindx <- 0.00
write.table(gaugings,
            paste0(DIRPATH, BAMWS, "data/Gaugings.txt"),
            sep = "\t",
            row.names = F,
            quote = F)

### Write out the spaghettis that make up the WY2024 rating curve ####
sdrc_sampledParameters <- sdrc_WALK$sdrc_sampledParameters[
  sdrc_WALK$sdrc_sampledParameters$curveID=="WALK.2024",
]
stageQCurve::txt.out.spag.data(spagDataIn=sdrc_sampledParameters,
                               spagOutPath=paste0(DIRPATH, BAMWS, 
                                                  "Results_MCMC_Cooked.txt"))

### Write configuration and data files to the BaM folder for the water year ####
Config_Data <- readLines(paste0(DIRPATH, BAMWS, "Config_Data.txt"))
Config_Data[3] <- gsub("[0-9]{1,6}",nrow(gaugings),Config_Data[3])
writeLines(Config_Data, paste0(DIRPATH, BAMWS, "Config_Data.txt"))

### Run the BaM model ####
BaM_path <- gsub("/$","",DIRPATH)
if(!file.exists(BaM_path)){
  failureMessage <- "Path to BaM executable not found"
  stop(failureMessage)
}
setwd(BaM_path)
system2("BaM_MiniDMSL.exe")

## Read in model outputs and format for plotting ####

### Predicted Max Post Discharge ####
Qrc_Maxpost_spag <- read.table(paste0(DIRPATH,BAMWS, "Qrc_Maxpost.spag"),
                               header = F)
Qrc_Maxpost_spag$Hgrid <- Hgrid
totalUTop <- cbind.data.frame(Hgrid,Qrc_TotalU_env$Q_q2.5)
totalUBottom <- cbind.data.frame(Hgrid,Qrc_TotalU_env$Q_q97.5)
names(totalUTop) <- c("Hgrid","Q")
names(totalUBottom) <- c("Hgrid","Q")
totalUForPlotting <- rbind(totalUTop,totalUBottom[dim(totalUBottom)[1]:1,])

### Predicted Parametric Uncertainty ####
Qrc_ParamU_spag <- read.table(paste0(DIRPATH,BAMWS, "Qrc_ParamU.spag"),
                              header = F)
Qrc_ParamU_env <- read.table(paste0(DIRPATH,BAMWS, "Qrc_ParamU.env"),
                             header = T)
pramUForPlottingTop <- cbind.data.frame(Hgrid,Qrc_ParamU_env$Q_q2.5)
pramUForPlottingBottom <- cbind.data.frame(Hgrid,Qrc_ParamU_env$Q_q97.5)
names(pramUForPlottingTop) <- c("Hgrid","Q")
names(pramUForPlottingBottom) <- c("Hgrid","Q")
pramUForPlotting <- rbind(pramUForPlottingTop,
                          pramUForPlottingBottom[
                            dim(pramUForPlottingBottom)[1]:1,
                          ])

### Predicted Remnant Uncertainty ####
Qrc_TotalU_spag <- read.table(paste0(DIRPATH,BAMWS, "Qrc_TotalU.spag"),
                              header = F)
Qrc_TotalU_env <- read.table(paste0(DIRPATH,BAMWS, "Qrc_TotalU.env"),
                             header = T)

### Empirial Gauge and Discharge Pairs ####
gaugings <- read.table(paste0(DIRPATH,BAMWS, "data/Gaugings.txt"),
                       sep = "\t",header = T)
gaugings$Q <- as.numeric(gaugings$Q) #Convert to cms from lps
gaugings$uQ <- as.numeric(gaugings$uQ) #Convert to cms from lps

## Generate WALK WY 2024 posterior rating curve plot ####
rcdf <- data.frame(Hgrid = round(Hgrid,digits=3),
                   totalULower = Qrc_TotalU_env$Q_q2.5 * 1000,
                   totalUUpper = Qrc_TotalU_env$Q_q97.5 * 1000,
                   paraULower = Qrc_ParamU_env$Q_q2.5 * 1000,
                   paraUUpper = Qrc_ParamU_env$Q_q97.5 * 1000,
                   Q = Qrc_Maxpost_spag$V1 * 1000)
rcdf$totalULower[rcdf$totalULower<0.02] <- 0
rcdf$totalUUpper[rcdf$totalUUpper<0.02] <- 0
rcdf$paraULower[rcdf$paraULower<0.02] <- 0
rcdf$paraUUpper[rcdf$paraUUpper<0.02] <- 0
rcdf$Q[rcdf$H==0.005] <- 0.01
rcdf$Q[rcdf$H<0.005] <- NA
rcplot <- rcdf%>%
  ggplot2::ggplot()+
  ggplot2::geom_ribbon(aes(x = Hgrid,
                           ymin = totalULower, ymax = totalUUpper,
                           fill = "Total U"))+
  ggplot2::geom_ribbon(aes(x = Hgrid,
                           ymin = paraULower, ymax = paraUUpper,
                           fill = "Para U"))+
  ggplot2::geom_line(aes(x = Hgrid, y = Q), color="black")+
  ggplot2::geom_point(data=gaugings, aes(x=H, y=Q), color="black")+
  ggplot2::geom_point(aes(x = -999, y = -999, fill="RC"),
                      size=0, alpha=0)+
  ggplot2::geom_point(aes(x = -999, y = -999, fill="Field Q/H"),
                      size=0, alpha=0)+
  ggplot2::scale_fill_manual(name=NULL, 
    values=c("Total U"="#D55E00", 
             "Para U"="#E69F00", 
             "RC"="black", 
             "Field Q/H"="black"),
    breaks=c("Total U", "Para U", "RC", "Field Q/H"),
    guide = ggplot2::guide_legend(
      override.aes = list(
        linetype = c(1, 1, 1, 1),
        shape = c(22, 22, 95, 16),
        size = c(4, 4, 6, 2),
        alpha = c(1, 1, 1, 1)),
      order = 1))+
  ggplot2::scale_y_log10(n.breaks=6,
                         limits=c(0.01,450),
                         expand=c(0.01,NA))+
  ggplot2::scale_x_continuous(n.breaks=6,
                              limits=c(-0.05,NA),
                              expand=c(-0.05,NA))+
  ggplot2::labs(x = "Stage (m)",
                y = bquote("Discharge"~(L~s^-1))) +
  ggplot2::theme_bw(base_size = 12)+
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = c(0.57, 0.42),
    legend.justification = c(0, 1),
    legend.margin = margin(4, 6, 4, 6),
    legend.box.spacing = unit(0, "pt"),
    legend.text = ggplot2::element_text(size = 8))

# Panel 3: Continuous discharge timeseries ####
sdrc_gaugeDischargeMeas_WALK$collectDate <- as.POSIXct(paste(as.Date(gsub(
  "^.*\\.","",
  sdrc_gaugeDischargeMeas_WALK$gaugeEventID),
  format("%Y%m%d")),
  "16:00:00"),
  tz="UTC")
tsplot <- csd_WALK$csd_15_min%>%
  ggplot2::ggplot()+
  ggplot2::geom_ribbon(aes(x = endDateTime, 
                           ymin = dischargeLowerRemnUncert,
                           ymax = dischargeUpperRemnUncert,
                           fill = "Total U"))+
  ggplot2::geom_ribbon(aes(x = endDateTime, 
                           ymin = dischargeLowerParamUncert,
                           ymax = dischargeUpperParamUncert,
                           fill = "Para U"))+
  ggplot2::geom_line(aes(x=endDateTime,y=dischargeContinuous),
                     color="black")+
  ggplot2::geom_point(data=sdrc_gaugeDischargeMeas_WALK,
                      aes(x=collectDate,y=streamDischarge),
                      color="black")+
  ggplot2::geom_point(aes(x = as.POSIXct("2020-01-01"), 
                          y = 0, 
                          fill="Modeled Q"),
                      size=0, alpha=0)+
  ggplot2::geom_point(aes(x = as.POSIXct("2020-01-01"),
                          y = 0, fill="Field Q"),
                      size=0, alpha=0)+
  ggplot2::scale_fill_manual(name=NULL, 
    values=c("Total U"="#D55E00", 
             "Para U"="#E69F00", 
             "Modeled Q"="black", 
             "Field Q"="black"),
    breaks=c("Total U", "Para U", "Modeled Q", "Field Q"),
    guide = ggplot2::guide_legend(
      override.aes = list(
        linetype = c(1, 1, 1, 1),
        shape = c(22, 22, 95, 16),
        size = c(4, 4, 6, 2),
        alpha = c(1, 1, 1, 1)),
      order = 1))+
  ggplot2::labs(x = "Date",
                y = bquote("Discharge"~(L~s^-1))) +
  ggplot2::theme_bw(base_size = 12)+
  ggplot2::theme(
    axis.text.x = element_text(angle=10, hjust = 1),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = c(0.8, 0.98),
    legend.justification = c(0, 1),
    legend.margin = margin(4, 6, 4, 6),
    legend.box.spacing = unit(0, "pt"),
    legend.text = ggplot2::element_text(size = 8))+
  ggplot2::scale_y_continuous(limits=c(NA,310))+
  ggplot2::scale_x_datetime(date_breaks = "2 months",
                            date_labels = "%b %Y",
                            limits=c(as.POSIXct("2023-10-01",tz="UTC"),
                                     as.POSIXct("2024-10-01",tz="UTC")),
                            expand=c(0,0))

# Create multi-panel figure with cowplot ####
top_row <- cowplot::plot_grid(xsplot, rcplot, 
                              ncol = 2, 
                              labels = c("a)", "b)"),
                              label_size = 14)
bottom_row <- cowplot::plot_grid(tsplot, 
                                 ncol = 1, 
                                 labels = "c)",
                                 label_size = 14)
final_figure <- cowplot::plot_grid(top_row, bottom_row, 
                                   ncol = 1, 
                                   nrow = 2,
                                   rel_heights = c(1, 1))

# Write out the final plot ####
ggplot2::ggsave(plotSavePath,
                plot = final_figure,
                width = 6.5,height = 6.5,units = "in",
                dpi = 300)

# End ####