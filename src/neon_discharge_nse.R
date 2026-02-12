##############################################################################################
#' @title CALCULATE NASH-SUTCLIFFE EFFICIENCIES AND REGRESSION COEFFICIENTS FOR 
#' CORRECTED NEON DISCHARGE DATA

#' @author
#' Zachary Nickerson \email{nickerson@battelleecology.org} \cr

#' @description First, field and modeled discharge (15-min resolution) are 
#' downloaded from the NEON Data Portal. Each field discharge value is mapped to
#' the nearest modeled discharge timestamp. Plots are generated for each site
#' showing the linear fit of field vs. modeled discharge, then stacked into a 
#' single plot for output. Finally, an output table is generated that contains
#' the following summary statistics for each site: Nash-Sutcliffe Efficiency 
#' (NSE), slope of linear equation, coefficient of determination (R^2) of linear
#' equation, the percent of the period of record that contains NULL modeled
#' discharge, the channel slope of the site, and the watershed area of the site.

#' @details
#' Users should run this script within the .Rproj file in the root directory.
#' 
#' It is recommended that downloaders of NEON data have a NEON Data Portal
#' personal access token (PAT) set as an environment variable - 
#' Sys.setenv(NEON_PAT="YOUR PAT"). See instructions on obtaining a PAT at 
#' https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-tokens-tutorial 

#' @return NEON Portal downloads are saved as .rds to the /data subdirectory
#' @return Output plot (.png) and table (.csv) are saved to the /out subdirectory

# changelog and author contributions / copyrights
#   Zachary Nickerson (2026-01-14)
#     original creation
##############################################################################################

# Load libraries ####
library(neonUtilities)
library(scales)
library(ggplot2)
library(dplyr)
library(hydroGOF)
library(cowplot)

# Set options ####
options(stringsAsFactors=F)

# Set constants ####
queryStartDate <- "2021-10" # Beginning of WY 2022
queryEndDate <- "2023-09" # End of WY 2024
queryRelease <- "RELEASE-2026"

# Create subdirectories if they do not exist ####
if(!dir.exists("./data")){dir.create("./data")}
if(!dir.exists("./out")){dir.create("./out")}

# Download RELEASED modeled discharge data ####
csd <-neonUtilities::loadByProduct(dpID="DP4.00130.001",
                                   startdate=queryStartDate,
                                   enddate=queryEndDate,
                                   release=queryRelease,
                                   token=Sys.getenv("NEON_PAT"),
                                   check.size = F,
                                   package='basic')

# Save downloaded data as RDS for quicker access in future runs
saveRDS(csd,paste0("./data/NEON.DP4.00130.001_",queryRelease,".rds"))

# If data is saved locally, read it in
# csd <- readRDS(paste0("./data/NEON.DP4.00130.001_",queryRelease,".rds"))

# Download RELEASED field discharge data ####
dsc <-neonUtilities::loadByProduct(dpID="DP1.20048.001",
                                   startdate=queryStartDate,
                                   enddate=queryEndDate,
                                   release=queryRelease,
                                   token=Sys.getenv("NEON_PAT"),
                                   check.size = F,
                                   package='basic')

# Save downloaded data as RDS for quicker access in future runs
saveRDS(dsc,paste0("./data/NEON.DP1.20048.001_",queryRelease,".rds"))

# If data is saved locally, read it in
# dsc <- readRDS(paste0("./data/NEON.DP1.20048.001_",queryRelease,".rds"))

# Unpack the modeled and field discharge data frames ####
csd_15_min <- csd$csd_15_min
dsc_fieldData <- dsc$dsc_fieldData

# Subset out TOMB from field discharge records ####
dsc_fieldData <- dsc_fieldData[dsc_fieldData$siteID!="TOMB",]

# Associate the nearest modeled discharge value to each field record ####
dsc_fieldData$modeledQ <- NA
for(i in 1:nrow(dsc_fieldData)){
  # Subset modeled discharge by location
  if(dsc_fieldData$namedLocation[i]=="TOOK.AOS.discharge.inflow"){
    csd_15_min_loc <- csd_15_min[csd_15_min$siteID=="TOOK"
                                 &csd_15_min$horizontalPosition=="150",]
  }else{
    if(dsc_fieldData$namedLocation[i]=="TOOK.AOS.discharge.outflow"){
      csd_15_min_loc <- csd_15_min[csd_15_min$siteID=="TOOK"
                                   &csd_15_min$horizontalPosition=="160",]
    }else{
      csd_15_min_loc <- csd_15_min[csd_15_min$siteID==dsc_fieldData$siteID[i],]
    }
  }
  
  # Get nearest modeled discharge value
  dsc_fieldData$modeledQ[i] <- csd_15_min_loc$dischargeContinuous[
    csd_15_min_loc$startDateTime==csd_15_min_loc$startDateTime[
      which.min(abs(difftime(csd_15_min_loc$startDateTime,
                             dsc_fieldData$collectDate[i],
                             units = "secs")
                    )
                )
      ]
    ]
}

# Plot of field vs. continuous discharge relationship for each location ####

# Assign a plot label showing domain ID and location
dsc_fieldData$plotLabel <- dsc_fieldData$siteID
dsc_fieldData$plotLabel[
  dsc_fieldData$namedLocation=="TOOK.AOS.discharge.inflow"
] <- paste(dsc_fieldData$plotLabel[
  dsc_fieldData$namedLocation=="TOOK.AOS.discharge.inflow"
],"IN",sep = "-")
dsc_fieldData$plotLabel[
  dsc_fieldData$namedLocation=="TOOK.AOS.discharge.outflow"
] <- paste(dsc_fieldData$plotLabel[
  dsc_fieldData$namedLocation=="TOOK.AOS.discharge.outflow"
],"OT",sep = "-")


# Create a list to store individual plots
{
  plot_list <- list()
  dsc_fieldData$plotLabel <- gsub("^D[0-9]{2}\\-","",dsc_fieldData$plotLabel)
  labels <- unique(dsc_fieldData$plotLabel)
  labels <- labels[order(labels)]
  for (lbl in labels) {
    df_sub <- dsc_fieldData[dsc_fieldData$plotLabel == lbl, ]
    df_sub <- df_sub[is.finite(df_sub$finalDischarge)
                     &is.finite(df_sub$modeledQ),]
    if (nrow(df_sub) > 0) {
      # Calculate min/max for both axes for this plot only
      plot_x <- df_sub$finalDischarge
      plot_y <- df_sub$modeledQ
      axis_min <- min(c(plot_x, plot_y), na.rm = TRUE)
      axis_max <- max(c(plot_x, plot_y), na.rm = TRUE)
      axis_mid <- (axis_min + axis_max) / 2
      axis_breaks <- c(axis_min, axis_mid, axis_max)
  
      # Add a 5% buffer to both ends
      axis_range <- axis_max - axis_min
      buffer <- axis_range * 0.05
      axis_lim_min <- axis_min - buffer
      axis_lim_max <- axis_max + buffer
  
      p <- ggplot2::ggplot(df_sub,aes(x = finalDischarge, y = modeledQ)) +
        ggplot2::geom_point(size = 2, alpha = 0.7, color = "#0072B2") +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#D55E00", 
                             linetype = "solid", size = 1) +
        ggplot2::scale_x_continuous(breaks = axis_breaks, 
                                    limits = c(axis_lim_min, axis_lim_max), 
                                    labels = scales::label_scientific(digits = 0),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = axis_breaks, 
                                    limits = c(axis_lim_min, axis_lim_max), 
                                    labels = scales::label_scientific(digits = 0), 
                                    expand = c(0, 0)) +
        ggplot2::labs(
          x = NULL,
          y = NULL,
          title = lbl
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(face = "bold", 
                                             hjust = 0.5, 
                                             size = 10),
          axis.title = ggplot2::element_text(face = "bold",
                                             size = 8),
          axis.text.x = ggplot2::element_text(angle = 30,
                                              hjust = 1),
          aspect.ratio = 1
        ) +
        ggplot2::coord_fixed(xlim = c(axis_lim_min, axis_lim_max), 
                             ylim = c(axis_lim_min, axis_lim_max))
      plot_list[[lbl]] <- p
    }
  }
  
  # Combine all plots into a grid
  qRelPlot <- cowplot::plot_grid(plotlist = plot_list, ncol = 5)
  # Add more white space to the bottom for axis title
  qRelPlot <- cowplot::plot_grid(
    qRelPlot,
    NULL,
    rel_heights = c(1,0.05),
    ncol = 1
  )
  qRelPlot <- cowplot::plot_grid(
    NULL,
    qRelPlot,
    rel_widths = c(0.02,1),
    ncol = 2)+
    cowplot::draw_plot_label(c("Field~Discharge~(L~s^-1)", 
                               "Modeled~Discharge~(L~s^-1)"), 
                             x = c(0.25, -0.02), 
                             y = c(0.07, 0.25), 
                             angle=c(0,90), 
                             size=12, 
                             parse = TRUE)
}
ggplot2::ggsave(paste0("./out/discharge_fig_",queryRelease,".png"),
                plot = qRelPlot,
                width = 6.5,height = 8,units = "in",
                dpi = 300)

# Generate summary table with statistics for each location ####
summaryDF <- dsc_fieldData%>%
  dplyr::mutate(siteID_loc=gsub("^D[0-9]{2}\\-","",plotLabel))%>%
  dplyr::group_by(siteID_loc)%>%
  dplyr::summarise(
    NSE=abs(round(hydroGOF::NSE(modeledQ,finalDischarge,
                                na.rm=TRUE,fun=NULL,
                                epsilon.type = "none", epsilon.value = NA),
                  digits=2)),
    Slope=round(summary(lm(modeledQ~finalDischarge))$coefficients[2,1],
                digits=2),
    R2=round(summary(lm(modeledQ~finalDischarge))$r.squared,
             digits=2)
  )
summaryDF$perGap <- NA
for(i in 1:nrow(summaryDF)){
  if(!grepl("TOOK",summaryDF$siteID_loc[i])){
    curveIDSite <- summaryDF$siteID_loc[i]
  }else{
    if(grepl("\\-IN",summaryDF$siteID_loc[i])){
      curveIDSite <- "TKIN"
    }else{
      curveIDSite <- "TKOT"
    }
  }
  summaryDF$perGap[i] <- round(
    sum(is.na(csd_15_min$dischargeContinuous[
      grepl(curveIDSite,csd_15_min$curveID)]))/
      nrow(csd_15_min[
        grepl(curveIDSite,csd_15_min$curveID),])
  ,digits = 2)
}

# Add channel slope to each site ####

# Download RELEASED morphology map data - stream sites ####
geo <-neonUtilities::loadByProduct(dpID="DP4.00131.001",
                                   release=queryRelease,
                                   token=Sys.getenv("NEON_PAT"),
                                   check.size = F,
                                   package='basic')

# Save downloaded data as RDS for quicker access in future runs
saveRDS(geo,paste0("./data/NEON.DP4.00131.001_",queryRelease,".rds"))

# If data is saved locally, read it in
# geo <- readRDS(paste0("./data/NEON.DP4.00131.001_",queryRelease,".rds"))

# Assign the most recent slope value to each site ####
geo_surveySummary <- geo$geo_surveySummary
summaryDF$channelSlope <- NA
for(i in 1:nrow(summaryDF)){
  maxDate <- max(
    geo_surveySummary$surveyEndDate[
      geo_surveySummary$siteID==summaryDF$siteID_loc[i]
      &geo_surveySummary$surveyBoutTypeID=="geomorphology"],
    na.rm = T)
  if(maxDate!=-Inf){
    summaryDF$channelSlope[i] <- geo_surveySummary$thalwegSlope[
      geo_surveySummary$siteID==summaryDF$siteID_loc[i]
      &geo_surveySummary$surveyEndDate==maxDate]
  }
  rm(maxDate)
}

# Download RELEASED rating curve data - river and lake sites only ####
sdrc <-neonUtilities::loadByProduct(dpID="DP4.00133.001",
                                    site = c("BLWA","FLNT","TOOK"),
                                    release=queryRelease,
                                    token=Sys.getenv("NEON_PAT"),
                                    check.size = F,
                                    package='basic')

# Save downloaded data as RDS for quicker access in future runs
saveRDS(sdrc,paste0("./data/NEON.DP4.00133.001_",queryRelease,".rds"))

# If data is saved locally, read it in
# sdrc <- readRDS(paste0("./data/NEON.DP4.00133.001_",queryRelease,".rds"))

# Assign the most recent slope value to each site ####
sdrc_controlType <- sdrc$sdrc_controlType
sdrc_controlType$siteID_loc <- gsub("\\..*$","",sdrc_controlType$namedLocation)
sdrc_controlType$siteID_loc[
  grepl("inflow",sdrc_controlType$namedLocation)
] <- "TOOK-IN"
sdrc_controlType$siteID_loc[
  grepl("outflow",sdrc_controlType$namedLocation)
] <- "TOOK-OT"
for(i in 1:nrow(summaryDF)){
  if(summaryDF$siteID_loc[i]%in%c("BLWA","FLNT","TOOK-IN","TOOK-OT")){
    maxDate <- max(sdrc_controlType$endDate[
      grepl(summaryDF$siteID_loc[i],sdrc_controlType$siteID_loc)
    ])
    if(maxDate!=-Inf){
      summaryDF$channelSlope[i] <- unique(sdrc_controlType$channelSlope[
        sdrc_controlType$siteID_loc==summaryDF$siteID_loc[i]
        &sdrc_controlType$endDate==maxDate
        &!is.na(sdrc_controlType$channelSlope)])
    }
    rm(maxDate)
  }
}

# Add watershed area from the site metadata on the NEON Data Portal ####
# https://www.neonscience.org/field-sites/explore-field-sites
summaryDF$watershedArea <- NA
summaryDF$watershedArea[summaryDF$siteID_loc=='WALK'] <- 1
summaryDF$watershedArea[summaryDF$siteID_loc=='POSE'] <- 2
summaryDF$watershedArea[summaryDF$siteID_loc=='TECR'] <- 3
summaryDF$watershedArea[summaryDF$siteID_loc=='CUPE'] <- 4
summaryDF$watershedArea[summaryDF$siteID_loc=='COMO'] <- 4
summaryDF$watershedArea[summaryDF$siteID_loc=='MCRA'] <- 4
summaryDF$watershedArea[summaryDF$siteID_loc=='WLOU'] <- 5
summaryDF$watershedArea[summaryDF$siteID_loc=='MART'] <- 6
summaryDF$watershedArea[summaryDF$siteID_loc=='LECO'] <- 9
summaryDF$watershedArea[summaryDF$siteID_loc=='GUIL'] <- 10
summaryDF$watershedArea[summaryDF$siteID_loc=='BIGC'] <- 11
summaryDF$watershedArea[summaryDF$siteID_loc=='LEWI'] <- 12
summaryDF$watershedArea[summaryDF$siteID_loc=='HOPB'] <- 12
summaryDF$watershedArea[summaryDF$siteID_loc=='KING'] <- 13
summaryDF$watershedArea[summaryDF$siteID_loc=='MAYF'] <- 14
summaryDF$watershedArea[summaryDF$siteID_loc=='REDB'] <- 17
summaryDF$watershedArea[summaryDF$siteID_loc=='MCDI'] <- 23
summaryDF$watershedArea[summaryDF$siteID_loc=='CARI'] <- 31
summaryDF$watershedArea[summaryDF$siteID_loc=='BLDE'] <- 38
summaryDF$watershedArea[summaryDF$siteID_loc=='PRIN'] <- 49
summaryDF$watershedArea[summaryDF$siteID_loc=='OKSR'] <- 58
summaryDF$watershedArea[summaryDF$siteID_loc=='TOOK-OT'] <- 68
summaryDF$watershedArea[summaryDF$siteID_loc=='TOOK-IN'] <- 68
summaryDF$watershedArea[summaryDF$siteID_loc=='SYCA'] <- 280
summaryDF$watershedArea[summaryDF$siteID_loc=='BLUE'] <- 322
summaryDF$watershedArea[summaryDF$siteID_loc=='ARIK'] <- 2632
summaryDF$watershedArea[summaryDF$siteID_loc=='FLNT'] <- 14999
summaryDF$watershedArea[summaryDF$siteID_loc=='BLWA'] <- 16159

# Write out the table ####
write.csv(summaryDF,paste0("out/discharge_table_",queryRelease,".csv"),row.names = F)

# End