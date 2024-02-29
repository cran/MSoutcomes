#' Identification of confirmed disability improvement events
#' @description Identify disability improvement events confirmed over a specified time period. The identification of events is based on clinical visit records, with each record including entries for patient code, visit date, EDSS score, and days since the most recent relapse.
#' If a baseline EDSS score is not provided, it is determined as the first EDSS score recorded in the dataset outside 30 days (the default) of a relapse.
#' Following a confirmed disability improvement event, the minimum EDSS score within the confirmation period, regardless of the recency of a relapse, becomes the new baseline EDSS score.
#' By default, only identify those improvement events that are sustained for the remainder of the follow-up.
#' @references Kalincik, et al. Brain 2015;138(11):3287-3298.
#' @param Visits A data frame consisting of 6 columns: ID, dateEDSS, EDSS, daysPostRelapse (days since most recent relapse), bEDSS (baseline EDSS score), base.date (date of bEDSS).
#' @param mconf Confirmation period (days) for EDSS improvement.
#' @param tRelapse Minimum time in days since the most recent relapse to EDSS assessment.
#' @param sustained If TRUE, the default, identifies only those EDSS improvement events sustained for the remaining recorded follow-up.
#' @examples
#' data(SampleData)
#' output<-CDI(SampleData)
#' @return A data frame.
#' @export
CDI <- function(Visits, mconf=3*30.25, tRelapse=30, sustained=TRUE) {

  if ( !all(with(Visits, is.element(c("ID", "dateEDSS", "EDSS", "daysPostRelapse"), names(Visits)))) ) {
    stop("Input requires a 'Visits' data frame with columns 'ID', 'dateEDSS', 'EDSS', 'daysPostRelapse'.")
  }

  Visits$ID <- as.character(Visits$ID)
  Visits$dateEDSS <- as.Date(Visits$dateEDSS, "%Y-%m-%d")
  Visits$daysPostRelapse <- abs(Visits$daysPostRelapse)
  Visits <- Visits[!is.na(Visits$EDSS),]

  Pats <- Visits$ID[!duplicated(Visits$ID)]
  n <- length(Pats)

  CDIlist <- NULL

  for (i in 1:n) {

    if ((i/n)==0.2|(i/n)==0.4|(i/n)==0.6|(i/n)==0.8|(i/n)==1.00) message("Processing patient #", i, "/", n, "  (", round(i*100/n, 2), "%)")

    Visits.sav <- subset(Visits, Visits$ID==Pats[i])

    # Define bEDSS
    if(all(is.na(Visits.sav$bEDSS))) {
      Visits.sn.base <- subset(Visits.sav, Visits.sav$daysPostRelapse>tRelapse | is.na(Visits.sav$daysPostRelapse))
      Visits.sn.base <- Visits.sn.base[order(Visits.sn.base$ID, as.numeric(Visits.sn.base$dateEDSS)), ]
      Visits.sav$bEDSS <- Visits.sn.base$EDSS[1]
      Visits.sav$base.date <- as.Date(Visits.sn.base$dateEDSS[1])
    } else {
      Visits.sav$base.date <- as.Date(Visits.sav$base.date, "%Y-%m-%d")
    }

    # Remove visits recorded prior to the date of baseline
    Visits.sav <- subset(Visits.sav, Visits.sav$dateEDSS>=Visits.sav$base.date)
    Visits.sav <- Visits.sav[order(Visits.sav$ID, as.numeric(Visits.sav$dateEDSS)), ]

    # Define timepoint
    Visits.sav$timepoint <- difftime(Visits.sav$dateEDSS, Visits.sav$base.date, units=(c="days"))
    Visits.sav$timepoint <- round(Visits.sav$timepoint, 0)

    # EDSS improvement
    Visits.s <- Visits.sav

    repeat {

      if (all(is.na(Visits.s$bEDSS))) {break}

      Visits.s$dEDSS <- Visits.s$EDSS-Visits.s$bEDSS

      # determine step EDSS improvement
      Visits.s$regression <- 0
      Visits.s$regression <- ifelse(Visits.s$bEDSS<=1.5, ifelse(Visits.s$dEDSS<=-1.5, 1, 0), ifelse(Visits.s$bEDSS<=6.0, ifelse(Visits.s$dEDSS<=-1, 1, 0), ifelse(Visits.s$dEDSS<=-0.5, 1, 0)))
      CL <- Visits.s[!is.na(Visits.s$dEDSS),]
      CL <- rbind(CL, CL[1,])

      # determine sustained EDSS improvement
      CL$rownr <- c(1:nrow(CL))
      zeros <- which(CL$regression %in% 0)
      for(j in 1:nrow(CL)) {
        if (CL$regression[j]==0) (CL$last.1[j] <- NA) else
          (CL$last.1[j] <- max(CL$rownr[CL$rownr<min(zeros[zeros>CL$rownr[j]])]))
      }
      rm(zeros, j)
      CL$sust.reg <- CL$timepoint[match(CL$last.1, CL$rownr)]-CL$timepoint
      if(sustained==TRUE){
        CL <- subset(CL, CL$sust.reg>=max(CL$timepoint)-CL$timepoint & CL$sust.reg>=mconf)
      } else {
        CL <- subset(CL, CL$sust.reg>=mconf)
      }
      CP <- subset(CL[1,], select=-c(rownr, last.1))
      if (nrow(CL)==0) {break}
      CDIlist <- rbind(CDIlist, CP)
      Visits.s <- subset(Visits.s, Visits.s$timepoint>=CP$timepoint)

      # re-set EDSS baseline (minimum EDSS within the confirmation period)
      CL$reg.time <- CL$timepoint-CL$timepoint[1]
      min2 <- pmin(CL$EDSS[1], CL$EDSS[2])
      CL <- subset(CL, CL$reg.time<=mconf)
      min.reg <- min(CL$EDSS)
      Visits.s$bEDSS <- pmin(min2, min.reg, na.rm=T)
      Visits.s$base.date <- as.Date(Visits.s$dateEDSS[match(Visits.s$bEDSS, Visits.s$EDSS)])
      Visits.s <- subset(Visits.s, Visits.s$dateEDSS>=Visits.s$base.date)
    }
  }
  rm(n, i, Visits.sav)
  if (exists("CL")) { rm(CL) }
  if (exists("min2")) { rm(min2) }
  if (exists("min.reg")) { rm(min.reg) }

  # Output: CDI dataframe
  CDIlist$conc <- paste(CDIlist$ID, CDIlist$dateEDSS)
  CDIlist <- CDIlist[!duplicated(CDIlist$conc),]
  CDIlist <- subset(CDIlist, select=c(ID, bEDSS, base.date, EDSS, dateEDSS, dEDSS, timepoint, sust.reg))
  CDIlist <- CDIlist[which(as.numeric(as.character(CDIlist$timepoint))>0),]
  CDIlist <- subset(CDIlist, select= -timepoint)
  rownames(CDIlist) <- NULL

  attr(CDIlist, 'mconf') <- mconf
  attr(CDIlist, 'tRelapse') <- tRelapse
  attr(CDIlist, 'sustained') <- sustained
  attr(CDIlist, 'timestamp') <- paste(format(Sys.time(),"%Y%m%d%H%M"))

  return(CDIlist)
}

globalVariables(c("ID", "dateEDSS", "EDSS", "daysPostRelapse", "bEDSS", "base.date", "timepoint", "dEDSS",
                  "regression", "sust.reg", "rownr", "last.1"))
