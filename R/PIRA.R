#' Identification of progression independent of relapse activity (PIRA) events
#' @description Identify progression independent of relapse activity (PIRA) events confirmed over a specified time period. The identification of events is based on clinical visit records, each record including entries for patient code, visit date, EDSS score, and days since most recent relapse.
#' If a baseline EDSS score is not provided, it is determined as the first EDSS score recorded in the dataset outside 30 days (the default) of a relapse.
#' Following a confirmed PIRA event, the minimum EDSS score within the confirmation period, regardless of the recency of a relapse, becomes the new baseline EDSS score.
#' Following a relapse, the first EDSS score recorded in the dataset outside 30 days (the default) of a relapse, becomes the new baseline EDSS score.
#' By default, identify PIRA events that are sustained for the remainder of the follow-up.
#' @references Dzau, et al. Journal of Neurology, Neurosurgery & Psychiatry 2023;94(12):984-991.
#' @param Visits A data frame consisting of 6 columns: ID, dateEDSS, EDSS, daysPostRelapse (days since most recent relapse), bEDSS (baseline EDSS score), base.date (date of bEDSS).
#' @param mconf Confirmation period (days) for EDSS progression.
#' @param tRelapse Minimum time in days since the most recent relapse to EDSS assessment.
#' @param sustained If TRUE, the default, identifies only those PIRA events sustained for the remaining recorded follow-up.
#' @importFrom dplyr %>% mutate case_when row_number
#' @examples
#' data(SampleData)
#' output<-PIRA(SampleData)
#' @return A data frame.
#' @export
PIRA <- function(Visits, mconf=3*30.25, tRelapse=30, sustained=TRUE) {

  if ( !all(with(Visits, is.element(c("ID", "dateEDSS", "EDSS", "daysPostRelapse"), names(Visits)))) ) {
    stop("Input requires a 'Visits' data frame with columns 'ID', 'dateEDSS', 'EDSS', 'daysPostRelapse'.")
  }

  Visits$ID <- as.character(Visits$ID)
  Visits$dateEDSS <- as.Date(Visits$dateEDSS, "%Y-%m-%d")
  Visits$daysPostRelapse <- abs(Visits$daysPostRelapse)
  Visits <- Visits[!is.na(Visits$EDSS),]

  Pats <- Visits$ID[!duplicated(Visits$ID)]
  n <- length(Pats)

  PIRAlist <- NULL

  for (i in 1:n) {

    if ((i/n)==0.2|(i/n)==0.4|(i/n)==0.6|(i/n)==0.8|(i/n)==1.00) message("Processing patient #", i, "/", n, "  (", round(i*100/n, 2), "%)")

    Visits.sav <- subset(Visits, Visits$ID==Pats[i])

    # Define bEDSS
    if(all(is.na(Visits.sav$bEDSS))) {
      Visits.sn <- subset(Visits.sav, Visits.sav$daysPostRelapse>tRelapse | is.na(Visits.sav$daysPostRelapse))
      Visits.sn <- Visits.sn[order(Visits.sn$ID, as.numeric(Visits.sn$dateEDSS)), ]
      Visits.sav$bEDSS <- Visits.sn$EDSS[1]
      Visits.sav$base.date <- as.Date(Visits.sn$dateEDSS[1])
    } else {
      Visits.sav$base.date <- as.Date(Visits.sav$base.date, "%Y-%m-%d")
      Visits.sav <- subset(Visits.sav, Visits.sav$dateEDSS>=Visits.sav$base.date)
      Visits.sav<-rbind(Visits.sav[1,], Visits.sav)
      Visits.sav<-Visits.sav %>% dplyr::mutate(dateEDSS=dplyr::case_when(row_number()==1~base.date, TRUE ~ dateEDSS))
      Visits.sav<-Visits.sav %>% dplyr::mutate(EDSS=dplyr::case_when(row_number()==1~bEDSS, TRUE ~ EDSS))
      Visits.sav<-Visits.sav %>% dplyr::mutate(daysPostRelapse=dplyr::case_when(row_number()==1~NA, TRUE ~ daysPostRelapse))
    }

    Visits.sav <- Visits.sav[order(Visits.sav$ID, as.numeric(Visits.sav$dateEDSS)), ]

    # Define timepoint
    Visits.sav$timepoint <- difftime(Visits.sav$dateEDSS, Visits.sav$base.date, units=(c="days"))
    Visits.sav$timepoint <- round(Visits.sav$timepoint, 0)

    # EDSS progression
    Visits.s <- Visits.sav

    Visits.s$daysPostRelapse.timepoint <- as.numeric(Visits.s$timepoint)-Visits.s$daysPostRelapse
    Visits.s$daysPostRelapse.timepoint[is.na(Visits.s$daysPostRelapse.timepoint)] <- 0

    repeat {

      if (all(is.na(Visits.s$bEDSS))) {break}

      Visits.s$dEDSS <- Visits.s$EDSS-Visits.s$bEDSS

      # determine step EDSS progression
      Visits.s$progression <- 0
      Visits.s$progression <- ifelse(Visits.s$bEDSS==0.0, ifelse(Visits.s$dEDSS>=1.5, 1, 0), ifelse(Visits.s$bEDSS<6.0, ifelse(Visits.s$dEDSS>=1, 1, 0), ifelse(Visits.s$dEDSS>=0.5, 1, 0)))
      CL <- Visits.s[!is.na(Visits.s$dEDSS),]
      CL <- rbind(CL, CL[1,])

      # determine sustained EDSS progression
      CL$rownr <- c(1:nrow(CL))
      zeros <- which(CL$progression %in% 0)
      norel <- which(CL$daysPostRelapse>tRelapse| is.na(CL$daysPostRelapse))
      for(j in 1:nrow(CL)) {
        if (CL$progression[j]==0) (CL$last.1[j] <- NA) else
          (CL$last.1[j] <- max(norel[norel<min(zeros[zeros>CL$rownr[j]])]))
      }
      rm(zeros, norel, j)
      CL$sust.prog <- CL$timepoint[match(CL$last.1, CL$rownr)]-CL$timepoint
      if(sustained==TRUE){
        CL <- subset(CL, CL$sust.prog>=max(CL$timepoint)-CL$timepoint & CL$sust.prog>=mconf)
      } else {
        CL <- subset(CL, CL$sust.prog>=mconf)
      }
      CP <- subset(CL[1,], select=-last.1)
      if (nrow(CL)==0) {break}

      # re-baseline following a relapse
      if(as.numeric(CP$dateEDSS[1]-Visits.s$dateEDSS[1])>=CP$daysPostRelapse[1] & !is.na(CP$daysPostRelapse[1])){
        x <- min(Visits.s$daysPostRelapse.timepoint[Visits.s$daysPostRelapse.timepoint < CP$daysPostRelapse.timepoint])
        Visits.s <- subset(Visits.s, Visits.s$daysPostRelapse.timepoint>x)
        if (nrow(Visits.s)==0) {break}
        repeat{
          if(isTRUE((Visits.s$daysPostRelapse[1]<tRelapse & !is.na(Visits.s$daysPostRelapse[1])) | is.na(Visits.s$EDSS[1]))){
            Visits.s<-Visits.s[-1,]
          } else{break}
        }
        Visits.sn_rebaseline <- subset(Visits.s, Visits.s$daysPostRelapse>tRelapse | is.na(Visits.s$daysPostRelapse))
        Visits.s$bEDSS <- Visits.sn_rebaseline$EDSS[1]
        Visits.s$base.date <- Visits.sn_rebaseline$dateEDSS[1]
      } else{
        CP <- subset(CP, select=-c(rownr))
        PIRAlist <- rbind(PIRAlist, CP)
        Visits.s <- subset(Visits.s, Visits.s$timepoint>=CP$timepoint)

        # re-baseline following a PIRA event (minimum EDSS within the confirmation period)
        CL$prog.time <- CL$timepoint-CL$timepoint[1]
        min2 <- pmin(CL$EDSS[1], CL$EDSS[2])
        CL <- subset(CL, CL$prog.time<=mconf)
        min.prog <- min(CL$EDSS)
        Visits.s$bEDSS <- pmin(min2, min.prog, na.rm=T)
        Visits.s$base.date <- as.Date(Visits.s$dateEDSS[match(Visits.s$bEDSS, Visits.s$EDSS)])
        Visits.s <- subset(Visits.s, Visits.s$dateEDSS>=Visits.s$base.date)
      }
    }
  }
  rm(n, i, Visits.sav)
  if (exists("CL")) { rm(CL) }
  if (exists("min2")) { rm(min2) }
  if (exists("min.prog")) { rm(min.prog) }

  # Output: PIRA dataframe
  if (is.null(PIRAlist)) {
    PIRAlist <- data.frame()
  } else {
    PIRAlist$conc <- paste(PIRAlist$ID, PIRAlist$dateEDSS)
    PIRAlist <- PIRAlist[!duplicated(PIRAlist$conc),]
    PIRAlist <- subset(PIRAlist, select=c(ID, bEDSS, base.date, EDSS, dateEDSS, dEDSS, timepoint, sust.prog))
    PIRAlist <- PIRAlist[which(as.numeric(as.character(PIRAlist$timepoint))>0),]
    rownames(PIRAlist) <- NULL

    attr(PIRAlist, 'mconf') <- mconf
    attr(PIRAlist, 'tRelapse') <- tRelapse
    attr(PIRAlist, 'sustained') <- sustained
    attr(PIRAlist, 'timestamp') <- paste(format(Sys.time(),"%Y%m%d%H%M"))
  }

  return(PIRAlist)
}

globalVariables(c("ID", "dateEDSS", "EDSS", "daysPostRelapse", "base.date", "bEDSS", "timepoint", "sust.prog", "dEDSS",
                  "progression", "rownr", "last.1"))

