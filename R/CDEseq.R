#' Identification of confirmed and sequential disability worsening and improvement events
#' @description Identify sequential disability worsening and improvement events confirmed over a specified time period, using roving baseline EDSS. The identification of events is based on clinical visit records, with each record including entries for patient code, visit date, EDSS score, and days since the most recent relapse.
#' If a baseline EDSS score is not provided, it is determined as the first EDSS score recorded in the dataset outside 30 days (the default) of a relapse.
#' Following a confirmed disability worsening or improvement event, the minimum EDSS score within the confirmation period, regardless of the recency of a relapse, becomes the new baseline EDSS score.
#' @references Sharmin, et al. European Journal of Neurology 2022;29(8):2321-2334.
#' @param Visits A data frame consisting of 6 columns: ID, dateEDSS, EDSS, daysPostRelapse (days since most recent relapse), bEDSS (baseline EDSS score), base.date (date of bEDSS).
#' @param mconf Confirmation period (days) for EDSS worsening or improvement.
#' @param tRelapse Minimum time in days since the most recent relapse to EDSS assessment.
#' @importFrom dplyr %>% mutate case_when row_number
#' @examples
#' data(SampleData)
#' output<-CDEseq(SampleData)
#' @return A data frame.
#' @export
CDEseq <- function(Visits, mconf=3*30.25, tRelapse=30) {

  if ( !all(with(Visits, is.element(c("ID", "dateEDSS", "EDSS", "daysPostRelapse"), names(Visits)))) ) {
    stop("Input requires a 'Visits' data frame with columns 'ID', 'dateEDSS', 'EDSS', 'daysPostRelapse'.")
  }

  Visits$ID <- as.character(Visits$ID)
  Visits$dateEDSS <- as.Date(Visits$dateEDSS, "%Y-%m-%d")
  Visits$daysPostRelapse <- abs(Visits$daysPostRelapse)
  Visits <- Visits[!is.na(Visits$EDSS),]

  Pats <- Visits$ID[!duplicated(Visits$ID)]
  n <- length(Pats)

  CDWseqlist <- NULL
  CDIseqlist <- NULL

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

    # Define disability event
    Visits.s <- Visits.sav

    repeat {

      if (all(is.na(Visits.s$bEDSS))) {break}

      Visits.s$dEDSS <- Visits.s$EDSS-Visits.s$bEDSS

      Visits.s$progression <- ifelse(Visits.s$bEDSS==0.0, ifelse(Visits.s$dEDSS>=1.5, 1, 0), ifelse(Visits.s$bEDSS<6.0, ifelse(Visits.s$dEDSS>=1, 1, 0), ifelse(Visits.s$dEDSS>=0.5, 1, 0)))
      Visits.s$regression <- ifelse(Visits.s$bEDSS<=1.5, ifelse(Visits.s$dEDSS<=-1.5, 1, 0), ifelse(Visits.s$bEDSS<=6.0, ifelse(Visits.s$dEDSS<=-1, 1, 0), ifelse(Visits.s$dEDSS<=-0.5, 1, 0)))

      if(Visits.s$progression[1]==1) {
        disability.event <- 1
      } else if (Visits.s$regression[1]==1){
        disability.event <- -1
      } else if(Visits.s$progression[1]==0 & Visits.s$regression[1]==0) {
        disability.event <- 0
      }

      if(disability.event>0){ ################# EDSS worsening  #################

        CL <- Visits.s[!is.na(Visits.s$dEDSS),]
        CL <- rbind(CL, CL[1,])
        CL$progression[length(CL$progression)] <- 0
        CL$rownr <- c(1:nrow(CL))
        zeros <- which(CL$progression %in% 0)
        norel <- which(CL$daysPostRelapse>tRelapse | is.na(CL$daysPostRelapse))
        for(j in 1:nrow(CL)) {
          if (CL$progression[j]==0) (CL$last.1[j] <- NA) else
            (CL$last.1[j] <- max(norel[norel<min(zeros[zeros>CL$rownr[j]])]))
        }
        rm(zeros, norel, j)
        CL$sust.prog <- CL$timepoint[match(CL$last.1, CL$rownr)]-CL$timepoint
        CL <- subset(CL, CL$sust.prog>=mconf)
        CP <- subset(CL[1,], select=-c(rownr, last.1))
        if (nrow(CL)==0) {break}
        CDWseqlist <- rbind(CDWseqlist, CP)
        Visits.s <- subset(Visits.s, Visits.s$timepoint>=CP$timepoint)

        # re-set EDSS baseline (minimum EDSS within the confirmation period)
        CL$prog.time <- CL$timepoint-CL$timepoint[1]
        min2 <- pmin(CL$EDSS[1], CL$EDSS[2])
        CL <- subset(CL, CL$prog.time<=mconf)
        min.prog <- min(CL$EDSS)
        Visits.s$bEDSS <- pmin(min2, min.prog, na.rm=T)
        Visits.s$base.date <- as.Date(Visits.s$dateEDSS[match(Visits.s$bEDSS, Visits.s$EDSS)])
        Visits.s <- subset(Visits.s, Visits.s$dateEDSS>=Visits.s$base.date)

      } else if (disability.event<0) { #################  EDSS improvement  #################

        CL <- Visits.s[!is.na(Visits.s$dEDSS),]
        CL <- rbind(CL, CL[1,])
        CL$regression[length(CL$regression)] <- 0
        CL$rownr <- c(1:nrow(CL))
        zeros <- which(CL$regression %in% 0)
        for(j in 1:nrow(CL)) {
          if (CL$regression[j]==0) (CL$last.1[j] <- NA) else
            (CL$last.1[j] <- max(CL$rownr[CL$rownr<min(zeros[zeros>CL$rownr[j]])]))
        }
        rm(zeros, j)
        CL$sust.reg <- CL$timepoint[match(CL$last.1, CL$rownr)]-CL$timepoint
        CL <- subset(CL, CL$sust.reg>=mconf)
        CP <- subset(CL[1,], select=-c(rownr, last.1))
        if (nrow(CL)==0) {break}
        CDIseqlist <- rbind(CDIseqlist, CP)
        Visits.s <- subset(Visits.s, Visits.s$timepoint>=CP$timepoint)

        # re-set EDSS baseline (minimum EDSS within the confirmation period)
        CL$reg.time <- CL$timepoint-CL$timepoint[1]
        min2 <- pmin(CL$EDSS[1], CL$EDSS[2])
        CL <- subset(CL, CL$reg.time<=mconf)
        min.reg <- min(CL$EDSS)
        Visits.s$bEDSS <- pmin(min2, min.reg, na.rm=T)
        Visits.s$base.date <- as.Date(Visits.s$dateEDSS[match(Visits.s$bEDSS, Visits.s$EDSS)])
        Visits.s <- subset(Visits.s, Visits.s$dateEDSS>=Visits.s$base.date)

      } else { #################  NEITHER worsening NOR improvement  #################
        Visits.s <- Visits.s[-1,]
      }
      if (nrow(Visits.s)<=1) {break}
    }
  }

  rm(n, i, Visits.sav)
  if (exists("CL")) { rm(CL) }
  if (exists("min2")) { rm(min2) }
  if (exists("min.prog")) { rm(min.prog) }
  if (exists("min.reg")) { rm(min.reg) }

  # Output: CDEseq dataframe
  if (is.null(CDWseqlist)) {
    CDWseqlist <- data.frame()
  } else {
    CDWseqlist$conc <- paste(CDWseqlist$ID, CDWseqlist$dateEDSS)
    CDWseqlist <- CDWseqlist[!duplicated(CDWseqlist$conc),]
    CDWseqlist <- subset(CDWseqlist, select=c(ID, bEDSS, base.date, EDSS, dateEDSS, dEDSS, timepoint, sust.prog))
    CDWseqlist <- CDWseqlist[which(as.numeric(as.character(CDWseqlist$timepoint))>0),]
    CDWseqlist <- subset(CDWseqlist, select= -timepoint)
    rownames(CDWseqlist) <- NULL

    attr(CDWseqlist, 'mconf') <- mconf
    attr(CDWseqlist, 'tRelapse') <- tRelapse
    attr(CDWseqlist, 'timestamp') <- paste(format(Sys.time(),"%Y%m%d%H%M"))
  }

  if (is.null(CDIseqlist)) {
    CDIseqlist <- data.frame()
  } else {
    CDIseqlist$conc <- paste(CDIseqlist$ID, CDIseqlist$dateEDSS)
    CDIseqlist <- CDIseqlist[!duplicated(CDIseqlist$conc),]
    CDIseqlist <- subset(CDIseqlist, select=c(ID, bEDSS, base.date, EDSS, dateEDSS, dEDSS, timepoint, sust.reg))
    CDIseqlist <- CDIseqlist[which(as.numeric(as.character(CDIseqlist$timepoint))>0),]
    CDIseqlist <- subset(CDIseqlist, select= -timepoint)
    rownames(CDIseqlist) <- NULL

    attr(CDIseqlist, 'mconf') <- mconf
    attr(CDIseqlist, 'tRelapse') <- tRelapse
    attr(CDIseqlist, 'timestamp') <- paste(format(Sys.time(),"%Y%m%d%H%M"))
  }

  CDEseqlist <- list(CDWseqlist, CDIseqlist)

  return(CDEseqlist)
}

globalVariables(c("ID", "dateEDSS", "EDSS", "daysPostRelapse", "bEDSS", "base.date", "timepoint", "dEDSS", "sust.prog", "sust.reg",
                  "progression", "regression", "rownr", "last.1"))

