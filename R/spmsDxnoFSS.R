#' Diagnosis of secondary progressive multiple sclerosis without functional system scores and ambulation score
#' @description Diagnosis of conversion from relapsing-remitting multiple sclerosis (RRMS) to secondary progressive multiple sclerosis (SPMS), using the CORe definition without Functional System Scores (FSS) of Expanded Disability Status Scale (EDSS). Diagnosis is based on clinical visit records, each record including entries for patient code, visit date, EDSS score, and days since most recent relapse.
#' @references Lorscheider J, et al. Brain 2016; 139 (9): 2395-2405.
#' @references Brown JW, et al. JAMA 2019; 321 (2): 175-87.
#' @references Lizak N, et al. JAMA neurology 2020; 77 (11): 1398-407.
#' @param visits A dataframe consisting of 4 columns: ID, dateEDSS, EDSS, daysPostRelapse (days since most recent relapse).
#' @param minEDSS Minimum EDSS score to reach SPMS conversion.
#' @param tRelapse Minimum time in days from prior relapse to confirmation of EDSS progression.
#' @param tProgression SPMS confirmation period in days.
#' @param tRegression Confirmation period for EDSS improvement in days.
#' @param tRelProg Confirmation period (days) for rebaselining EDSS (after a relapse led to non-confirmed increase in EDSS).
#' @importFrom stats aggregate
#' @importFrom utils head tail
#' @importFrom dplyr %>% group_by mutate
#' @examples
#' data(SampleData)
#' output<-spmsDx_no_fss(SampleData)
#' @return A data frame.
#' @export
spmsDx_no_fss <- function(visits, minEDSS = 4, tRelapse = 30, tProgression = 3*30.25, tRegression = 9*30.25, tRelProg = 6*30.25) {

  if ( !is.element(minEDSS, c(seq(1, 6.0, 0.5))) | tRelapse < 0 | tProgression < 0 | tRegression < 0 | tRelProg < 0 ) {
    stop("Input requires minimum EDSS (minEDSS) 1.0 to 6.0 in 0.5-step increments (default = 4) and non-negative time intervals specified in days (tProgression, tRegression, tRelProg, tRelapse).")
  }

  if ( !all(with(visits, is.element(c("ID", "dateEDSS", "EDSS", "daysPostRelapse"), names(visits)))) ) {
    stop("Input requires a 'visits' data frame with columns 'ID', 'dateEDSS', 'EDSS', 'daysPostRelapse'.")
  }


  patients <- as.data.frame(unique(visits$ID))
  colnames(patients)[1] <- "ID"

  patients$SPMS <- 0
  patients$date.SPMS <- NA
  patients$EDSS.SPMS <- NA
  patients$interveningRelapse <- NA

  patients$EDSS.preDxFinal <- NA

  patients$date.obsFinal <- NA

  SPDXlist <- NULL

  visits$dateEDSS <- as.Date(visits$dateEDSS, "%Y-%m-%d")

  vis <- visits[which(!is.na(visits$EDSS)), ]
  vis <- vis[order(vis$ID, as.numeric(vis$dateEDSS)), ]

  rownames(vis) <- NULL

  tempDateObsFinal <- aggregate(vis$dateEDSS, by = list(vis$ID), FUN = max)
  patients$date.obsFinal <- tempDateObsFinal$x[match(patients$ID, tempDateObsFinal$Group.1)]
  rm(tempDateObsFinal)

  KurtzkeVars <- "EDSS"
  n <- nrow(patients)

  for (p in 1:n) {

    if ((p/n)==0.2||(p/n)==0.4|(p/n)==0.6|(p/n)==0.8|(p/n)==1.00) message("Processing patient #", p, "/", n, "  (", round(p*100/n, 2), "%)")

    pvisAll <- vis[which(vis$ID == patients$ID[p]), ]

    if (nrow(pvisAll) >= 3) {

      pvis <- pvisAll

      pvis<-pvis[order(pvis$ID, pvis$dateEDSS),]
      pvis$dateBlineVisit <- pvis$dateEDSS[1]
      pvis$daysPostBline <- as.numeric(pvis$dateEDSS - pvis$dateBlineVisit)
      pvis$daysPostBlineFinal <- as.numeric(max(pvis$daysPostBline))

      # set baseline values

      for (k in 1:length(KurtzkeVars)) {
        eval(parse(text = paste("pvis$b", KurtzkeVars[k], " <- pvis[which(colnames(pvis) == KurtzkeVars[k])][1, ]", sep = "")))
      }
      rm(k)

      if (nrow(pvis) >= 3) {

        daysPostPrior <- c(0, diff(pvis$daysPostBline))
        pvis$interveningRelapse <- (pvis$daysPostRelapse <= daysPostPrior)
        rm(daysPostPrior)

        pvis$dEDSSvPrior <- c(NA, diff(pvis$EDSS))

        pvisB <- pvis[which(pvis$dateEDSS >= pvis$dateBlineVisit), ]

        repeat {

          pvisB$dEDSSvBline <- pvisB$EDSS - pvisB$bEDSS

          # determine step EDSS progression & regression

          pvisB$progression <- 0
          pvisB$progression <- ifelse(pvisB$bEDSS == 0.0, ifelse(pvisB$dEDSSvBline >= 1.5 & pvisB$EDSS >= minEDSS, 1, 0), ifelse(pvisB$bEDSS < 6.0, ifelse(pvisB$dEDSSvBline >= 1 & pvisB$EDSS >= minEDSS, 1, 0), ifelse(pvisB$dEDSSvBline >= 0.5 & pvisB$EDSS >= minEDSS, 1, 0)))
          pvisB<- pvisB %>%
            group_by(ID) %>%
            mutate(flag_0 = cumsum(progression == 1)) %>%
            group_by(ID, flag_0) %>%
            mutate(conseq = ifelse(progression == 0, sum(progression == 0), 0)) %>%
            mutate(progression = ifelse(progression == 0 & conseq>1, 0, 1)) %>%
            ungroup()
          pvisB$regression <- 0
          pvisB$regression <- ifelse(pvisB$dEDSSvBline <= -0.5, 1, 0)

          pvisBP <- pvisB[!is.na(pvisB$dEDSSvBline), ]
          pvisBP[nrow(pvisBP) + 1, ] <- NA
          pvisBP$ID[nrow(pvisBP)] <- pvisBP$ID[1]
          pvisBP$progression[nrow(pvisBP)] <- 0
          pvisBP$regression[nrow(pvisBP)] <- 0
          pvisBP$rownumber <- c(1:nrow(pvisBP))
          pvisBP$finalProgVisnum <- NA
          pvisBP$finalRegVisnum <- NA

          # determine sustained EDSS progression

          rNoProg <- which(pvisBP$progression %in% 0)
          rNoRel <- which( (pvisBP$daysPostRelapse >= tRelapse) | is.na(pvisBP$daysPostRelapse) )
          for (v in 1:nrow(pvisBP)) {
            if (pvisBP$progression[v] == 0) {
              pvisBP$finalProgVisnum[v] <- NA
            } else {
              pvisBP$finalProgVisnum[v] <- max(rNoRel[ rNoRel < min( rNoProg[rNoProg > pvisBP$rownumber[v]])])
            }
          }
          rm(rNoProg, rNoRel, v)
          pvisBP$daysProgression <- pvisBP$daysPostBline[match(pvisBP$finalProgVisnum, pvisBP$rownumber)] - pvisBP$daysPostBline
          pvisBP$daysProgression[which(is.na(pvisBP$daysProgression))] <- 0

          # determine sustained EDSS regression

          rNoReg <- which(pvisBP$regression %in% 0)
          for (v in 1:nrow(pvisBP)) {
            if (pvisBP$regression[v] == 0) {
              pvisBP$finalRegVisnum[v] <- NA
            } else {
              pvisBP$finalRegVisnum[v] <- max(pvisBP$rownumber[pvisBP$rownumber < min(rNoReg[rNoReg > pvisBP$rownumber[v]])])
            }
          }
          rm(rNoReg, v)
          pvisBP$daysRegression <- pvisBP$daysPostBline[match(pvisBP$finalRegVisnum, pvisBP$rownumber)] - pvisBP$daysPostBline
          pvisBP$daysRegression[which(is.na(pvisBP$daysRegression))] <- 0

          iProg <- ifelse(pvisBP$daysProgression > tProgression, 1, 0)
          iReg <- ifelse(pvisBP$daysRegression > tRegression, 1, 0)
          pvisBP$progOrReg <- pmax(iProg, iReg)
          rm(iProg, iReg)


          if (sum(pvisBP$progOrReg == 1) == 0) { break }

          # identify relapses occurring after baseline and before the next progression or regression, and identify which are diagnostically significant

          pvisBPR <- pvisBP[1 : (which(pvisBP$progOrReg == 1)[1] - 1), ]
          pvisBPR <- pvisBPR[which(pvisBPR$dateEDSS != pvisBPR$dateBlineVisit[1]), ]
          pvisBPR$postRelapseSig <- NA
          if (nrow(pvisBPR) >= 1 & pvisBP$progOrReg[1] != 1) {
            pvisBPR$postRelapseSig <- FALSE
            for (v in 1:nrow(pvisBPR)) {
              if (is.na(pvisBPR$interveningRelapse[v])) {
                pvisBPR$postRelapseSig[v] <- FALSE & next
              }
              if (pvisBPR$interveningRelapse[v] == TRUE & pvisBPR$EDSS[v] > 3.0 & pvisBPR$dEDSSvBline[v] >= 0.5 & pvisBPR$dEDSSvPrior[v] >= 0.5 & pvisBPR$bEDSS[v] < 6.0) {
                pvisBPR$postRelapseSig[v] <- TRUE
              }
            }
            rm(v)
          }

          # analyze the diagnostically significant relapses

          if (sum(which(pvisBPR$postRelapseSig == TRUE)) >= 1) {

            if (max(which(pvisBPR$postRelapseSig == T)) == nrow(pvisBPR)) {
              blineVisit <- pvisBP[which(pvisBP$progOrReg == 1)[1], ]
              pvisB <- pvis[which(pvis$dateEDSS > blineVisit$dateEDSS), ]
              pvisB$dateBlineVisit <- blineVisit$dateEDSS
              pvisB$bEDSS <- blineVisit$EDSS
              rm(blineVisit)
            } else {
              pvisBPRN <- pvisBPR[(max(which(pvisBPR$postRelapseSig == T)) + 1) : nrow(pvisBPR), ]
              pvisBPRN[nrow(pvisBPRN) + 1, ] <- NA
              pvisBPRN$ID[nrow(pvisBPRN)] <- pvisBPRN$ID[1]
              pvisBPRN$blineSustained[nrow(pvisBPRN)] <- 0
              pvisBPRN$blineSustainedVisFinal <- NA
              for (v in 1:nrow(pvisBPRN)) {
                pvisBPRN$newbEDSS <- pvisBPRN$EDSS[v]
                pvisBPRN$newdEDSSvBline <- abs(pvisBPRN$EDSS - pvisBPRN$newbEDSS)
                pvisBPRN$blineSustained <- ifelse(pvisBPRN$newbEDSS == 0.0, ifelse(pvisBPRN$newdEDSSvBline >= 1.5, 0, 1), ifelse(pvisBPRN$bEDSS < 6.0, ifelse(pvisBPRN$newdEDSSvBline >= 1, 0, 1), ifelse(pvisBPRN$newdEDSSvBline >= 0.5, 0, 1)))
                pvisBPRN$blineSustained[nrow(pvisBPRN)] <- 0
                rBlineNotSustained <- which(pvisBPRN$blineSustained %in% 0)
                pvisBPRN$rownumber <- c(1:nrow(pvisBPRN))
                pvisBPRN$blineSustainedVisFinal[v] <- max(pvisBPRN$rownumber[pvisBPRN$rownumber < min(rBlineNotSustained[rBlineNotSustained > pvisBPRN$rownumber[v]])])
                rm(rBlineNotSustained)
              }
              rm(v)
              pvisBPRN$daysBlineSustained <- pvisBPRN$daysPostBline[match(pvisBPRN$blineSustainedVisFinal, pvisBPRN$rownumber)] - pvisBPRN$daysPostBline
              pvisBPRN$daysBlineSustained[which(is.na(pvisBPRN$daysBlineSustained))] <- 0
              pvisBPRN$blineValid <- ifelse(pvisBPRN$daysBlineSustained > tRelProg, 1, 0)
              if (sum(pvisBPRN$blineValid) == 0) {
                blineVisit <- pvisBP[which(pvisBP$progOrReg == 1)[1], ]
              } else {
                blineDate <- pvisBPRN$dateEDSS[min(which(pvisBPRN$blineValid == 1))]
                blineVisit <- pvisBP[which(pvisBP$dateEDSS == blineDate), ]
                rm(blineDate)
              }
              pvisB <- pvis[which(pvis$dateEDSS > blineVisit$dateEDSS), ]
              pvisB$dateBlineVisit <- blineVisit$dateEDSS
              pvisB$bEDSS <- blineVisit$EDSS
              rm(blineVisit, pvisBPRN)
            }
            if (nrow(pvisB) <= 1) { break }
            next
          }

          pvisBP <- pvisBP[which(pvisBP$progOrReg == 1)[1]:nrow(pvisBP), ]

          pvisBPall <- pvisBP

          # record whether an intervening relapse occured

          pvisBP1 <- pvisBP[1, ]
          pvisBP1$interveningRelapse[which(is.na(pvisBP1$interveningRelapse))] <- FALSE

          # confirmatory EDSS

          pvisBPall <- pvisBPall[which(pvisBPall$dateEDSS > pvisBP1$dateEDSS), ]
          pvisBPall$confirm <- ifelse((pvisBPall$daysPostBline - pvisBP1$daysPostBline) > tProgression, 1, 0)
          pvisBPall <- pvisBPall[1 : min(which(pvisBPall$confirm == 1)), ]
          pvisBP1$confirmedEDSS <- min(pvisBPall$EDSS)
          rm(pvisBPall)

          # determine whether (progression) + (confirmed over >= tProgression days) + (sustained throughout follow-up) + (no preceding relapses) + (EDSS >= minEDSS)

          if ( with( pvisBP1,
                     (progression == 1) &
                     (daysProgression >= tProgression) &
                     (daysPostBline + daysProgression >= daysPostBlineFinal) &
                     (interveningRelapse==FALSE) &
                     (EDSS >= minEDSS) &
                     (confirmedEDSS >= minEDSS))) {
            newrow <- subset(pvisBP1, select = c(ID, dateEDSS, EDSS, interveningRelapse, daysPostRelapse, dateBlineVisit, bEDSS, daysPostBline, daysProgression, dEDSSvBline))

            datePreDxFinal <- max(pvisAll$dateEDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)])
            EDSSPreDxFinal <- data.frame()
            EDSSPreDxFinal <- pvisAll$EDSS[which(pvisAll$dateEDSS == datePreDxFinal)]

            visitsRRMS <- length(pvisAll$dateEDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)])
            visitsSPMS <- length(pvisAll$dateEDSS[which(pvisAll$dateEDSS >= newrow$dateEDSS)])

            daysVisitIntervalRRMSMax <- ifelse(length(diff(pvisAll$dateEDSS[which(pvisAll$dateEDSS <= newrow$dateEDSS)]) > 0), max(diff(pvisAll$dateEDSS[which(pvisAll$dateEDSS <= newrow$dateEDSS)])), NA)
            daysVisitIntervalPreDxFinal <- as.numeric(difftime(newrow$dateEDSS, max(pvisAll$dateEDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)]), units = "days") )

            EDSSMinRRMS <- min(pvisAll$EDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)])
            EDSSMinSPMS <- min(pvisAll$EDSS[which(pvisAll$dateEDSS >= newrow$dateEDSS)])

            newrow <- cbind(newrow, datePreDxFinal, EDSSPreDxFinal, visitsRRMS, visitsSPMS, daysVisitIntervalRRMSMax, daysVisitIntervalPreDxFinal, EDSSMinRRMS, EDSSMinSPMS)
            rm(datePreDxFinal, EDSSPreDxFinal, visitsRRMS, visitsSPMS, daysVisitIntervalRRMSMax, daysVisitIntervalPreDxFinal, EDSSMinRRMS, EDSSMinSPMS)

            SPDXlist <- rbind(SPDXlist, newrow)
            rm(newrow)
            break
          } else {
            if (with(pvisBP1, (progression == 1) & (daysPostBline + daysProgression < daysPostBlineFinal) )) {
              pvisB <- pvisB[which(pvisB$daysPostBline > pvisBP1$daysPostBline + pvisBP1$daysProgression), ]
            } else {
              pvisB <- pvisB[which(pvisB$dateEDSS > pvisBP1$dateEDSS), ]
              if ( with( pvisBP1,
                         ( (progression == 1 & interveningRelapse == TRUE) | regression == 1 ) ) ) {
                blineVisit <- pvisBP1
                pvisB$dateBlineVisit <- blineVisit$dateEDSS
                pvisB$bEDSS <- blineVisit$EDSS
                rm(blineVisit)
              }
            }
          }
          if (nrow(pvisB) <= 1) { break }
        }
        pvisBP1<-NULL
      }
    }
  }
  rm(vis, KurtzkeVars, n, p, pvisAll, pvis)
  if (exists("pvisB")) { rm(pvisB) }
  if (exists("pvisBP")) { rm(pvisBP) }
  if (exists("pvisBPR")) { rm(pvisBPR) }
  if (exists("pvisBPRN")) { rm(pvisBPRN) }
  if (exists("pvisBP1")) { rm(pvisBP1) }

  # output: SPDX dataframe

  if (!is.null(SPDXlist)) {

    patients$date.SPMS <- as.Date(SPDXlist$dateEDSS[match(patients$ID, SPDXlist$ID)])
    patients$SPMS[which(!is.na(patients$date.SPMS))] <- 1
    patients$SPMS[which(is.na(patients$date.obsFinal))] <- NA

    patients$EDSS.SPMS <- SPDXlist$EDSS[match(patients$ID, SPDXlist$ID)]

    patients$EDSS.preDxFinal <- SPDXlist$EDSSPreDxFinal[match(patients$ID, SPDXlist$ID)]

  }

  SPDX <- patients
  SPDX<- subset(SPDX, select=-c(interveningRelapse, date.obsFinal))
  SPDX$dateEDSSFinal <- NULL
  rownames(SPDX) <- NULL

  attr(SPDX, 'minEDSS') <- minEDSS
  attr(SPDX, 'tRelapse') <- tRelapse
  attr(SPDX, 'tProgression') <- tProgression
  attr(SPDX, 'tRegression') <- tRegression
  attr(SPDX, 'tRelProg') <- tRelProg
  attr(SPDX, 'timestamp') <- paste(format(Sys.time(),"%Y%m%d%H%M"))

  return(SPDX)
}
globalVariables(c("ID", "dateEDSS", "EDSS", "FSpyr", "FSLeadConfirmed", "interveningRelapse", "daysPostRelapse", "dateBlineVisit", "bEDSS", "daysPostBline", "daysProgression", "dEDSSvBline", "date.obsFinal",
                  "progression", "flag_0", "conseq", "pvisBPR", "pvisBPRN", "%>%", "group_by", "mutate", "ungroup"))

