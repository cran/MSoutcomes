#' Identification of secondary progressive multiple sclerosis
#' @description Identify conversion from relapsing-remitting multiple sclerosis (RRMS) to secondary progressive multiple sclerosis (SPMS), using the CORe definition, including Functional System Scores (FSS) of Expanded Disability Status Scale (EDSS). The identification of SPMS is based on clinical visit records, each record including entries for patient code, visit date, EDSS score, FSS, ambulation score, and days since most recent relapse.
#' If a baseline EDSS score and corresponding FSS are not provided, these are determined as the first EDSS score and corresponding FSS recorded in the dataset, outside 30 days (the default) of a relapse.
#' Following a relapse, the first EDSS score recorded in the dataset outside 30 days (the default) of a relapse, becomes the new baseline EDSS score.
#' SPMS is sustained for the remainder of the follow-up, unless followed by two consecutive improvements in EDSS scores.
#' @references Lorscheider J, et al. Brain 2016; 139 (9): 2395-2405.
#' @param visits A data frame consisting of 22 columns: ID, dateEDSS, EDSS, FSpyr (pyramidal FSS), FScrbl (cerebellar FSS), FSbstem (brainstem FSS), FSsens (sensory FSS), FSsph (bowel bladder FSS), FSvis (visual FSS), FScereb (cerebral FSS), FSamb (ambulation score), dateBlineVisit, bEDSS (baseline EDSS), bFSpyr (baseline pyramidal FSS), bFScrbl (baseline cerebellar FSS), bFSbstem (baseline brainstem FSS), bFSsens (baseline sensory FSS), bFSsph (baseline bowel bladder FSS), bFSvis (baseline visual FSS), bFScereb (baseline cerebral FSS), bFSamb (baseline ambulation score), daysPostRelapse (days since most recent relapse).
#' @param minEDSS Minimum EDSS score required to reach SPMS conversion.
#' @param minFSpyr Minimum pyramidal FSS to reach SPMS conversion.
#' @param tRelapse Minimum time in days since the most recent relapse to EDSS assessment.
#' @param tProgression SPMS confirmation period in days.
#' @param tRegression Confirmation period for EDSS improvement in days.
#' @param tRelProg Confirmation period (days) for re-baselining EDSS (after a relapse led to non-confirmed increase in EDSS).
#' @importFrom stats aggregate
#' @importFrom dplyr %>% group_by mutate bind_cols bind_rows
#' @examples
#' data(SampleData)
#' output<-spmsDx(SampleData)
#' @return A data frame.
#' @export
spmsDx <- function(visits, minEDSS = 4, minFSpyr = 2, tRelapse = 30, tProgression = 3*30.25, tRegression = 9*30.25, tRelProg = 6*30.25) {

  if ( !is.element(minEDSS, c(seq(1, 6.0, 0.5))) | !is.element(minFSpyr, c(0, 2)) | tRelapse < 0 | tProgression < 0 | tRegression < 0 | tRelProg < 0 ) {
    stop("Input requires minimum EDSS (minEDSS) 1.0 to 6.0 in 0.5-step increments (default = 4); minimum pyramidal functional system score (minFSpyr) of 2; and non-negative time intervals specified in days (tProgression, tRegression, tRelProg, tRelapse).")
  }

  if ( !all(with(visits, is.element(c("ID", "dateEDSS", "EDSS", "FSpyr", "FScrbl", "FSbstem", "FSsens", "FSsph", "FSvis", "FScereb", "FSamb", "daysPostRelapse"), names(visits)))) ) {
    stop("Input requires a 'visits' data frame with columns 'ID', 'dateEDSS', 'EDSS', 'FSpyr', 'FScrbl', 'FSbstem', 'FSsens', 'FSsph', 'FSvis', 'FScereb', 'FSamb', 'daysPostRelapse'.")
  }


  patients <- as.data.frame(unique(visits$ID))
  colnames(patients)[1] <- "ID"

  patients$SPMS <- 0
  patients$date.SPMS <- NA
  patients$EDSS.SPMS <- NA
  patients$FSLeadConfirmed <- NA
  patients$interveningRelapse <- NA

  patients$EDSS.preDxFinal <- NA

  patients$date.obsFinal <- NA

  SPDXlist <- NULL

  visits$dateEDSS <- as.Date(visits$dateEDSS, "%Y-%m-%d")
  visits$daysPostRelapse <- abs(visits$daysPostRelapse)

  vis <- visits[which(!is.na(visits$EDSS)), ]
  vis <- vis[order(vis$ID, as.numeric(vis$dateEDSS)), ]
  vis$FSNA <- ifelse(is.na(vis$FSpyr) & is.na(vis$FScrbl) & is.na(vis$FSbstem) & is.na(vis$FSsens) & is.na(vis$FSsph) & is.na(vis$FSvis) & is.na(vis$FScereb) & is.na(vis$FSamb), 1, 0)

  rownames(vis) <- NULL

  tempDateObsFinal <- aggregate(vis$dateEDSS, by = list(vis$ID), FUN = max)
  patients$date.obsFinal <- tempDateObsFinal$x[match(patients$ID, tempDateObsFinal$Group.1)]
  rm(tempDateObsFinal)

  KurtzkeVars <- c("EDSS", "FSpyr", "FScrbl", "FSbstem", "FSsens", "FSsph", "FSvis", "FScereb", "FSamb")
  n <- nrow(patients)

  for (p in 1:n) {

    if ((p/n)==0.2|(p/n)==0.4|(p/n)==0.6|(p/n)==0.8|(p/n)==1.00) message("Processing patient #", p, "/", n, "  (", round(p*100/n, 2), "%)")

    pvisAll <- vis[which(vis$ID == patients$ID[p]), ]

    # Removing baseline visits with no FS data
    if ( (nrow(pvisAll) >= 3) & (length(which(pvisAll$FSNA == 0)) > 0)) {

      pvis <- pvisAll[min(which(pvisAll$FSNA == 0)) : nrow(pvisAll), ]

      pvis<-pvis[order(pvis$ID, pvis$dateEDSS),]

      # Define bEDSS
      if(all(is.na(pvis$bEDSS))) {
        Visits.sn <- subset(pvis, pvis$daysPostRelapse>tRelapse | is.na(pvis$daysPostRelapse))
        Visits.sn <- Visits.sn[order(Visits.sn$ID, as.numeric(Visits.sn$dateEDSS)), ]
        pvis$dateBlineVisit <- as.Date(Visits.sn$dateEDSS[1])
        for (k in 1:length(KurtzkeVars)) {
          eval(parse(text = paste("pvis$b", KurtzkeVars[k], " <- Visits.sn[which(colnames(Visits.sn) == KurtzkeVars[k])][1, ]", sep = "")))
        }
        rm(k)
      } else {
        pvis$dateBlineVisit <- as.Date(pvis$dateBlineVisit, "%Y-%m-%d")
      }

      pvis$daysPostBline <- as.numeric(pvis$dateEDSS - pvis$dateBlineVisit)
      pvis$daysPostBlineFinal <- as.numeric(max(pvis$daysPostBline))

      if (nrow(pvis) >= 3) {

        daysPostPrior <- c(0, diff(pvis$daysPostBline))
        pvis$interveningRelapse <- (pvis$daysPostRelapse <= daysPostPrior)
        rm(daysPostPrior)

        pvis$dEDSSvPrior <- c(NA, diff(pvis$EDSS))

        pvisB <- pvis[which(pvis$dateEDSS >= pvis$dateBlineVisit), ]

        repeat {

          pvisB$dEDSSvBline <- pvisB$EDSS - pvisB$bEDSS
          pvisB$dFSpyr <- pvisB$FSpyr - pvisB$bFSpyr
          pvisB$dFScrbl <- pvisB$FScrbl - pvisB$bFScrbl
          pvisB$dFSbstem <- pvisB$FSbstem - pvisB$bFSbstem
          pvisB$dFSsens <- pvisB$FSsens - pvisB$bFSsens
          pvisB$dFSsph <- pvisB$FSsph - pvisB$bFSsph
          pvisB$dFSvis <- pvisB$FSvis - pvisB$bFSvis
          pvisB$dFScereb <- pvisB$FScereb - pvisB$bFScereb
          pvisB$dFSamb <- pvisB$FSamb - pvisB$bFSamb

          # determine step EDSS progression & regression

          pvisB$progression <- 0
          pvisB$progression <- ifelse(pvisB$bEDSS == 0.0, ifelse(pvisB$dEDSSvBline >= 1.5 & pvisB$EDSS >= minEDSS, 1, 0), ifelse(pvisB$bEDSS < 6.0, ifelse(pvisB$dEDSSvBline >= 1 & pvisB$EDSS >= minEDSS, 1, 0), ifelse(pvisB$dEDSSvBline >= 0.5 & pvisB$EDSS >= minEDSS, 1, 0)))
          pvisB<- pvisB %>%
            dplyr::mutate(flag_0 = cumsum(progression == 1)) %>%
            dplyr::group_by(flag_0) %>%
            dplyr::mutate(conseq = ifelse(progression == 0, sum(progression == 0), 0)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(progression = ifelse(progression == 0 & conseq>1, 0, 1))

          pvisB$regression <- 0
          pvisB$regression <- ifelse(pvisB$dEDSSvBline <= -0.5, 1, 0)

          pvisBP <- pvisB[!is.na(pvisB$dEDSSvBline), ]
          pvisBP[nrow(pvisBP) + 1, ] <- NA
          pvisBP$ID[nrow(pvisBP)] <- pvisBP$ID[1]
          pvisBP$progression[nrow(pvisBP)] <- 0
          pvisBP$regression[nrow(pvisBP)] <- 0
          pvisBP$dFSpyr[nrow(pvisBP)] <- 0
          pvisBP$dFScrbl[nrow(pvisBP)] <- 0
          pvisBP$dFSbstem[nrow(pvisBP)] <- 0
          pvisBP$dFSsens[nrow(pvisBP)] <- 0
          pvisBP$dFSsph[nrow(pvisBP)] <- 0
          pvisBP$dFSvis[nrow(pvisBP)] <- 0
          pvisBP$dFScereb[nrow(pvisBP)] <- 0
          pvisBP$dFSamb[nrow(pvisBP)] <- 0
          pvisBP$rownumber <- c(1:nrow(pvisBP))
          pvisBP$finalProgVisnum <- NA
          pvisBP$finalRegVisnum <- NA

          # determine sustained EDSS progression

          rNoProg <- which(pvisBP$progression %in% 0)
          rNoRel <- which( (pvisBP$daysPostRelapse > tRelapse) | is.na(pvisBP$daysPostRelapse) )
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

          pvisBP$dFSpyr[which(is.na(pvisBP$dFSpyr))] <- 1

          if (sum(pvisBP$progOrReg == 1) == 0) { break }

          # identify relapses occuring after baseline and before the next progression or regression, and identify which are diagnostically significant

          pvisBPR <- pvisBP[1 : (which(pvisBP$progOrReg == 1)[1] - 1), ]
          pvisBPR <- pvisBPR[which(pvisBPR$dateEDSS != pvisBPR$dateBlineVisit[1]), ]
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
              pvisB$bFSpyr <- blineVisit$FSpyr
              pvisB$bFScrbl <- blineVisit$FScrbl
              pvisB$bFSbstem <- blineVisit$FSbstem
              pvisB$bFSsens <- blineVisit$FSsens
              pvisB$bFSsph <- blineVisit$FSsph
              pvisB$bFSvis <- blineVisit$FSvis
              pvisB$bFScereb <- blineVisit$FScereb
              pvisB$bFSamb <- blineVisit$FSamb
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
              if (sum(pvisBPRN$blineValid[which(pvisBPRN$FSNA == 0)]) == 0) {
                blineVisit <- pvisBP[which(pvisBP$progOrReg == 1)[1], ]
              } else {
                blineDate <- pvisBPRN$dateEDSS[min(which(pvisBPRN$blineValid == 1 & pvisBPRN$FSNA == 0))]
                blineVisit <- pvisBP[which(pvisBP$dateEDSS == blineDate), ]
                rm(blineDate)
              }
              pvisB <- pvis[which(pvis$dateEDSS > blineVisit$dateEDSS), ]
              pvisB$dateBlineVisit <- blineVisit$dateEDSS
              pvisB$bEDSS <- blineVisit$EDSS
              pvisB$bFSpyr <- blineVisit$FSpyr
              pvisB$bFScrbl <- blineVisit$FScrbl
              pvisB$bFSbstem <- blineVisit$FSbstem
              pvisB$bFSsens <- blineVisit$FSsens
              pvisB$bFSsph <- blineVisit$FSsph
              pvisB$bFSvis <- blineVisit$FSvis
              pvisB$bFScereb <- blineVisit$FScereb
              pvisB$bFSamb <- blineVisit$FSamb
              rm(blineVisit, pvisBPRN)
            }
            if (nrow(pvisB) <= 1) { break }
            next
          }

          # determine the leading FS

          pvisBP <- pvisBP[which(pvisBP$progOrReg == 1)[1]:nrow(pvisBP), ]

          tempFS <- pvisBP[ , c("dFSpyr", "dFScrbl", "dFSbstem", "dFSsens", "dFSsph", "dFSvis", "dFScereb", "dFSamb", "daysProgression", "daysPostBline")]
          tempFS1 <- tempFS[1, ]
          tempFS1[which(is.na(tempFS1))] <- 0

          cFSlead <- colnames(tempFS1[which(tempFS1[1:8] == max(tempFS1[1:8]))])

          tempFS <- tempFS[c(cFSlead, "daysProgression", "daysPostBline")]
          tempFS[nrow(tempFS) + 1, ] <- NA
          tempFS[nrow(tempFS), ] <- 0
          daysFSlead <- 0
          if (pvisBP$progression[1] > 0 & max(tempFS1[1, 1:8]) > 0) {
            for (i in 1:length(cFSlead)) {
              tempFS[is.na(tempFS[i]), i] <- 1
              daysFSlead[i] <- tempFS$daysPostBline[which(tempFS[i] == 0)[1] - 1] - tempFS$daysPostBline[1]
            }
            rm(i)
          }

          rm(tempFS, tempFS1, cFSlead)

          pvisBPall <- pvisBP

          # record whether the increase in the leading FS score was confirmed

          pvisBP$FSLeadConfirmed <- NA
          pvisBP$FSLeadConfirmed[1] <- length(which(daysFSlead >= tProgression)) > 0
          rm(daysFSlead)

          # record whether an intervening relapse occurred

          pvisBP1 <- pvisBP[1, ]
          pvisBP1$interveningRelapse[which(is.na(pvisBP1$interveningRelapse))] <- FALSE

          # confirmatory EDSS

          pvisBPall <- pvisBPall[which(pvisBPall$dateEDSS > pvisBP1$dateEDSS), ]
          pvisBPall$confirm <- ifelse((pvisBPall$daysPostBline - pvisBP1$daysPostBline) > tProgression, 1, 0)
          pvisBPall <- pvisBPall[1 : min(which(pvisBPall$confirm == 1)), ]
          pvisBP1$confirmedEDSS <- min(pvisBPall$EDSS)
          rm(pvisBPall)

          # determine whether (progression) + (confirmed over >= tProgression days) + (sustained throughout follow-up) + (no preceding relapses) + (EDSS >= minEDSS) + (pyramidal FS >= 2) + (confirmation in leading FS)

          pvisBP1$FSpyr[which(is.na(pvisBP1$FSpyr))] <- ifelse(pvisBP1$EDSS >= 5.5, 2, 0)
          pvisBP1$FSpyr[which(pvisBP1$FSpyr < 2 & pvisBP1$EDSS >= 5.5)] <- 2

          if ( with( pvisBP1,
                     (progression == 1) &
                     (daysProgression >= tProgression) &
                     (daysPostBline + daysProgression >= daysPostBlineFinal) &
                     (interveningRelapse==FALSE) &
                     (EDSS >= minEDSS) &
                     (confirmedEDSS >= minEDSS) &
                     (FSpyr >= minFSpyr) &
                     (FSLeadConfirmed == TRUE))) {
            newrow <- subset(pvisBP1, select = c(ID, dateEDSS, EDSS, FSpyr, FSLeadConfirmed, interveningRelapse, daysPostRelapse, dateBlineVisit, bEDSS, daysPostBline, daysProgression, dEDSSvBline))

            datePreDxFinal <- max(pvisAll$dateEDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)])
            EDSSPreDxFinal <- data.frame()
            EDSSPreDxFinal <- pvisAll$EDSS[which(pvisAll$dateEDSS == datePreDxFinal)]

            visitsRRMS <- length(pvisAll$dateEDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)])
            visitsSPMS <- length(pvisAll$dateEDSS[which(pvisAll$dateEDSS >= newrow$dateEDSS)])

            daysVisitIntervalRRMSMax <- ifelse(length(diff(pvisAll$dateEDSS[which(pvisAll$dateEDSS <= newrow$dateEDSS)]) > 0), max(diff(pvisAll$dateEDSS[which(pvisAll$dateEDSS <= newrow$dateEDSS)])), NA)
            daysVisitIntervalPreDxFinal <- as.numeric(difftime(newrow$dateEDSS, max(pvisAll$dateEDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)]), units = "days") )

            EDSSMinRRMS <- min(pvisAll$EDSS[which(pvisAll$dateEDSS < newrow$dateEDSS)])
            EDSSMinSPMS <- min(pvisAll$EDSS[which(pvisAll$dateEDSS >= newrow$dateEDSS)])

            newvars<-cbind(datePreDxFinal, EDSSPreDxFinal, visitsRRMS, visitsSPMS, daysVisitIntervalRRMSMax, daysVisitIntervalPreDxFinal, EDSSMinRRMS, EDSSMinSPMS)
            newvars<-data.frame(newvars)

            newrow <- dplyr::bind_cols(newrow, newvars)
            rm(newvars, datePreDxFinal, EDSSPreDxFinal, visitsRRMS, visitsSPMS, daysVisitIntervalRRMSMax, daysVisitIntervalPreDxFinal, EDSSMinRRMS, EDSSMinSPMS)

            SPDXlist<-data.frame(SPDXlist)
            SPDXlist <- dplyr::bind_rows(SPDXlist, newrow)
            rm(newrow)
            break
          } else {
            if (with(pvisBP1, (progression == 1) & (daysPostBline + daysProgression < daysPostBlineFinal) )) {
              pvisB <- pvis[which(pvis$daysPostBline > pvisBP1$daysPostBline + pvisBP1$daysProgression), ]
            } else {
              pvisB <- pvis[which(pvis$dateEDSS > pvisBP1$dateEDSS), ]
              if ( with( pvisBP1,
                         ( (progression == 1 & interveningRelapse == TRUE) | regression == 1 ) ) ) {
                if (pvisBP1$FSNA == 1) {
                  pvisBP1 <- pvisB[1, ]
                }
                repeat {
                  if (nrow(pvisB) <= 1) { break }
                  ifelse( (pvisBP1$FSNA == 1),
                          pvisB <- pvisB[2:nrow(pvisB), ],
                          break )
                  pvisBP1 <- pvisB[1, ]
                }
                blineVisit <- pvisBP1
                pvisB$dateBlineVisit <- blineVisit$dateEDSS
                pvisB$bEDSS <- blineVisit$EDSS
                pvisB$bFSpyr <- blineVisit$FSpyr
                pvisB$bFScrbl <- blineVisit$FScrbl
                pvisB$bFSbstem <- blineVisit$FSbstem
                pvisB$bFSsens <- blineVisit$FSsens
                pvisB$bFSsph <- blineVisit$FSsph
                pvisB$bFSvis <- blineVisit$FSvis
                pvisB$bFScereb <- blineVisit$FScereb
                pvisB$bFSamb <- blineVisit$FSamb
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
  SPDX<- subset(SPDX, select=-c(FSLeadConfirmed, interveningRelapse, date.obsFinal))
  SPDX$dateEDSSFinal <- NULL
  rownames(SPDX) <- NULL

  attr(SPDX, 'minEDSS') <- minEDSS
  attr(SPDX, 'minFSpyr') <- minFSpyr
  attr(SPDX, 'tRelapse') <- tRelapse
  attr(SPDX, 'tProgression') <- tProgression
  attr(SPDX, 'tRegression') <- tRegression
  attr(SPDX, 'tRelProg') <- tRelProg
  attr(SPDX, 'timestamp') <- paste(format(Sys.time(),"%Y%m%d%H%M"))

  return(SPDX)
}
globalVariables(c("ID", "dateEDSS", "EDSS", "FSpyr", "FSLeadConfirmed", "interveningRelapse", "daysPostRelapse", "dateBlineVisit", "bEDSS", "daysPostBline", "daysProgression", "dEDSSvBline", "date.obsFinal",
                  "progression", "flag_0", "conseq", "pvisBPR", "pvisBPRN", "%>%", "mutate", "ungroup", "bind_cols", "bind_rows"))
