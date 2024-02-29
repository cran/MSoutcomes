#' Expanded disability status scale (EDSS) score and functional system score (FSS) recorded at each visit
#'
#' A long data frame containing 12 variables 'ID', 'dateEDSS', 'EDSS', 'FSpyr', 'FScrbl', 'FSbstem', 'FSsens', 'FSsph', 'FSvis', 'FScereb', 'FSamb', 'daysPostRelapse'.
#'
#' @format A long data frame with 798 rows and 12 variables:
#' \describe{
#'   \item{ID}{(character) patient ID}
#'   \item{dateEDSS}{(date YYYY-mm-dd) date of disability score}
#'   \item{EDSS}{(numeric) disability score (Expanded Disability Status Scale; EDSS)}
#'   \item{FSpyr}{(numeric) pyramidal functional system score}
#'   \item{FScrbl}{(numeric) cerebellar functional system score}
#'   \item{FSbstem}{(numeric) brainstem functional system score}
#'   \item{FSsens}{(numeric) sensory functional system score}
#'   \item{FSsph}{(numeric) bowel & bladder functional system score}
#'   \item{FSvis}{(numeric) visual functional system score}
#'   \item{FScereb}{(numeric) cerebral functional system score}
#'   \item{FSamb}{(numeric) ambulation functional system score}
#'   \item{daysPostRelapse}{(numeric) days since most recent relapse}
#' }
"SampleData"
