#'@title Microarray T-cell acitivation
#'
#'@description A data set with 44 replications of four genes involved in T-cells activation. For each gene there are ten repeated measurements recorded at 0, 2, 4, 6, 8, 18, 24, 32, 48, 72 hours after the treatement.
#' @format  A data frame with 44 rows and 40 columns.
#' \describe{
#' \item{rows}{independent and identical replications.}
#' \item{FYB.0}{expression levels for gene FYB at the time of treatment.}
#' \item{CD69.0}{expression levels for CD 69 at the time of treatment.}
#' \item{IL2RG.0}{expression levels for interleuking gene IL-2Rgamma at the time of treatment.}
#' \item{CDC2.0}{expression levels for cyclin A2 gene at the time of treatment.}
#' \item{FYB.2}{expression levels for gene FYB 2 hours after treatment.}
#' \item{CD69.2}{expression levels for CD 69 2 hours after treatment.}
#' \item{IL2RG.2}{expression levels for interleuking gene IL-2Rgamma 2 hours after treatment.}
#' \item{CDC2.2}{expression levels for cyclin A2 gene 2 hours after treatment.}
#' \item{FYB.4}{expression levels for gene FYB 4 hours after treatment.}
#' \item{CD69.4}{expression levels for CD 69 4 hours after treatment.}
#' \item{IL2RG.4}{expression levels for interleuking gene IL-2Rgamma 4 hours after treatment.}
#' \item{CDC2.4}{expression levels for cyclin A2 gene 4 hours after treatment.}
#' \item{FYB.6}{expression levels for gene FYB 6 hours after treatment.}
#' \item{CD69.6}{expression levels for CD 69 6 hours after treatment.}
#' \item{IL2RG.6}{expression levels for interleuking gene IL-2Rgamma 6 hours after treatment.}
#' \item{CDC2.6}{expression levels for cyclin A2 gene 6 hours after treatment.}
#' \item{FYB.8}{expression levels for gene FYB 8 hours after treatment.}
#' \item{CD69.8}{expression levels for CD 69 8 hours after treatment.}
#' \item{IL2RG.8}{expression levels for interleuking gene IL-2Rgamma 8 hours after treatmentt.}
#' \item{CDC2.8}{expression levels for cyclin A2 gene 8 hours after treatment.}
#' \item{FYB.18}{expression levels for gene FYB 18 hours after treatment.}
#' \item{CD69.18}{expression levels for CD 69 18 hours after treatment.}
#' \item{IL2RG.18}{expression levels for interleuking gene IL-2Rgamma 18 hours after treatment.}
#' \item{CDC2.18}{expression levels for cyclin A2 gene 18 hours after treatment.}
#' \item{FYB.24}{expression levels for gene FYB 24 hours after treatment.}
#' \item{CD69.24}{expression levels for CD 69 24 hours after treatment.}
#' \item{IL2RG.24}{expression levels for interleuking gene IL-2Rgamma 24 hours after treatment.}
#' \item{CDC2.24}{expression levels for cyclin A2 gene 24 hours after treatment.}
#' \item{FYB.32}{expression levels for gene FYB 32 hours after treatment.}
#' \item{CD69.32}{expression levels for CD 69 32 hours after treatment.}
#' \item{IL2RG.32}{expression levels for interleuking gene IL-2Rgamma 32 hours after treatment.}
#' \item{CDC2.32}{expression levels for cyclin A2 gene 32 hours after treatment.}
#' \item{FYB.48}{expression levels for gene FYB 48 hours after treatment.}
#' \item{CD69.48}{expression levels for CD 69 48 hours after treatment.}
#' \item{IL2RG.48}{expression levels for interleuking gene IL-2Rgamma 48 hours after treatment.}
#' \item{CDC2.48}{expression levels for cyclin A2 gene 48 hours after treatment.}
#' \item{FYB.72}{expression levels for gene FYB 72 hours after treatment.}
#' \item{CD69.72}{expression levels for CD 69 72 hours after treatment.}
#' \item{IL2RG.72}{expression levels for interleuking gene IL-2Rgamma 72 hours after treatment.}
#' \item{CDC2.72}{expression levels for cyclin A2 gene 72 hours after treatment.}
#' }
#' @source {R package longitudinal}<https://cran.r-project.org/web/packages/longitudinal/>
"Tcells"

#'@title Primary Biliary Cirrhosis for placebo and drug D-pencillamine (DPCA)
#'
#'@description A data set with three biomarkers recorded at seven time points for 40 patients in the placebo group and 42 in the DPCA group.
#' @format  A data frame with 82 rows and 21 columns.
#' \describe{
#' \item{rows}{first 40 rows are for the patients in placebo group and next 42 for those in DPCA group.}
#' \item{sr1}{serum bilirubin  measurements in mg/dl at the first visit.}
#' \item{al1}{serum albumin measurements in mg/dl at the first visit.}
#' \item{pr1}{prothrombin measurementrs in seconds at the first visit.}
#' \item{sr2}{serum bilirubin  measurements in mg/dl at the second visit.}
#' \item{al2}{serum albumin measurements in mg/dl at the second visit.}
#' \item{pr2}{prothrombin measurementrs in seconds at the second visit.}
#' \item{sr3}{serum bilirubin  measurements in mg/dl at the third visit.}
#' \item{al3}{serum albumin measurements in mg/dl at the third visit.}
#' \item{pr3}{prothrombin measurementrs in seconds at the third visit.}
#' \item{sr4}{serum bilirubin  measurements in mg/dl at the fourth visit.}
#' \item{al4}{serum albumin measurements in mg/dl at the fourth visit.}
#' \item{pr4}{prothrombin measurementrs in seconds at the fourth visit.}
#' \item{sr5}{serum bilirubin  measurements in mg/dl at the fifth visit.}
#' \item{al5}{serum albumin measurements in mg/dl at the fifth visit.}
#' \item{pr5}{prothrombin measurementrs in seconds at the fifth visit.}
#' \item{sr6}{serum bilirubin  measurements in mg/dl at the sixth visit.}
#' \item{al6}{serum albumin measurements in mg/dl at the sixth visit.}
#' \item{pr6}{prothrombin measurementrs in seconds at the sixth visit.}
#' \item{sr7}{serum bilirubin  measurements in mg/dl at the seventh visit.}
#' \item{al7}{serum albumin measurements in mg/dl at the seventh visit.}
#' \item{pr7}{prothrombin measurementrs in seconds at the seventh visit.}
#' }
#' @source {R package joineRML} <https://cran.r-project.org/web/packages/joineRML/index.html>
"pbc2"
