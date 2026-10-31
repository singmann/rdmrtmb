#' Multi-source interference task (MSIT) data
#'
#' Responses and response times on a three-alternative MSIT, pooled over the two
#' sessions reported by Stevenson et al. (2024). On each trial three digits are
#' shown and participants report the digit that occurs only once. Two sources of
#' interference act on the target: flanker interference from the two distractor
#' digits and Simon interference from the position of the target in the
#' left-to-right array. The five trial types cross the two sources, following
#' the adapted MSIT design described in Appendix 1 of Stevenson et al. (2024).
#'
#' @format A data frame with 37,720 rows and 9 variables:
#' \describe{
#'   \item{subject}{Factor, participant identifier (79 participants).}
#'   \item{session}{Factor, session of testing, `"1"` or `"2"`.}
#'   \item{condition}{Factor, the five trial types: `"1_none"` (neither
#'     interference), `"2_simon"` (Simon interference only), `"3_flanker"`
#'     (flanker interference only), `"4_both_different"` (both, pointing to
#'     different incorrect responses) and `"5_both_same"` (both, pointing to
#'     the same incorrect response).}
#'   \item{R}{Factor, the observed response, `"1"`, `"2"` or `"3"`.}
#'   \item{RT}{Numeric, response time in seconds from stimulus onset.}
#'   \item{target_resp}{Character, the response supported by the target, i.e.,
#'     the correct response, `"1"` to `"3"`.}
#'   \item{flanker_resp}{Character, the response supported by the two identical
#'     flankers, `"1"` to `"3"`, or `"0"` if there were no flankers (conditions
#'     1 and 2).}
#'   \item{simon_resp}{Character, the response supported by the position of the
#'     target, `"1"` to `"3"`, or `"0"` if the position was consistent with the
#'     identity of the target (conditions 1 and 3).}
#'   \item{target_pos}{Character, the position of the target in the array,
#'     `"1"` (left) to `"3"` (right).}
#' }
#'
#' @source <doi:10.1007/s42113-023-00192-3>
#'
#' @references
#' Stevenson, N., Innes, R. J., Boag, R. J., Miletic, S., Isherwood, S. J. S.,
#' Trutti, A. C., Heathcote, A., & Forstmann, B. U. (2024). Joint modelling of
#' latent cognitive mechanisms shared across decision-making domains.
#' *Computational Brain & Behavior*, 7, 1--22.
#' \doi{10.1007/s42113-023-00192-3}
#'
#' Bush, G., Shin, L. M., Holmes, J., Rosen, B. R., & Vogt, B. A. (2003). The
#' Multi-Source Interference Task: Validation study with fMRI in individual
#' subjects. *Molecular Psychiatry*, 8, 60--70.
"msit"
