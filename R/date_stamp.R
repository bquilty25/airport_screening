#' Prepare a date and attribution statement
#'
#' @return A two element vector giving the date and attribution, and a link to
#' the referenced preprint.
#' @keywords internal
get_date_stamp <- function() {
  c(
    sprintf(
      paste0(
        "Last built at %s by B. Quilty, S. Clifford, S. Flasche, R. Eggo",
        "and other members of CMMID at LSHTM"
      ),
      format(
        Sys.time(),
        "%d %b %Y at %H:%M:%S"
      )
    ),
    paste0(
      "<br/><br/>Download the preprint of our analysis ",
      '<a href="https://github.com/cmmid/cmmid.github.io/raw/master',
      "/ncov/airport_screening_report/airport_screening_preprint_2020_01_28.",
      'pdf">here</a>.<br>Download the code for this app ',
      '<a href="https://github.com/bquilty25/airport_screening">on GitHub</a>',
      "<br/><br/>"
    )
  )
}
