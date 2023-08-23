## code to prepare `pathogen_parameters` dataset

# nolint begin
# bind pathogen parameters taken from the literature
pathogen_parameters <- do.call(
  rbind,
  list(
    data.frame(
      name = "nCoV-2019",
      # (Li et al. (2020) NEJM)
      mu_inc = 5.2,
      sigma_inc = 4.1,
      mu_inf = 9.1,
      sigma_inf = 14.7,
      prop.asy = 0.17
    ),
    data.frame(
      name = "SARS-like (2002)",
      mu_inc = 6.4,
      sigma_inc = 16.7,
      mu_inf = 3.8,
      sigma_inf = 6.0,
      prop.asy = 0.0
    ),
    data.frame(
      name = "Flu A/H1N1-like (2009)",
      mu_inc = 4.3,
      sigma_inc = 1.05,
      mu_inf = 9.3,
      sigma_inf = 0.7,
      prop.asy = 0.16 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586318/
    ),
    data.frame(
      name = "MERS-like (2012)",
      mu_inc = 5.5,
      sigma_inc = 6.25,
      # nolint begin
      # https://www.sciencedirect.com/science/article/
      # pii/S1473309913703049?via%3Dihub#sec1
      # nolint end
      mu_inf = 5.0, # https://www.nejm.org/doi/10.1056/NEJMoa1306742
      sigma_inf = 7.5,
      prop.asy = 0.21
      # nolint begin
      # https://doi.org/10.1016/j.tmaid.2018.12.003
      # citing https://www.who.int/csr/disease/coronavirus_infections/
      # risk-assessment-august-2018.pdf?ua=1
      # nolint end
    ),
    data.frame(
      name = "Custom",
      mu_inc = 5.0,
      sigma_inc = 5.0,
      mu_inf = 5.0,
      sigma_inf = 5.0,
      prop.asy = 0.5
    )
  )
)
# nolint end

usethis::use_data(pathogen_parameters, overwrite = TRUE)
