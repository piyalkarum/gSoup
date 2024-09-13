################# DEV SCRIPTS #####################
###################################################


#### checking script performance
Rprof(tf <- "../performanceLog.log", memory.profiling = TRUE)
tt<-depthVsSample1(cov.len = 100,sam.len = 100)## the script
Rprof (NULL) ; print(summaryRprof(tf))




########### compiling for github and CRAN ###############
###################################################
## codecov
library(devtools)
use_testthat()
use_test()
use_coverage(type=c("codecov"))
use_github_action("test-coverage")
library(covr)
#codecov(token = "22ec927c-09d6-48ec-8d62-1e1a3d873f6e")
package_coverage()


######## vignette #######################
library(markdown)
usethis::use_vignette("rCNV")
usethis::use_github_actions()
usethis::use_readme_md()
#update version
usethis::use_version("major")
#add the license
usethis::use_agpl3_license()
## badges
usethis::use_github_actions_badge(name = "R-CMD-check", repo_spec = NULL)
usethis::use_lifecycle_badge("stable")
usethis::use_roxygen_md()
devtools::document()
usethis::use_logo("vignettes/logo.png")
#Github checks
usethis::use_github_action_check_standard()
#add news md
usethis::use_news_md()
## website
usethis::use_pkgdown()
pkgdown::build_site()
pkgdown::build_site_github_pages()
devtools::show_news()
# functions manual
devtools::build_manual()
# citation
usethis::use_citation()
### check with CRAN
devtools::check_win_release() # for Windows
devtools::check()
devtools::spell_check()
devtools::check_rhub()
devtools::check_win_devel()
####################################################################
################ RELEASE ########################
#################################################
devtools::check_rhub()
devtools::release(check = T)
# Resubmission
devtools::submit_cran()
# post submission
usethis::use_cran_badge()

####################################################################
# CPP implementation
Sys.setenv(R_MAKEVARS_USER = "~/.R/Makevars")
Sys.setenv(RCPP_VERBOSE = 1)



Rcpp::compileAttributes()

devtools::build()
devtools::install()
