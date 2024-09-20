################# DEV SCRIPTS #####################
###################################################

## saving data for examples -----------------------
col.names<-c("seq","annotator","type","start","end","score","strand","phase","attribute")
aug_out<-data.table::fread("/Users/piyalkaru/Desktop/DDORF/Ann/HMA4/augustus_out/Alyr/Aug5_24/Alyrata_all_AL3G52820_blast_extr_scaff3_7.4_intron_fixed_mesq_fixed.fas.aug.txt",col.names = col.names)
annotation<-aug_out[grep("F1-14",aug_out$seq),c(3:5,7)]
annotation$seq<-"gene1"
save(annotation,file="/Users/piyalkaru/Desktop/DDORF/R/gSoup/data/annotation.rda",compress = "xz")

pi_and_theta<-readRDS("/Users/piyalkaru/Desktop/DDORF/Ann/HMA4/stats/diversity/alyr/Alyr_scaff3_full_gene_pi_gen_Aug26_24_v3.rds")
save(pi_and_theta,file="/Users/piyalkaru/Desktop/DDORF/R/gSoup/data/pi_and_theta.rda",compress="xz")

piA_and_S<-readRDS("/Users/piyalkaru/Desktop/DDORF/Ann/HMA4/stats/diversity/ahal/Alyr_scaff3_HMA4_full_gene_pi_AS_v5.rds")
save(piA_and_S,file="/Users/piyalkaru/Desktop/DDORF/R/gSoup/data/piA_and_S.rda",compress="xz")

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
devtools::check()
devtools::install()
