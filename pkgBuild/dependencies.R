
update_dependencies <- function(){
	devtools::use_package("rmarkdown", type="Imports")
	devtools::use_package("testthat", type="Imports")
	devtools::use_package("lintr", type="Imports")
	devtools::use_package("rticles", type="Imports")
	devtools::use_package("tidyverse", type="Suggests")
	
	devtools::use_package("stats", type="Imports")
	devtools::use_package("methods", type="Imports")
	devtools::use_package("rootSolve", type="Imports")
	devtools::use_package("deSolve", type="Imports")
}