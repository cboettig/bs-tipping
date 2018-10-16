
all:
		R -e 'lapply(list.files(pattern=".*.Rmd", recursive=TRUE), rmarkdown::render)'

