TITLE=handbook

.PHONY: all
all: book

book:
	sh ./_build.sh

style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

code:
	rm R/*.R
	R CMD BATCH purl.R
	rm purl.Rout .RData
