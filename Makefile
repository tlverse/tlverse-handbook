.PHONY: all
all: style book code
crc: style pdf

book:
	sh ./_build.sh

style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

code:
	rm R/*.R
	R CMD BATCH purl.R
	rm purl.Rout .RData

pdf:
	Rscript -e "bookdown::clean_book(TRUE)"
	Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
