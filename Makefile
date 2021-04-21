.PHONY: all
all: style book code

book:
	bash ./_build.sh

pdf:
	Rscript -e "bookdown::clean_book(TRUE)"
	Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"

style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

code:
	rm -f R_code/*.R
	R CMD BATCH purl.R
	rm purl.Rout .RData
