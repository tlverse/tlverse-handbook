TITLE=handbook

.PHONY: all
all: book

book:
	sh ./_build.sh

style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"
