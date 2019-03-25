TITLE=handbook

.PHONY: all
all: book

book:
	Rscript -e "bookdown::clean_book(TRUE)"
	Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"

