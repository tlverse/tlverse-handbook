TITLE=handbook

.PHONY: all
all: book

book:
	Rscript -e "bookdown::render_book('.')"

