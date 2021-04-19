#!/bin/sh

Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book('index.Rmd', quiet = FALSE)"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
