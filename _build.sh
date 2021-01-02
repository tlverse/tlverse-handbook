#!/bin/sh

Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book('index.Rmd', quiet = TRUE)"
