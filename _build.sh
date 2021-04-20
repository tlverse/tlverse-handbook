#!/bin/bash

echo "Building bookdown $BOOKDOWN_FORMAT"
case "$BOOKDOWN_FORMAT" in
   PDF) BOOKDOWN_OUTPUT="pdf_book" ;;
  SITE) BOOKDOWN_OUTPUT="bookdown_site" ;;
esac

Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::$(BOOKDOWN_OUTPUT)')"
