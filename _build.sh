#!/bin/bash

echo "Building bookdown $BOOKDOWN_FORMAT"
if [ "$BOOKDOWN_FORMAT" == "SITE" ]
then
    BOOKDOWN_OUTPUT=bookdown_site
else
    BOOKDOWN_OUTPUT=pdf_book
fi

Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::$(BOOKDOWN_OUTPUT)')"
