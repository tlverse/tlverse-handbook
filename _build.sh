#!/bin/bash

echo "Building bookdown $BOOKDOWN_FORMAT"
Rscript -e "bookdown::clean_book(TRUE)"

# run bookdown conditionally
if [ "$BOOKDOWN_FORMAT" == "SITE" ]
then
    Rscript -e "bookdown::render_book('index.Rmd', quiet = FALSE)"
else
    Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
fi
