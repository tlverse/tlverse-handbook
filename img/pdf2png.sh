#!/bin/bash
# requires imageMagick

for f in pdf/*.pdf
do
	echo "Converting PDF file: $f to ${f%.pdf}.png"
	convert -density 300 $f -quality 100 ${f%.pdf}.png
done

mv pdf/*.png png/
