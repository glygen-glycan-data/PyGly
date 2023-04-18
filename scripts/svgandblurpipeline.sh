#!/bin/sh
###### shell script for svg -> training pipeline

python2 randsvgs.py

#makes png, randomizes monosaccharide colors within a hard-coded range, creates bounding boxes from svg
for file in *.svg
do
	filename=${file%.*}
	outfile="${filename}.png"
	python2 svg2imagemap.py $file
  python3 svg2png.py $file $outfile
	python3 randomcolors.py $outfile
	python3 imagemap2boundingbox.py $outfile "${filename}.html" "${filename}.txt"
done

#makes directory for training images, deletes svg and html (you might want to keep the svg)
mkdir trainingimgs
mv *.png trainingimgs
mv *.txt trainingimgs
mv *.log trainingimgs

rm *.svg
rm *.html

