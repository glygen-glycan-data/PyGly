#!/bin/sh
###### shell script for svg -> training pipeline

mkdir trainingimgs
cd trainingimgs

#this will give 1000 random svgs
# randimgs.py now takes 3 command line arguments
# in order: number of total images you want, mode (svg/png), file with accessions to avoid
python2 ../randimgs.py 100 svg

#makes png, randomizes monosaccharide colors within a hard-coded range, creates bounding boxes from svg
for file in *.svg
do
	filename=${file%.*}
	outfile="${filename}.png"
	python2 ../svg2imagemap.py $file
    python3 ../svg2png.py $file $outfile
	python3 ../randomcolors.py $outfile
	python3 ../imagemap2boundingbox.py $outfile "${filename}.html" "${filename}.txt"
done

#makes directory for training images, deletes svg and html (you might want to keep the svg)
#mkdir trainingimgs
#mv *.png trainingimgs
#mv *.txt trainingimgs
#mv *.log trainingimgs

rm *.svg
rm *.html

