#!/bin/sh
###### shell script for svg -> training pipeline

#give this the list of glycans you trained on
python testingpngs.py orientationtrainingglycanslist

#create image map and bounding boxes
for file in *.
do
	filename=${file%.*}
	outfile="${filename}.png"
	python svg2imagemap.py $file
	python3 svg2png.py $file $outfile
	python3 imagemap2boundingbox.py $outfile "${filename}.html" "${filename}.txt"
done

mkdir testingimgs
mv *.png testingimgs
mv *.txt testingimgs
mv *.log testingimgs

rm *.svg
rm *.html

