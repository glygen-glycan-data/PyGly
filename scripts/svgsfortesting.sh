#!/bin/sh
###### shell script for svg -> training pipeline


#give this the filename for the file with your trained glycans listed
python testingsvgs.py 1000trainingglycanslist


#loop through the created testing svgs, get bounding boxes and pngs
for file in *.svg
do
	filename=${file%.*}
	outfile="${filename}.png"
	python svg2imagemap.py $file
	python3 svg2png.py $file $outfile
	python3 imagemap2boundingbox.py $outfile "${filename}.html" "${filename}.txt"
done


#create directory for testing images; delete svgs and htmls (you might not want to delete the svgs but the htmls are not useful anymore)
mkdir testingimgs
mv *.png testingimgs
mv *.txt testingimgs
mv *.log testingimgs

rm *.svg
rm *.html

