
GlycoCT2ImageBundle.jar: pygly/GlycoCT2Image/GlycoCT2Image.java pygly/GlycoCT2FragImages/GlycoCT2FragImages.java
	ant -f build.xml

test: GlycoCT2ImageBundle.jar
	java -cp GlycoCT2ImageBundle.jar GlycoCT2Image format png orient BT scale 2.0 ../data/Named/*.txt
	java -cp GlycoCT2ImageBundle.jar GlycoCT2Image format jpg scale 4.0 notation cfglink display compact out test1.jpg ../data/Named/A2F1G1S1a.txt
	java -cp GlycoCT2ImageBundle.jar GlycoCT2Image notation cfgbw out test2.png ../data/Named/A2F1G1S1a.txt
	java -cp GlycoCT2ImageBundle.jar GlycoCT2Image display compact notation uoxfcol out test3.png ../data/Named/A2F1G1S1a.txt
	java -cp GlycoCT2ImageBundle.jar GlycoCT2Image display normal notation uoxf out test4.png - < ../data/Named/A2F1G1S1a.txt

test1: GlycoCT2ImageBundle.jar
	-rm -f A2F1G1S1a-[BY]*.png
	java -cp GlycoCT2ImageBundle.jar GlycoCT2FragImages ion BY ../data/Named/A2F1G1S1a.txt
	-rm -f A2F1G1S1b-[BY]*.png
	java -cp GlycoCT2ImageBundle.jar GlycoCT2FragImages ion BY display compact orient BT scale 1.0 ../data/Named/A2F1G1S1b.txt

clean: 
	- rm -f GlycoCT2Image*.jar *.png *.jpg

install: GlycoCT2ImageBundle.jar
	cp -f GlycoCT2ImageBundle.jar ../pygly
