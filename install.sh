#!/bin/bash

mkdir Outputs/

mkdir Files/

mkdir Programs/bin

cd Programs/

tar xvf muscle-5.1.0.tar.gz

cd muscle-5.1.0/src

make

cd ..

cd ..

mv muscle-5.1.0/src/Linux/muscle bin/muscle5

cd TeloPort-master/

make

cd ..

mv TeloPort-master/build/apps/telomereFinder bin/

mv TeloPort-master/build/apps/junctionFinder bin/

mv TeloPort-master/build/apps/sequenceQuality bin/

mv TeloPort-master/build/apps/wcdInterrogate bin/

cd wcdest-source

bash ./configure

touch doc/wcd.info doc/wcd.pdf

make

cd ..

mv wcdest-source/src/wcd bin/

curl https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-src.tar.gz -o ncbi-blast-2.14.0+-src.tar.gz

tar xvf ncbi-blast-2.14.0+-src.tar.gz

cd ncbi-blast-2.14.0+-src/c++

bash ./configure

cd ReleaseMT/build

make all_r
#build
cd .. 
#ReleaseMt
cd ..
#c++
cd ..
#ncbi-blast
cd ..
#Programs
mv ncbi-blast-2.14.0+-src/c++/ReleaseMT/bin/blastn bin/

curl ftp://emboss.open-bio.org:21/pub/EMBOSS/EMBOSS-6.6.0.tar.gz -o EMBOSS-6.6.0.tar.gz

tar xvf EMBOSS-6.6.0.tar.gz

cd EMBOSS-6.6.0


bash ./configure

make

cd ..

mv EMBOSS-6.6.0/emboss/cons bin/

rm -r EMBOSS-6.6.0

rm EMBOSS-6.6.0.tar.gz

rm ncbi-blast-2.14.0+-src.tar.gz

rm -r ncbi-blast-2.14.0+-src

rm -r muscle-5.1.0

rm muscle-5.1.0.tar.gz

rm -r TeloPort-master

rm -r wcdest-source
