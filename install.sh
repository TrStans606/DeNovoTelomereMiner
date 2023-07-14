#!/bin/bash

mkdir Programs/bin

tar xzf Programs/muscle-5.1.0.tar.gz

cd Programs/muscle-5.1.0/src

make

cd ..

cd ..

mv Programs/muscle-5.1.0/src/Linux/muscle Programs/bin/muscle5

cd Programs/TeloPort-master/

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

curl ftp://emboss.open-bio.org:21/pub/EMBOSS/EMBOSS-6.6.0.tar.gz -o EMBOSS-6.6.0.tar.gz

tar xvf EMBOSS-6.6.0.tar.gz

cd EMBOSS-6.6.0


bash ./configure

make


