#!/bin/bash

(
cd Programs/TeloPort-master/ || exit
make clean
make
)
(
mkdir -p Programs/bin
)
mv Programs/TeloPort-master/build/apps/telomereFinder Programs/bin/

mv Programs/TeloPort-master/build/apps/junctionFinder Programs/bin/

mv Programs/TeloPort-master/build/apps/sequenceQuality Programs/bin/