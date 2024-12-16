#!/bin/bash

(
cd Programs/TeloPort-master/ || exit
make
)

mv Programs/TeloPort-master/build/apps/telomereFinder bin/

mv Programs/TeloPort-master/build/apps/junctionFinder bin/

mv Programs/TeloPort-master/build/apps/sequenceQuality bin/