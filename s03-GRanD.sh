#!/bin/sh
# allocate gauge (ID, lat, lon, uparea) on 1min MERIT river network

ONEMIN="./1min_v400"
rm -f 1min_river
ln -sf $ONEMIN 1min_river

#./src/allocate_gauge  GRDC_input.txt
#mv gauge_alloc.txt    GRDC_alloc.txt

./src/allocate_gauge  GRanD_input.csv
mv gauge_alloc.txt    GRanD_alloc.txt

