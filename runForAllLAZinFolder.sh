#!/bin/bash

for i in $1/*.LAZ; do
    fullpath=$i
    file=$(basename $fullpath)
    filename="${file%.*}"
    ./SpatialCoherence $fullpath $2/$filename.tif $3 $4
done
