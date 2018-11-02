#!/bin/bash
# Get all *location.txt files in the first argument and copy them to the second argument
# E.g. ./copygenelists.sh /cygdrive/q/labo\ biol\ mol/Oncoscan/OncoScan_results /cygdrive/q/labo\ biol\ mol/Oncoscan/genelists-wLocation

wordkdir=$1
destdir=$2
pushd `pwd` > /dev/null
cd "$1"
for f in */*location.txt; do
  d=${f#*/}
  df="$destdir/${d// /_}"
  if [ ! -e "$df" ]; then
    echo 'copying '$f
    cp "$f" "$df";
  fi
done
popd > /dev/null