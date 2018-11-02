#!/bin/bash
#Copy all 'genelist_fulllocation.txt' files from OncoScan_results and Oncoscan Reports to the switchdrive folder
#and call the python script compHRD.py (output to resultats-all.txt in the switchdrive folder).

#copy "full location" files
./copygenelists.sh /cygdrive/q/labo\ biol\ mol/Oncoscan/OncoScan_results /cygdrive/q/labo\ biol\ mol/Oncoscan/genelists-wLocation
./copygenelists.sh /cygdrive/q/labo\ biol\ mol/Oncoscan/Oncoscan\ Reports /cygdrive/q/labo\ biol\ mol/Oncoscan/genelists-wLocation

#compute HRD
export PYTHONPATH=Q:/1.BIOINFORMATIQUE/GitCentralRepo/CNV-tools/CNV-tools
python -m oncoscan_tools.hrd.compHRD /cygdrive/q/labo\ biol\ mol/Oncoscan/genelists-wLocation > /cygdrive/q/labo\ biol\ mol/Oncoscan/HRD-resultats_all.txt