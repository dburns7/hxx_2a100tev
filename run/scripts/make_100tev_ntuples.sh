#!/bin/sh

source COMPILE

rm ../data/analysis/*.root 

# Using lumi in fb, xsec in fb^-1
# H -> llll ................................. 2.813706
# ZZ -> llll ................................ 71.637508096
# ZH, Z -> vv, H -> llll .................... 0.011097
# ZH, Z -> ll, H -> llvv .................... 0.3139155
# WH, W -> lv, H -> llll .................... 0.020300112

./bin/MakeHxxTree --maxevent=100 --sample=1 --lumi=20 --xsec=1 ../data/analysis/Higgs_hhxx_combined_1GeV_8TeV_ntuple.root ../data/preprocessed/Higgs_hhxx_combined_1GeV_8TeV.root

pushd ../data/analysis/
if [ -e "all.root" ] 
then 
  rm all.root
fi
hadd all.root *.root
popd
