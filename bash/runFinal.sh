#!/bin/bash

DATE=`date +%Y%m%d`

inNames=(output/firstRAAHist_nBayes8_DoSvd0_nSvd11_akCs3PU3PF_20180407_11.root output/firstRAAHist_nBayes8_DoSvd0_nSvd6_akCs10PU3PF_20180407_11.root output/firstRAAHist_nBayes8_DoSvd0_nSvd11_akCs4PU3PF_20180407_11.root output/firstRAAHist_nBayes8_DoSvd0_nSvd6_akCs8PU3PF_20180407_11.root output/firstRAAHist_nBayes8_DoSvd0_nSvd11_akCs6PU3PF_20180407_11.root)

tempStr=""
pos=0
for i in "${inNames[@]}"
do
    ./bin/makePlotValidation_FromTree.exe $i >& valid_$pos.log &

    pos=$((pos + 1))

    tempStr="$tempStr $i"
done