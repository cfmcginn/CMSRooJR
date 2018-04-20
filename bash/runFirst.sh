#!/bin/bash

DATE=`date +%Y%m%d`


algoNames=(akCs3PU3PF akCs4PU3PF akCs6PU3PF akCs8PU3PF akCs10PU3PF)
pbpbName=/data/cmcginn/RAA/HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20180218_11LumiPer_RawRAAHIST_NJet5_20180324.root
ppName=/data/cmcginn/RAA/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20180219_5Lumi_Run262081to262328Fix_RawRAAHIST_NJet5_20180324.root
pbpbMCName=/data/cmcginn/RAA/Dijet30And80And170And280And370PYTHIAHYDJET_RESPONSE_20180406.root
ppMCName=/data/cmcginn/RAA/Dijet30And80And170And280PYTHIA_RESPONSE_20180325.root

for i in "${algoNames[@]}"
do
    ./bin/makeFirstRAAHist_FromTree.exe $pbpbName $ppName $pbpbMCName $ppMCName $i >& runFirst_$i\_$DATE.log &
    echo $i
done

