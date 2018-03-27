To build, run 'make'

If on svmit machine, a public build of RooUnfold should already be linked in Makefile at following path:
/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold/


In addition, please do the following:

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold/

I recommend above be put in your .bashrc file so it is not necessary every session

Notes: Using gcc 6.3.0, root 6.10.09, (comes for free with cmsenv in CMSSW_10_0_5)

To test:

./bin/buildAndTest.exe <inData> <inMC> <isPP>

If on Svmit, try:

./bin/buildAndTest.exe /data/cmcginn/RAA/HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20180218_11LumiPer_RawRAAHIST_NJet5_20180324.root /data/cmcginn/RAA/Dijet30And80And170And280PYTHIAHYDJET_RESPONSE_20180325.root 0 >& run.log