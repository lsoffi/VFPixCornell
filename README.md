#VFPixCornell

scram project -n CMSSW620-VFpix CMSSW CMSSW_6_2_0_SLHC23_patch1 
cd CMSSW_6_2_0_SLHC23_patch1/src/
git clone git@github.com:lsoffi/VFPixCornell.git VFPixCornell
cd TrkJetAnalyzer/python/
cmsRun ConfFile_cfg.py
