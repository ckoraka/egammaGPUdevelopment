**EGamma seeding efficiency studies**


This framework is used in order to check and compare the


**Install instructions**
```
cmsrel CMSSW_12_6_0_pre4
cd CMSSW_12_6_0_pre4/src
cmsenv
git clone git@github.com:ckoraka/egammaGPUdevelopment.git
scram b -j 8
```

**Run instructions**
```
cd egammaGPUdevelopment/egammaSeedEfficiency/test
cmsRun egSeedEff_cfg.py
```

**Note:** Remember to run 
```voms-proxy-init --voms cms``` 
to obtain grid certificate before trying to access any non-locan input .root file.
