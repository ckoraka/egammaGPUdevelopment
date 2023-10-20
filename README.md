**EGamma seeding efficiency studies**


This framework is used in order to check and compare the doublets and triplets that can be constructed using gfs electron seed information


**Install instructions**
```
cmsrel CMSSW_13_0_0
cd CMSSW_13_0_0/src
cmsenv
git clone git@github.com:ckoraka/egammaGPUdevelopment.git
scram b -j 8
```

**Run instructions**
```
cd egammaGPUdevelopment/egammaSeedEfficiency/test
cmsRun hlt.py &> out.log
```

**Notes:** 

1. Remember to run ```voms-proxy-init --voms cms``` to obtain grid certificate before trying to access any non-locan input .root file.
2. To check the variables/collections stored in an edm root tree, under any cmssw release & after cmsenv do : 

```
edmDumpEventContent /store/relval/CMSSW_12_6_0_pre4/RelValZEE_14/GEN-SIM-DIGI-RAW/125X_mcRun3_2022_realistic_v4-v1/2580000/94a518d8-73e4-4cba-a97f-23f8c2c9834c.root
```
3. To access the e/g seeds collection from the hltgsfelectrons one has to re-run the hlt via ```hlt.py```. 
