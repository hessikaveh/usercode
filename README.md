## Tool for estimating the noise in ECAL crystals

#### Noise Plotter setup
Setup a working area for example in `CMSSW_9_0_2`. Any release `>=9XY`
should work just fine, contact me in case not. For older releases please refer to the branches
`cmssw_8x` and `cmssw_7x` for `CMSSW_8X` and `CMSSW_7X` respectively.
```bash
cmsrel CMSSW_9_0_2
cd CMSSW_9_0_2/src
cmsenv
git cms-init
git clone git@github.com:hessikaveh/usercode.git
git cms-merge-topic -u ferriff:ecal_calib_tools # works with some apparently harmless errors
cd usercode/
scram b -j10
```

Example: how to draw plots:
```bash
noise 1 -1 EcalADCToGeVConstant_2017_V1_Bon_mc EcalIntercalibConstants_V1_hlt  EcalLaserAlphas_2017_mc EcalPedestals_2017extrap_25fb_mc EcalLaserAPDPNRatios_weekly_hlt 1530432369

where the last number(1530432369) is the unix timestamp for 1st of July 2018, you can use the following website to change your desired time to unix timestamp:
https://www.epochconverter.com/

To draw Noise profile plots do the following at the end:
root -l plotProduct.C



```

#### Additional documentation
   * https://twiki.cern.ch/twiki/bin/view/CMS/DBDump
