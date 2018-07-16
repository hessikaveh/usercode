## Tools for dumping CMS database payloads

#### Package contents
   * a `conddb_dumper` executable that can be used to dump (mostly ECAL)
     conditions in txt files from a given object-specific tag, with several options
   * a general plugin for dumping several conditions at the same time/run
     from a given Global Tag, via cmsRun + configuration file.
   * a general plugin for dumping ECAL events for closer inspection
   * some utilities developed for the monitoring of the ECAL laser monitoring corrections

#### Dumper list
   * `conddb_dumper.cpp`: self-explicative via the `-h/--help` options. Currently supported objects: 
      * `BeamSpotObjects`
      * `ESEEIntercalibConstants`
      * `ESGain`
      * `ESIntercalibConstants`
      * `EcalADCToGeVConstant`
      * `EcalChannelStatus`
      * `EcalClusterLocalContCorrParameters`
      * `EcalGainRatios`
      * `EcalIntercalibConstants`
      * `EcalIntercalibConstantsMC`
      * `EcalLaserAlphas`
      * `EcalPedestals`
      * `EcalPulseShapes`
      * `EcalTPGLinearizationConst`
      * `EcalTimeCalibConstants`
      * `EcalTimeOffsetConstant`
      * `RunInfo`
   * `lava_db.cpp`: validate a tag of the ECAL monitoring corrections
   * `lava_text.cpp`: validate a set of ECAL monitoring corrections starting from a txt file with dumped-corrections
   * `lava_db_compare.cpp`: compare two tags of ECAL monitoring corrections, by
                            matching the IOV of one tag to the closest IOV of
                            the other tag
   * `lava_db_cond.cpp`: (superseded) for ECAL monitoring corrections, compare the content
                         of a tag to the content of another tag, by matching
                         the IOV of one tag to the closest IOV of the other tag
   * `lava_db_dumpId.cpp`: dump ECAL monitoring correction histories for a given set of DetId (or all of them)
   * `lava_db2txt.cpp`: convert the content of ECAL monitoring corrections to a txt file
   * `lava_db2root.cpp`: convert the content of ECAL monitoring corrections to a root file
   * `merge_dump.cpp`: merges different txt files of ECAL monitoring
                       corrections with overlapping IOV into one single
                       coherent txt file (N.B. does not need the CMSSW environment to work)

#### Dumper setup
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

```

#### Additional documentation
   * https://twiki.cern.ch/twiki/bin/view/CMS/DBDump
