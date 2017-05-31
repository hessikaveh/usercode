#include "CondCore/CondDB/interface/ConnectionPool.h"
#include "CondCore/CondDB/interface/IOVProxy.h"
#include "CondCore/Utilities/interface/Utilities.h"
#include "CondFormats/Common/interface/TimeConversions.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "usercode/DBDump/interface/EcalLaserPlotter.h"
#include "usercode/DBDump/interface/EcalLaserDumper.h"

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <time.h>

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string>
#include <vector>

using namespace std;


namespace cond {
        class LaserValidation : public Utilities {
                public:
                        typedef EcalLaserAPDPNRatios A;
                        typedef EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap AMap;
                        typedef EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap::const_iterator AMapCit;
                        typedef EcalLaserAPDPNRatios::EcalLaserAPDPNpair AP;
                        typedef EcalLaserAPDPNRatios::EcalLaserTimeStamp AT;

                        LaserValidation();
                        ~LaserValidation();
                        int execute();
                        void initRoot();
                        int initRootToAppend();
                        void closeRoot();
                        void dbToRoot(const EcalLaserAPDPNRatios & corr);

                private: 
                        TTree *tree;
                        TFile *tfile;
                        int time[92];
                        float cor[75848];
                        int x[75848];
                        int y[75848];
                        int z[75848];
                        int eta[75848];
                        int phi[75848];
                        int detId[75848];
                        static std::string timeToString(time_t t);
                        std::vector<DetId> ecalDetIds_;
                        bool _verbose;
        };

}


std::string cond::LaserValidation::timeToString(time_t t)
{
        char buf[256];
        struct tm lt;
        localtime_r(&t, &lt);
        strftime(buf, sizeof(buf), "%F %R:%S", &lt);
        buf[sizeof(buf)-1] = 0;
        return string(buf);
}


cond::LaserValidation::LaserValidation():Utilities("cmscond_list_iov")
{
        addAuthenticationOptions();
        addOption<bool>("verbose","v","verbose");
        addOption<bool>("all","a","list all tags (default mode)");
        addOption<cond::Time_t>("beginTime","b","begin time (first since) (optional)");
        addOption<cond::Time_t>("endTime","e","end time (last till) (optional)");
        addOption<std::string>("tag","t","list info of the specified tag");
        addOption<std::string>("geom","g","geometry file (default: detid_geom.dat)");
        addOption<std::string>("output","o","output file (default: ecallaserplotter.root)");
        addOption<int>("niov","n","number of IOV to write");
        addOption<int>("skipiov","S","number of IOV to skip");
        addOption<int>("appendiov","N","number of IOV to append");
        addOption<int>("updateTree","U","update existing tree");

        ecalDetIds_.resize(EBDetId::MAX_HASH + 1 + EEDetId::kSizeForDenseIndexing);
        int idx = -1;
        for (int hi = EBDetId::MIN_HASH; hi <= EBDetId::MAX_HASH; ++hi ) {
                EBDetId ebId = EBDetId::unhashIndex(hi);
                if (ebId != EBDetId()) {
                        ecalDetIds_[hi] = ebId;
                }
        }
        for ( int hi = 0; hi < EEDetId::kSizeForDenseIndexing; ++hi ) {
                EEDetId eeId = EEDetId::unhashIndex(hi);
                if (eeId != EEDetId()) {
                        idx = EBDetId::MAX_HASH + 1 + hi;
                        ecalDetIds_[idx] = eeId;
                }
        }
        assert(ecalDetIds_.size() == 75848);
        assert(ecalDetIds_.size() == EBDetId::MAX_HASH + 1 + EEDetId::kSizeForDenseIndexing);



}


cond::LaserValidation::~LaserValidation()
{
}


void cond::LaserValidation::initRoot()
{
        std::stringstream title;
        std::stringstream fname;
        title << "Dump of Laser data from file";
        fname << "/tmp/DumpLaserDB";
        fname << ".root";
        std::cout << "Building tree " << title.str() << " on file " << fname.str()
                << std::endl;

        tfile = new TFile(fname.str().c_str(), "RECREATE", title.str().c_str());
        tree = new TTree("LDB", title.str().c_str());
        tree->Branch("time",          time,          "time[92]/I");
        tree->Branch("detId",         detId,         "detId[75848]/I");
        tree->Branch("eta",           eta,           "eta[75848]/I");
        tree->Branch("phi",           phi,           "phi[75848]/I");
        tree->Branch("x",             x,             "x[75848]/I");
        tree->Branch("y",             y,             "y[75848]/I");
        tree->Branch("z",             z,             "z[75848]/I");
        tree->Branch("cor",           cor,           "cor[75848]/F");
        std::cout << "Tree created" << std::endl;
}


int cond::LaserValidation::initRootToAppend()
{
        int result=0;
        std::stringstream title;
        std::stringstream fname;
        title << "Dump of Laser data from file";
        fname << "DumpLaserDB";
        fname << ".root";
        std::cout << "Building tree " << title.str() << " on file " << fname.str()
                << std::endl;

        tfile = new TFile(fname.str().c_str(), "UPDATE", title.str().c_str());
        // tree = new TTree("LDB", title.str().c_str());

        tree =(TTree*)gDirectory->Get("LDB");

        TBranch        *b_time;   //!
        TBranch        *b_detId;   //!
        TBranch        *b_eta;   //!
        TBranch        *b_phi;   //!
        TBranch        *b_x;   //!
        TBranch        *b_y;   //!
        TBranch        *b_z;   //!
        TBranch        *b_cor;   //!

        tree->SetBranchAddress("time", time, &b_time);
        tree->SetBranchAddress("detId", detId, &b_detId);
        tree->SetBranchAddress("eta", eta, &b_eta);
        tree->SetBranchAddress("phi", phi, &b_phi);
        tree->SetBranchAddress("x", x, &b_x);
        tree->SetBranchAddress("y", y, &b_y);
        tree->SetBranchAddress("z", z, &b_z);
        tree->SetBranchAddress("cor", cor, &b_cor);

        std::cout << "Tree re-loaded" << std::endl;

        result = (int) tree->GetEntriesFast();
        std::cout << ">> number of entries in tree: " << result<< std::endl;
        return result; 
}


void cond::LaserValidation::closeRoot()
{
        tree->Write();
        tfile->Close();
}


int cond::LaserValidation::execute()
{
        _verbose = hasOptionValue("verbose");

        // connects to the DB 
        std::string connect = getOptionValue<std::string>("connect" );
        cond::persistency::ConnectionPool connPool;
        if( hasOptionValue("authPath") ){
                connPool.setAuthenticationPath( getOptionValue<std::string>( "authPath") );
        }
        connPool.configure();
        cond::persistency::Session condDbSession = connPool.createSession( connect );
        condDbSession.transaction().start( true ); 

        // loads the tag with its IOV list
        std::string tag = getOptionValue<std::string>("tag");
        cond::persistency::IOVProxy iov = condDbSession.readIov(tag, true);

        // loads the geometry file 
        std::string geom = hasOptionValue("geom") ? getOptionValue<std::string>("geom") : "detid_geom.dat";
        std::string output = hasOptionValue("output") ? getOptionValue<std::string>("output") : "ecallaserplotter.root";


        int nIOVTree = 0;
        bool updateTree = false;
        int upd = 0;
        if (hasOptionValue("updateTree")) upd = getOptionValue<int>("updateTree");
        if(upd == 1) updateTree = true;
        /* start the root Tree */
        if(!updateTree) {
                initRoot();
        } else {
                nIOVTree = initRootToAppend();
                if(_verbose) std::cout << "the tree already contains "<<nIOVTree<<" entries"<<std::endl; 
        }

        cond::Time_t since = std::numeric_limits<cond::Time_t>::min();
        if( hasOptionValue("beginTime" )) since = getOptionValue<cond::Time_t>("beginTime");
        cond::Time_t till = std::numeric_limits<cond::Time_t>::max();
        if( hasOptionValue("endTime" )) till = getOptionValue<cond::Time_t>("endTime");

        since = std::max((cond::Time_t)2, since ); 
        till  = std::min(till,  iov.getLast().since);

        if(_verbose) std::cout << "tag " << tag << std::endl ;
        if(_verbose) std::cout << "since/till  " << since << "/"<< till<< std::endl ;

        int niov = -1;
        if (hasOptionValue("niov")) niov = getOptionValue<int>("niov");

        int skipiov = -1;
        if (hasOptionValue("skipiov")) skipiov = getOptionValue<int>("skipiov");

        int NiovAppend = -1;
        if (hasOptionValue("appendiov")) NiovAppend = getOptionValue<int>("appendiov");

        int cnt = -1;
        int cntAppend=0;

        // iterate over the IOVs 
        if(_verbose) std::cout << "iterate over the IOVs "  << std::endl ;

        for (cond::persistency::IOVProxy::Iterator ita = iov.begin(); ita != iov.end() ; ++ita ) {
                cnt++;
                if(_verbose) std::cout << "cnt= "<<cnt  << std::endl ;
                if (cnt == 0) continue;

                if((updateTree && ( (niov>0 && cnt <nIOVTree ) || (NiovAppend>0 && cnt<= nIOVTree) ))||(skipiov>0 && cnt<=skipiov)) {
                        std::cout <<"IOV not written in rootuple (already present or to be skipped )"<< "\n";
                } else {
                        std::cout <<"IOV "<<cnt<<" will be written in rootuple"<< "\n";
                        std::shared_ptr<EcalLaserAPDPNRatios> my_laser 
                                = condDbSession.fetchPayload<EcalLaserAPDPNRatios>( (*ita).payloadId );
                        dbToRoot(*my_laser);
                        ++cntAppend;
                }
                if (niov > 0 && cnt >= niov) break;
                if (NiovAppend > 0 && cntAppend >= NiovAppend) break;
        }

        if(_verbose) std::cout << "number of IOVs appended ="<<cntAppend<<std::endl; 
        if(_verbose) std::cout << "number of IOVs read ="<<cnt<<std::endl; 

        condDbSession.transaction().commit();
        closeRoot();
        return 0;
}


void cond::LaserValidation::dbToRoot(const EcalLaserAPDPNRatios & corr)
{
        const EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap& p = corr.getLaserMap();
        const EcalLaserAPDPNRatios::EcalLaserTimeStampMap& t = corr.getTimeMap();

        if(t.size() != 92) throw cms::Exception("LasCor") << "Unexpected number time parameter triplets\n";

        unsigned t1 = t[0].t1.unixTime();
        unsigned t3 = t[0].t3.unixTime();
        if(_verbose) cout << "Processing IOV " << t1 << " - " << t3
                << "(" << timeToString(t1) << " - " << timeToString(t3) << "\n";

        for(unsigned i = 0; i < t.size(); ++i) {
                time[i]=t[i].t2.unixTime();
        }

        int icrys=0;
        for(int ieta = -EBDetId::MAX_IETA; ieta <= EBDetId::MAX_IETA; ++ieta) {
                if(ieta==0) continue;
                for(int iphi = EBDetId::MIN_IPHI; iphi <= EBDetId::MAX_IPHI; ++iphi) {
                        if (EBDetId::validDetId(ieta,iphi)) {
                                EBDetId detid(ieta,iphi);
                                EcalLaserAPDPNRatios::EcalLaserAPDPNpair corr =  p.barrel(detid.hashedIndex());
                                cor[icrys]   = corr.p2;
                                detId[icrys] = detid.hashedIndex();
                                eta[icrys]   = detid.ieta();
                                phi[icrys]   = detid.iphi();
                                x[icrys]     = 0;
                                y[icrys]     = 0;
                                z[icrys]     = 0;
                                icrys++;
                        }
                }
        }

        for(int iZ = 1; iZ >= -1; --iZ){
                for(int iX = EEDetId::IX_MIN; iX <= EEDetId::IX_MAX; ++iX) {
                        for(int iY = EEDetId::IY_MIN; iY <= EEDetId::IY_MAX; ++iY) {
                                if (EEDetId::validDetId(iX,iY,iZ)) {
                                        EEDetId detid(iX,iY,iZ);
                                        EcalLaserAPDPNRatios::EcalLaserAPDPNpair corr =  p.endcap(detid.hashedIndex());
                                        cor[icrys]   = corr.p2;
                                        detId[icrys] = detid.hashedIndex();
                                        eta[icrys]   = 0;
                                        phi[icrys]   = 0;
                                        x[icrys]     = detid.ix();
                                        y[icrys]     = detid.iy();
                                        z[icrys]     = detid.zside();
                                        icrys++;
                                }
                        }
                }
        }
        tree->Fill();
        std::cout << "filled tree.\n";
}


int main( int argc, char** argv )
{
        cond::LaserValidation valida;
        return valida.run(argc,argv);
}
