#ifndef NOISE_LASER_DUMPER
#define NOISE_LASER_DUMPER
#include "CondCore/Utilities/interface/Utilities.h"
#include "CondCore/CondDB/interface/ConnectionPool.h"
#include "CondCore/CondDB/interface/IOVProxy.h"

#include "CondFormats/Common/interface/TimeConversions.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"

#include "usercode/DBDump/interface/EcalLaserPlotter.h"

#include <boost/program_options.hpp>
#include <iterator>
#include <iostream>

using namespace std;

namespace cond {
        class NoiseLaserDumper : public Utilities {
                public:
            cond::persistency::Session sessionMain;

            struct APDPN{
                double p1;
                double p2;
                double p3;
                double ix;
                double iy;
                double iz;

                int id;
            };
            struct NoiseStruct{
                vector<APDPN> noise_apdpn;
                double IOVBegin;
                double IOVEnd;
                //int IOV;
                vector<double> Alpha;
            };
            APDPN struct_apdpn;
            NoiseStruct a;
            vector<NoiseStruct> NoiseVariables;
            vector<APDPN> v_apdpn;
            cond::Time_t TimeLimit;
            cond::Time_t TimeLimitStart;
	    cond::Time_t MyTime;

            struct Coord {
                    int ix_;
                    int iy_;
                    int iz_;
            } _c;





                        typedef EcalLaserAPDPNRatios A;
                        typedef EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap AMap;
                        typedef EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap::const_iterator AMapCit;
                        typedef EcalLaserAPDPNRatios::EcalLaserAPDPNpair AP;
                        typedef EcalLaserAPDPNRatios::EcalLaserTimeStamp AT;

                        NoiseLaserDumper();
                        ~NoiseLaserDumper();
                        int execute();
                        int executeClone(std::string UserTag, cond::persistency::Session session);
                        void dump_txt(char * filename, const A & obja, time_t tb, time_t te, bool flat);
                        void coord(DetId id);
                private:
                        std::vector<DetId> ecalDetIds_;
        };

}

cond::NoiseLaserDumper::NoiseLaserDumper():Utilities("cmscond_list_iov")
{
        addAuthenticationOptions();
        addOption<bool>("verbose","v","verbose");
        addOption<bool>("all","a","list all tags (default mode)");
        addOption<std::string>("connect","c","connect to a specific db");
        addOption<cond::Time_t>("beginTime","b","begin time (first since) (optional)");
        addOption<cond::Time_t>("endTime","e","end time (last till) (optional)");
        addOption<std::string>("tag","t","first tag to be compared");
        addOption<int>("niov","n","number of IOV");
        addOption<int>("prescale","s","prescale factor");
        addOption<std::string>("dump", "d", "file for a txt dump");
        addOption<bool>("flat", "f", "flat iov, i.e. p1 = p2 = p3");

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

cond::NoiseLaserDumper::~NoiseLaserDumper(){
}

void cond::NoiseLaserDumper::coord(DetId id)
{
        if (id.subdetId() == EcalBarrel) {
                EBDetId eid(id);
                _c.ix_ = eid.ieta();
                _c.iy_ = eid.iphi();
                _c.iz_ = 0;
        } else if (id.subdetId() == EcalEndcap) {
                EEDetId eid(id);
                _c.ix_ = eid.ix();
                _c.iy_ = eid.iy();
                _c.iz_ = eid.zside();
        } else {
                fprintf(stderr, "[coord] ERROR: invalid DetId %d", id.rawId());
                assert(0);
        }
}
void cond::NoiseLaserDumper::dump_txt(char * filename, const A & obja, time_t tb, time_t te, bool flat)
{
        FILE * fd = fopen(filename, "w");
        if (!fd) {
                char err[256];
                sprintf(err, "[cond::NoiseLaserDumper::dump] Impossible to open file `%s' for dumping:", filename);
                perror(err);
                exit(1);
        }
        AT ta;
        AMapCit ita, itb;
        fprintf(fd, "T %ld %ld", tb, te);
        printf("T %ld %ld", tb, te);
        printf ("The current since time is: %s", ctime (&tb));
        printf ("The current till time is: %s", ctime (&te));
        for (size_t i = 0; i < 92; ++i) {
                ta = obja.getTimeMap()[i];
                fprintf(fd, " %d", ta.t2.unixTime());
                printf(" %d", ta.t2.unixTime());
        }
        fprintf(fd, "\n");
        for (size_t i = 0; i < ecalDetIds_.size(); ++i) {
                ita = obja.getLaserMap().find(ecalDetIds_[i]);
                if (!flat) fprintf(fd, "P %d %f %f %f\n", ecalDetIds_[i].rawId(), ita->p1, ita->p2, ita->p3);
                else       fprintf(fd, "P %d %f %f %f\n", ecalDetIds_[i].rawId(), ita->p2, ita->p2, ita->p2);
                coord(ecalDetIds_[i]);
                struct_apdpn.id = ecalDetIds_[i].rawId();
                struct_apdpn.p1 = ita->p1;
                struct_apdpn.p2 = ita->p2;
                struct_apdpn.p3 = ita->p3;
                struct_apdpn.ix = _c.ix_;
                struct_apdpn.iy = _c.iy_;
                struct_apdpn.iz = _c.iz_;
                v_apdpn.push_back(struct_apdpn);

        }
        a.noise_apdpn = v_apdpn;
}



int cond::NoiseLaserDumper::execute()
{
        std::string connect = getOptionValue<std::string>("connect" );
        cond::persistency::ConnectionPool connPool;
        if( hasOptionValue("authPath") ){
                connPool.setAuthenticationPath( getOptionValue<std::string>( "authPath") );
        }
        connPool.configure();
        sessionMain = connPool.createSession( connect );


        return 0;
}
int cond::NoiseLaserDumper::executeClone(std::string UserTag, cond::persistency::Session session)
{
    //std::string tag1 = getOptionValue<std::string>("tag");
    std::string tag1 = UserTag;
    session.transaction().start( true );
    const cond::persistency::IOVProxy & iov = session.readIov(tag1, true);

    cond::Time_t since = std::numeric_limits<cond::Time_t>::min();
    since =(time_t) TimeLimitStart;

    if(hasOptionValue("beginTime")) since = getOptionValue<cond::Time_t>("beginTime");
    cond::Time_t till = std::numeric_limits<cond::Time_t>::max();
    till = (time_t) TimeLimit;
    if(hasOptionValue("endTime")) till = getOptionValue<cond::Time_t>("endTime");

    FILE * fdump = NULL;
    if (hasOptionValue("dump")) {
            fdump = fopen(getOptionValue<std::string>("dump").c_str(), "w");
            assert(fdump);
    }

    //bool verbose = hasOptionValue("verbose");
    bool flat = hasOptionValue("flat");

    std::cout << "since: " << since << "   till: " << till << "\n";
    if(since < TimeLimitStart) return 0;
    if(since > TimeLimit  ) return 0;
    int niov = -1;
    if (hasOptionValue("niov")) niov = getOptionValue<int>("niov");

    int prescale = 1;
    if (hasOptionValue("prescale")) prescale = getOptionValue<int>("prescale");
    assert(prescale > 0);

    int cnt = 0, cnt_iov = 0;
    char filename[256];
    A res;

    for (const auto & i : iov) {
            ++cnt_iov;
            if (i.since < since || i.till > till) continue;
	    if (MyTime < i.since || MyTime > i.till) continue;
            if (cnt_iov % prescale != 0) continue;
            ++cnt;
            a.IOVBegin = i.since;
            a.IOVEnd = i.till;
            std::cout << cnt_iov << " " << i.since << " -> " << i.till << " " << cnt << "\n";
            std::shared_ptr<A> pa = session.fetchPayload<A>(i.payloadId);
            sprintf(filename, "dump_%s__since_%08ld_till_%08ld.dat", tag1.c_str(), (long int)i.since, (long int)i.till);
            fprintf(stderr, "going to dump on `%s'\n", filename);

            dump_txt(filename, *pa, (time_t)i.since>>32, (time_t)i.till>>32, flat);
            NoiseVariables.push_back(a);

            if (niov > 0 && cnt >= niov) break;
    }
    session.transaction().commit();
    return 0;
}

#endif


//int main( int argc, char** argv )
//{
//        cond::NoiseLaserDumper valida;
//        return valida.run(argc,argv);
//}
