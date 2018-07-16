#ifndef NOISE_DUMPER_hh
#define NOISE_DUMPER_hh

#include "CondCore/CondDB/interface/IOVProxy.h"
#include "CondCore/Utilities/interface/Utilities.h"

#include "CondCore/CondDB/interface/ConnectionPool.h"

#include <boost/program_options.hpp>
#include <iterator>
#include <iostream>

#include "CondFormats/BeamSpotObjects/interface/BeamSpotObjects.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/EcalObjects/interface/EcalClusterLocalContCorrParameters.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLinearizationConst.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLutGroup.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLutIdMap.h"
#include "CondFormats/EcalObjects/interface/EcalTPGPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightIdMap.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightGroup.h"
#include "CondFormats/EcalObjects/interface/EcalTPGSlidingWindow.h"
#include "CondFormats/EcalObjects/interface/EcalTPGSpike.h"
#include "CondFormats/ESObjects/interface/ESEEIntercalibConstants.h"
#include "CondFormats/ESObjects/interface/ESGain.h"
#include "CondFormats/ESObjects/interface/ESIntercalibConstants.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include <vector>

namespace cond {

        template<class C>
        class NoiseDumper : public cond::Utilities {

                public:

                struct IC{
                    int ix;
                    int iy;
                    int iz;
                    float ic;
                    int id;
                };
                struct Ped{
                    int ix;
                    int iy;
                    int iz;
                    float mean1;
                    float mean2;
                    float mean3;
                    float rms1;
                    float rms2;
                    float rms3;
                    int id;
                };

                struct NoiseStruct{
                    std::vector<IC> noise_ic;
                    std::vector<Ped> noise_ped;
                    int IOVBegin;
                    int IOVEnd;
                    int time;
                    float adcToGeVEB;
                    float adcToGeVEE;

                    //int IOV;
                };
                IC struct_ic;
                Ped struct_ped;
                NoiseStruct a;
                std::vector<NoiseStruct> NoiseVariables;
                std::vector<IC> v_ic;
                std::vector<Ped> v_ped;
                cond::Time_t TimeLimit;
                cond::Time_t TimeLimitStart;
		cond::Time_t MyTime;

                        NoiseDumper(std::string class_name) : Utilities("conddb_dumper"), _class_name(class_name)
                        {
                                addAuthenticationOptions();
                                addOption<bool>("join", "j", "produce one single output file, where IOVs are separated by double line break and a comment line starting with `#'");
                                addOption<bool>("verbose","v","verbose");
                                addOption<cond::Time_t>("beginTime","b","begin time (first since) (optional)");
                                addOption<cond::Time_t>("endTime","e","end time (last till) (optional)");
                                addOption<int>("niov","n","number of IOVs to be dumped");
                                addOption<int>("prescale","s","prescale factor, i.e. dump 1 in N IOVs");
                                addOption<std::string>("object","O","object to be dumped (required)");
                                addOption<std::string>("output","o","output file");
                                addOption<std::string>("tag","t","tag to be dumped (required)");
                                addOption<std::string>("db","d","database to run the command on;"
                                                       " arg can be an explicit db connection, e.g. for sqlite files,"
                                                       " or one of the following aliases:\n"
                                                       "     dev -> FrontierPrep\n"
                                                       "     pro -> PromptProd\n"
                                                       "     arc -> FrontierArc\n"
                                                       "     int -> FrontierInt\n"
                                                       );

                                _ids.resize(EBDetId::MAX_HASH + 1 + EEDetId::kSizeForDenseIndexing);
                                for (int hi = EBDetId::MIN_HASH; hi <= EBDetId::MAX_HASH; ++hi ) {
                                        EBDetId ebId = EBDetId::unhashIndex(hi);
                                        if (ebId != EBDetId()) {
                                                _ids[hi] = ebId;
                                        }
                                }
                                for ( int hi = 0; hi < EEDetId::kSizeForDenseIndexing; ++hi ) {
                                        EEDetId eeId = EEDetId::unhashIndex(hi);
                                        if (eeId != EEDetId()) {
                                                int idx = EBDetId::MAX_HASH + 1 + hi;
                                                _ids[idx] = eeId;
                                        }
                                }
                                assert(_ids.size() == 75848);
                                assert(_ids.size() == EBDetId::MAX_HASH + 1 + EEDetId::kSizeForDenseIndexing);
                               // NoiseVariables.clear();


                        }

                        ~NoiseDumper()
                        {
                        }


                        void print(int cnt, const cond::Iov_t & iov)
                        {
                                printf("%d  %llu %llu\n", cnt, iov.since, iov.till);

                                time_t s = iov.since>>32;
                                time_t t = iov.till>>32;
                                printf ("The current since time is: %s", ctime (&s));
                                printf ("The current till time is: %s", ctime (&t));
                                char ss[64], st[64];
                                ctime_r(&s, ss);
                                ctime_r(&s, st);

                                ss[strlen(ss) - 1] = '\0';
                                st[strlen(st) - 1] = '\0';
                                printf("%d  %s (%ld)  -->  %s (%ld)\n", cnt, ss, s, st, t);
                                return;
                        }


                        FILE * open_file(char * filename)
                        {
                                FILE * fd = fopen(filename, "a");
                                if (!fd) {
                                        char err[256];
                                        sprintf(err, "[cond::ICDump::dump] Impossible to open file `%s' for dumping:", filename);
                                        perror(err);
                                        exit(1);
                                }
                                fprintf(stderr, "going to dump on `%s'\n", filename);
                                return fd;
                        }
                        std::string connect;
                        std::string db;
                        bool join = false;
                         int niov = -1;
                         int prescale = 1;
                         cond::Time_t since = std::numeric_limits<cond::Time_t>::min();
                         cond::Time_t till = std::numeric_limits<cond::Time_t>::max();
                        cond::persistency::Session sessionMain;

                        // main loop
                        int execute()
                        {
                                 connect = "frontier://FrontierProd/CMS_CONDITIONS";
                                 db = "";
                                if (hasOptionValue("db")) {
                                        if      (db == "dev") connect = "frontier://FrontierPrep/CMS_CONDITIONS";
                                        else if (db == "pro") connect = "frontier://PromptProd/CMS_CONDITIONS";
                                        else if (db == "arc") connect = "frontier://FrontierArc/CMS_CONDITIONS";
                                        else if (db == "int") connect = "frontier://FrontierInt/CMS_CONDITIONS";
                                        else {
                                                connect = getOptionValue<std::string>("db" );
                                        }
                                }

                                cond::persistency::ConnectionPool connPool;

                                if(hasOptionValue("authPath")){
                                        connPool.setAuthenticationPath( getOptionValue<std::string>( "authPath") );
                                }



                                connPool.configure();
                                sessionMain = connPool.createSession(connect);

                                return 0;
                        }
                        int executeClone(std::string UserTag, cond::persistency::Session session )
                        {


                            if (hasOptionValue("niov")) niov = getOptionValue<int>("niov");


                            if (hasOptionValue("prescale")) prescale = getOptionValue<int>("prescale");
                            assert(prescale > 0);

                            if( hasOptionValue("beginTime" )) since = getOptionValue<cond::Time_t>("beginTime");
                            if( hasOptionValue("endTime" )) till = getOptionValue<cond::Time_t>("endTime");
                            since = TimeLimitStart;
                            till = TimeLimit;



                            if (hasOptionValue("join")) join = true;
                            //std::string tag = getOptionValue<std::string>("tag");
                            std::string tag = UserTag;
                            printf("%s %s\n", connect.c_str(), tag.c_str());

                            session.transaction().start( true );
                            cond::persistency::IOVProxy iov = session.readIov(tag, true);
                            std::string obj_type = iov.payloadObjectType();

                            auto first_iov = iov.begin();
                            auto last_iov  = iov.begin();
                            for (int i = 0; i < iov.loadedSize() - 1; ++i) ++last_iov;

                            std::cout << (*first_iov).since << " " << (*last_iov).since << "\n";


                            char filename[512];
                            //sprintf(filename, "dump_%s__since_%08llu_jtill_since_%08llu.dat", obj_type.c_str(), (*first_iov).since, (*last_iov).since);
                            sprintf(filename, "dump_%s__since_%08llu_jtill_since_%08llu.dat", _class_name.c_str(), (*first_iov).since, (*last_iov).since);
                            if (hasOptionValue("output")) sprintf(filename, getOptionValue<std::string>("output").c_str());

                            int cnt = 0, cnt_iov = -1;
                            FILE * fout;
                            if (join) remove(filename);
                            v_ic.clear();
                            v_ped.clear();

                            for (auto i : iov) {
                                    ++cnt_iov;
                                    if (i.since < since || i.till > till) continue;
				    if (MyTime < i.since || MyTime > i.till) continue;

				    // std::cout << "iov: " <<  (int) *i << std::endl;
                                    if (cnt_iov % prescale != 0) continue;
                                    ++cnt;
//                                    if(i.since < TimeLimitStart) continue;
//                                    if(i.since  > TimeLimit  ) continue;
                                    print(cnt_iov, i);
                                    a.IOVBegin = i.since;
                                    a.IOVEnd = i.till;

                                    //a.time = i.since + increment;
                                    //increment = increment + 1;
                                    //std::cout << "time:  "<< a.time << increment << std::endl;
                                    std::shared_ptr<C> pa = session.fetchPayload<C>(i.payloadId);
                                    if (!join) sprintf(filename, "dump_%s__since_%08llu_till_%08llu.dat", _class_name.c_str(), i.since, i.till);
                                    fout = open_file(filename);
                                    if (join) fprintf(fout, "# new IOV: since %llu  till %llu\n", i.since, i.till);

                                    dump(fout, *pa);
                                    NoiseVariables.push_back(a);
                                   // std::cout<< "We are here" << std::endl;

                                    if (join) fprintf(fout, "\n\n");
                                    fclose(fout);
                                    if (niov > 0 && cnt >= niov) break;
                            }

                            session.transaction().commit();
                            session.transaction().rollback();



                            return 0;
                        }


                        struct Coord {
                                int ix_;
                                int iy_;
                                int iz_;
                        } _c;


                        void coord(DetId id)
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


                        void dump(FILE * fd, bool) {}

                        void dump(FILE * fd, EcalCondObjectContainer<float> & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);

                                        EcalCondObjectContainer<float>::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                                fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
                                        fprintf(fd, "%d %d %d %f %u\n", _c.ix_, _c.iy_, _c.iz_, *it, id.rawId());

                                        struct_ic.ix = _c.ix_;
                                        struct_ic.iy = _c.iy_;
                                        struct_ic.iz = _c.iz_;

                                        struct_ic.ic = *it;
                                        struct_ic.id = id.rawId();
                                        v_ic.push_back(struct_ic);

                                }
                                a.noise_ic = v_ic;


                        }

                        void dump(FILE * fd, EcalChannelStatus & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);
                                        EcalChannelStatus::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                                fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
                                        fprintf(fd, "%d %d %d %d %d\n", _c.ix_, _c.iy_, _c.iz_, (*it).getStatusCode(), id.rawId());
                                }
                        }

                        void dump(FILE * fd, EcalTPGLinearizationConst & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);
                                        EcalTPGLinearizationConst::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                                fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
                                        fprintf(fd, "%d %d %d %u %u %u %u %u %u  %u\n", _c.ix_, _c.iy_, _c.iz_,
                                                it->mult_x12, it->shift_x12,
                                                it->mult_x6,  it->shift_x6,
                                                it->mult_x1,  it->shift_x1,
                                                id.rawId());
                                }
                        }

                        void dump(FILE * fd, EcalTPGLutGroup & o)
                        {
                                for (const auto & it : o.getMap()) {
                                        fprintf(fd, "%d %d\n", it.first, it.second);
                                }
                        }

                        void dump(FILE * fd, EcalTPGLutIdMap & o)
                        {
                                for (const auto & it : o.getMap()) {
                                        fprintf(fd, "%d  ", it.first);
                                        const unsigned int * lut = it.second.getLut();
                                        for (int i = 0; i < 1024; ++i) {
                                                fprintf(fd, " %d", lut[i]);
                                        }
                                        fprintf(fd, "\n");
                                }
                        }

                        void dump(FILE * fd, EcalTPGPedestals & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);
                                        EcalTPGPedestals::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                                fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
                                        fprintf(fd, "%d %d %d  %d  %d  %d  %u\n", _c.ix_, _c.iy_, _c.iz_,
                                                (*it).mean_x12,
                                                (*it).mean_x6,
                                                (*it).mean_x1,
                                                id.rawId());
                                }
                        }

                        void dump(FILE * fd, EcalTPGWeightGroup & o)
                        {
                                for (const auto & it : o.getMap()) {
                                        fprintf(fd, "%d %d\n", it.first, it.second);
                                }
                        }

                        void dump(FILE * fd, EcalTPGWeightIdMap & o)
                        {
                                for (const auto & it : o.getMap()) {
                                        fprintf(fd, "%d  ", it.first);
                                        uint32_t w0, w1, w2, w3, w4;
                                        it.second.getValues(w0, w1, w2, w3, w4);
                                        fprintf(fd, " %d %d %d %d %d\n", w0, w1, w2, w3, w4);
                                }
                        }

                        void dump(FILE * fd, EcalTPGSlidingWindow & o)
                        {
                                fprintf(fd, "#id value\n");
                                for (const auto & it : o.getMap()) {
                                        fprintf(fd, "%d %d\n", it.first, it.second);
                                }
                        }

                        void dump(FILE * fd, EcalTPGSpike & o)
                        {
                                fprintf(fd, "#stripId lut\n");
                                for (const auto & it : o.getMap()) {
                                        fprintf(fd, "%d %d\n", it.first, it.second);
                                }
                        }

                        void dump(FILE * fd, EcalPedestals & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);
                                        EcalPedestals::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                               fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
//                                        fprintf(fd, "%d %d %d  %f %f  %f %f  %f %f %u\n", _c.ix_, _c.iy_, _c.iz_,
//                                                (*it).mean(1), (*it).rms(1),
//                                                (*it).mean(2), (*it).rms(2),
//                                                (*it).mean(3), (*it).rms(3),
//                                                id.rawId());

                                        struct_ped.mean1 = (*it).mean(1);
                                        struct_ped.mean2 = (*it).mean(2);
                                        struct_ped.mean3 = (*it).mean(3);
                                        struct_ped.rms1 = (*it).rms(1);
                                        struct_ped.rms1 = (*it).rms(2);
                                        struct_ped.rms1 = (*it).rms(3);
                                        struct_ped.id = id.rawId();
                                        struct_ped.ix = _c.ix_;
                                        struct_ped.iy = _c.iy_;
                                        struct_ped.iz = _c.iz_;
                                        v_ped.push_back(struct_ped);


                                }
                                a.noise_ped = v_ped;
                        }

                        void dump(FILE * fd, EcalADCToGeVConstant & o)
                        {
                                fprintf(fd, "EB= %f  EE= %f\n", o.getEBValue(), o.getEEValue());
                                a.adcToGeVEB = o.getEBValue();
                                a.adcToGeVEE = o.getEEValue();

                        }

                        void dump(FILE * fd, EcalTimeOffsetConstant & o)
                        {
                                fprintf(fd, "EB= %f  EE= %f\n", o.getEBValue(), o.getEEValue());
                        }

                        void dump(FILE * fd, EcalClusterLocalContCorrParameters & a)
                        {
                                const EcalFunctionParameters & p = a.params();
                                fprintf(fd, "# %lu parameter(s)\n", p.size());
                                for (size_t s = 0; s < p.size(); ++s) {
                                        fprintf(fd, " %f", p[s]);
                                }
                                fprintf(fd, "\n");
                        }

                        void dump(FILE * fd, EcalGainRatios & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);
                                        EcalGainRatios::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                                fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
                                        fprintf(fd, "%d %d %d  %f %f  %u\n", _c.ix_, _c.iy_, _c.iz_,
                                                (*it).gain12Over6(), (*it).gain6Over1(),
                                                id.rawId());
                                }
                        }

                        void dump(FILE * fd, EcalPulseShapes & o)
                        {
                                for (size_t i = 0; i < _ids.size(); ++i) {
                                        DetId id(_ids[i]);
                                        EcalPulseShapes::const_iterator it = o.find(id);
                                        if (it == o.end()) {
                                                fprintf(stderr, "Cannot find value for DetId %u", id.rawId());
                                        }
                                        coord(_ids[i]);
                                        fprintf(fd, "%d %d %d  ", _c.ix_, _c.iy_, _c.iz_);
                                        for (int is = 0; is < EcalPulseShape::TEMPLATESAMPLES; ++is) {
                                                fprintf(fd, " %f", it->val(is));
                                        }
                                        fprintf(fd, "  %u\n", id.rawId());
                                }
                        }

                        void dump(FILE * fd, ESEEIntercalibConstants & ic)
                        {
                                fprintf(fd,
                                        "GammaLow0:  %f   AlphaLow0:  %f\n"
                                        "GammaLow1:  %f   AlphaLow1:  %f\n"
                                        "GammaLow2:  %f   AlphaLow2:  %f\n"
                                        "GammaLow3:  %f   AlphaLow3:  %f\n"
                                        "GammaHigh0: %f   AlphaHigh0: %f\n"
                                        "GammaHigh1: %f   AlphaHigh1: %f\n"
                                        "GammaHigh2: %f   AlphaHigh2: %f\n"
                                        "GammaHigh3: %f   AlphaHigh3: %f\n",
                                        ic.getGammaLow0(),  ic.getAlphaLow0(),
                                        ic.getGammaLow1(),  ic.getAlphaLow1(),
                                        ic.getGammaLow2(),  ic.getAlphaLow2(),
                                        ic.getGammaLow3(),  ic.getAlphaLow3(),
                                        ic.getGammaHigh0(), ic.getAlphaHigh0(),
                                        ic.getGammaHigh1(), ic.getAlphaHigh1(),
                                        ic.getGammaHigh2(), ic.getAlphaHigh2(),
                                        ic.getGammaHigh3(), ic.getAlphaHigh3()
                                       );


                        }

                        void dump(FILE * fd, ESGain & g)
                        {
                                fprintf(fd, "gain= %f\n", g.getESGain());
                        }

                        void dump(FILE * fd, ESIntercalibConstants & ic)
                        {

                                for (int i = 0; i < ESDetId::kSizeForDenseIndexing; ++i) {
                                        ESDetId id(ESDetId::detIdFromDenseIndex(i));
                                        ESIntercalibConstants::const_iterator it = ic.find(id);
                                        assert(it != ic.end());
                                        fprintf(fd, "%d %f\n", id.rawId(), *it);



                                }

                        }

                        void dump(FILE * fd, RunInfo & ri)
                        {
                                fprintf(fd, "%d %lld %lld (%s -> %s)\n", 
                                        ri.m_run,
                                        ri.m_start_time_ll,
                                        ri.m_stop_time_ll,
                                        ri.m_start_time_str.c_str(),
                                        ri.m_stop_time_str.c_str());
                        }

                        void dump(FILE * fd, BeamSpotObjects & bs)
                        {
                                fprintf(fd, "z= %f  z_err= %f  sigma_z= %f sigma_z_err= %f\n",
                                        bs.GetZ(), bs.GetZError(), bs.GetSigmaZ(), bs.GetSigmaZError());
                        }


                private:
                        std::string _class_name;
                        std::vector<DetId> _ids;

        };
}

#endif
