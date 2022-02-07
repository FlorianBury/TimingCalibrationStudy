#ifndef Phase2TrackerBXHistogram_h
#define Phase2TrackerBXHistogram_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

// DQM Histograming
class MonitorElement;
class PixelDigiSimLink;
class SimTrack;
class TrackingParticle;
class SimHit;
class TrackerTopology;
class PixelDigi;
class Phase2TrackerDigi;
class TrackerGeometry;
class TF1;
class TTree;

namespace CLHEP {
  class HepRandomEngine;
  class RandGaussQ;
  class RandFlat;
}  // namespace CLHEP


class Phase2TrackerBXHistogram : public DQMEDAnalyzer{

    public:

        explicit Phase2TrackerBXHistogram(const edm::ParameterSet&);
        ~Phase2TrackerBXHistogram() override;
        void bookHistograms(DQMStore::IBooker & ibooker, edm::Run const &  iRun ,
                edm::EventSetup const &  iSetup ) override;
        void dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) override;
        void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
        template <typename T>
        void runSimHit(T isim,double offset, const TrackerTopology* tTopo,const TrackerGeometry* tGeom);

        bool isPixel(const DetId& detId);
        bool isStrip(const DetId& detId);

        struct HistModes{
            MonitorElement* Sampled;
            MonitorElement* Latched;
        };

    private:

        float getSimTrackPt(EncodedEventId event_id, unsigned int tk_id);

        edm::ParameterSet config_;

        std::map<double, HistModes> offsetBX_;
        HistModes offsetBXMap_;
        HistModes attBXMap_;
        HistModes hitsTrueMap_;
        HistModes bookBXHistos(DQMStore::IBooker & ibooker,double offset);

        // Select Hit 

        enum { SquareWindow, SampledMode, LatchedMode, SampledOrLatchedMode, HIPFindingMode };
        double cbc3PulsePolarExpansion(double x) const;
        double signalShape(double x) const;
        double getSignalScale(double xval) const;
        void storeSignalShape();
        bool select_hit(float charge, int bx, float toa, DetId det_id, int hitDetectionMode);
        bool select_hit_sampledMode(float charge, int bx, float toa, DetId det_id, float threshold) const;
        bool select_hit_latchedMode(float charge, int bx, float toa, DetId det_id, float threshold) const;

        std::vector<double> pulseShapeVec_;
        static constexpr float bx_time{25};
        static constexpr size_t interpolationPoints{1000};
        static constexpr int interpolationStep{10};

        //

        std::string geomType_;

        edm::InputTag simTrackSrc_;
        edm::InputTag simVertexSrc_;
        std::vector<edm::InputTag> pSimHitSrc_;
        std::vector<edm::InputTag> pMixSimHitSrc_;
        const edm::InputTag puSummarySrc_;
        edm::InputTag tParticleSrc_;

        const edm::EDGetTokenT<std::vector<TrackingParticle> > tParticleToken_;
        edm::Handle<std::vector<TrackingParticle> > tParticleHandle_;
        std::vector< std::pair<edm::InputTag,edm::EDGetTokenT< edm::PSimHitContainer >> > simHitTokens_;
        std::vector< std::pair<edm::InputTag,edm::EDGetTokenT< CrossingFrame<PSimHit> >> > mixSimHitTokens_;
        const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
        edm::Handle<edm::SimTrackContainer> simTrackHandle_;

        edm::ESHandle<TrackerTopology> tTopoHandle_;
        edm::ESHandle<TrackerGeometry> geomHandle_;

        std::vector<double> pulseShapeParameters_;
        bool use_mixing_;
        std::string mode_;
        bool fireRandom_;
        int bx_range_;
        float deadTime_;
        double pTCut_;
        double theTofLowerCut_;
        double theTofUpperCut_;
        double offset_min_;
        double offset_max_;
        double offset_step_;
        double offset_emulate_;
        std::vector<double> offset_scan_;
        double theThresholdInE_Endcap_;
        double theThresholdInE_Barrel_;
        double theThresholdSmearing_Endcap_;
        double theThresholdSmearing_Barrel_;
        double tof_smearing_;

        float smearEndcapThresholdDetId(DetId);
        float smearBarrelThresholdDetId(DetId);
        float smearToFDetId(DetId);
        
        // Threshold gaussian smearing:
        std::unique_ptr<CLHEP::RandFlat> smearedThreshold_Endcap_;                          
        std::unique_ptr<CLHEP::RandFlat> smearedThreshold_Barrel_;
        std::unique_ptr<CLHEP::RandFlat> smearedTOF_;
        
        // smearing vectors 
        std::map<DetId,float> smearedThresholdFactors_Endcap_;
        std::map<DetId,float> smearedThresholdFactors_Barrel_;
        std::map<DetId,float> smearedTofFactors_;

        const float GeVperElectron; // 3.7E-09 
        int nEvent;
        int verbosity_;
};
#endif
