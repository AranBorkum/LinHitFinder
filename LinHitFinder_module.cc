////////////////////////////////////////////////////////////////////////
// Class:       LinHitFinder
// Plugin Type: producer (art v2_10_03)
// File:        LinHitFinder_module.cc
//
// Generated at Tue Jun  5 06:03:24 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include "dune/LinHitFinder/LinHitFinderAlgorithm.h"

#include <memory>

class LinHitFinder;

class LinHitFinder : public art::EDProducer {

public:
  explicit LinHitFinder(fhicl::ParameterSet const & p);

  LinHitFinder(LinHitFinder const &) = delete;
  LinHitFinder(LinHitFinder      &&) = delete;
  LinHitFinder & operator = (LinHitFinder const &) = delete;  
  LinHitFinder & operator = (LinHitFinder      &&) = delete;

  void produce(art::Event & e) override;

private:
  std::string                            fInputTag;
  std::unique_ptr<LinHitFinderAlgorithm> fFinder  ;

};

LinHitFinder::LinHitFinder(fhicl::ParameterSet const & p)
  : fInputTag(p.get<std::string>("InputTag", "daq")),
    fFinder{art::make_tool<LinHitFinderAlgorithm>(p.get<fhicl::ParameterSet>("finder"))} {

  produces<std::vector<recob::Hit>>();
  produces<art::Assns<raw::RawDigit, recob::Hit>>();
}

void LinHitFinder::produce(art::Event & e) {

  // Grabbing the data out of the file
  auto const& digits_handle = e.getValidHandle<std::vector<raw::RawDigit>>(fInputTag);
  auto& digits = * digits_handle;

  // Getting the geometry of the detector
  art::ServiceHandle<geo::Geometry>                geo             ;
  std::vector<std::vector<short>>                  inductionSamples;
  std::vector<unsigned int>                        channelNumbers  ;
  std::map<raw::ChannelID_t, const raw::RawDigit*> channelToDigit  ;

  for (auto&& digit: digits) {
    const geo::SigType_t signalType = geo->SignalType(digit.Channel());
    if (signalType == geo::kInduction) {
      channelToDigit[digit.Channel()] = &digit;
      channelNumbers.push_back(digit.Channel());
      inductionSamples.push_back(digit.ADCs());
    }
  }

  // Do the actuall hit finding with the induction hit finding algorithm
  std::vector<LinHitFinderAlgorithm::Hit> hits = fFinder->findHits(channelNumbers, inductionSamples);

  // Now actually make the hit objects
  recob::HitCollectionCreator hitCollection(*this,
					    e    ,
					    false,    // doWireAssns
					    true );   // doRawDigitAssns
  for (auto const& hit: hits) {
    const raw::RawDigit* digit=channelToDigit[hit.channel];
    if (!digit) 
      std::cout << "No digit with channel "
		<< hit.channel
		<< " found. Please set the channel correctly"
		<< std::endl;
    std::vector<geo::WireID> wireIDs = geo->ChannelToWire(hit.channel);
    geo::WireID wireID = wireIDs[0];

    recob::HitCreator lar_hit(*digit                                   ,  //RAW DIGIT REFERENCE
			      wireID                                   ,  //WIRE ID		 
			      hit.startTimePos                         ,  //START TICK	 
			      hit.startTimePos+hit.timeOverThresholdPos,  //END TICK 	 
			      hit.timeOverThresholdPos                 ,  //RMS		 
			      hit.startTimePos                         ,  //PEAK_TIME	 
			      0                                        ,  //SIGMA_PEAK_TIME	 
			      0                                        ,  //PEAK_AMPLITUDE	 
			      0                                        ,  //SIGMA_PEAK_AMPLITUDE
			      hit.chargePos                            ,  //HIT_INTEGRAL	 
			      0                                        ,  //HIT_SIGMA_INTEGRAL
			      hit.chargePos                            ,  //SUMMED CHARGE 	 
			      0                                        ,  //MULTIPLICITY	 
			      0                                        ,  //LOCAL_INDEX	 
			      0                                        ,  //WIRE ID		 
			      0                                        ); //DEGREES OF FREEDOM

    hitCollection.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
  }
  hitCollection.put_into(e);
}

DEFINE_ART_MODULE(LinHitFinder)
