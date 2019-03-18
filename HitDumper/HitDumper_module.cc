////////////////////////////////////////////////////////////////////////
// Class:       HitDumper
// Plugin Type: producer (art v2_10_03)
// File:        HitDumper_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "lardataobj/RawData/RawDigit.h"

#include <memory>
#include <fstream>

class HitDumper;

class HitDumper : public art::EDAnalyzer {

public:
  explicit HitDumper(fhicl::ParameterSet const & p);

  HitDumper(HitDumper const &) = delete;
  HitDumper(HitDumper      &&) = delete;
  HitDumper & operator = (HitDumper const &) = delete;
  HitDumper & operator = (HitDumper      &&) = delete; 

  void analyze(art::Event const& e) override;

private:
  std::string   fInputTag           ;
  std::string   fOutputFilenameData ;
  std::string   fOutputFilenamePrims;
  std::ofstream fOutputFileData     ;
  std::ofstream fOutputFilePrims    ;
};

HitDumper::HitDumper(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fInputTag           (p.get< std::string >("InputTag", "daq")),
    fOutputFilenameData (p.get< std::string >("OutputFileData" )),
    fOutputFilenamePrims(p.get< std::string >("OutputFilePrims")),
    fOutputFileData     (fOutputFilenameData),
    fOutputFilePrims    (fOutputFilenamePrims)
{}

void HitDumper::analyze(art::Event const& e) {

  // Definition of the threshold, pedestal and pedestal subtracted data values
  int threshold = 20;
  int pedestal = 0;
  std::vector<int> data;
  
  // Get the vector of the vectors of ADC values
  auto const& digits_handle = e.getValidHandle<std::vector<raw::RawDigit>>(fInputTag);
  auto&&      digit         = *digits_handle;

  // Make an alias for the geometry of the detector
  art::ServiceHandle<geo::Geometry> geo;

  // Set up storage for the trigger primatives we calculate
  // std::vector<int> Amplitudes;	
  // std::vector<int> PeakTimes ;

  // Looping over all of the events in the vector of vector of ADC values
  for (size_t event=0; event<digit.size(); ++event) {

    // clearing the vector of ADC values
    data.clear();
    
    // Get the plane that this hit is registered on; true represents Induction plane
    // Only wanting to look for Induction wire hit
    bool IsInduction = geo->SignalType(digit[event].Channel()) == geo::kInduction;

    if (!IsInduction) continue;

    int startTick = 0;
    int hitPeak   = 0;
    int peakTime  = 0;
    int endTick   = 0;
    bool ThereIsAPeak = false;
    bool ThereIsAnotherPeak = false;

    // Calculating the pedestal value
    for (size_t ADC=5; ADC<digit[event].ADCs().size(); ++ADC) 
      pedestal += digit[event].ADCs()[ADC];
    pedestal /= digit[event].ADCs().size() - 5;

    // Filling the dataset with the pedestal subtracted ADC values
    // Assessing whether there is a peak to look at in the dataset
    for (size_t ADC=5; ADC<digit[event].ADCs().size(); ++ADC) {
      int ADCValue = (int)digit[event].ADCs()[ADC] - (int)pedestal;
      data.push_back( ADCValue );
      if (ADCValue >=  threshold) ThereIsAPeak       = true;
      if (ADCValue <= -threshold) ThereIsAnotherPeak = true;
    }

    if (!ThereIsAPeak || !ThereIsAnotherPeak) continue;

    // Fill the data file output
    for (auto const& i: data) {
      fOutputFileData << i << " ";
    }
    fOutputFileData << std::endl;

    bool isHit  = false;
    bool wasHit = false;

    // Loop over the data in the waveform
    for (size_t adc=0; adc<data.size(); ++adc) {

      isHit = (int)data[adc] > (int)threshold;

      if (isHit && !wasHit) {
	startTick = adc;
      }
      if (isHit && wasHit) {
	if (data[adc] > hitPeak) {
	  hitPeak = (int)data[adc];
	  peakTime = (int)adc;
	}
      }
      if (!isHit && wasHit) {
	endTick = adc;
      }

      wasHit = isHit; 
      
    }

    fOutputFilePrims << startTick << " "
		     << endTick   << " "
		     << peakTime  << " "
		     << hitPeak   << std::endl;
    
  }
  
}


DEFINE_ART_MODULE(HitDumper)
