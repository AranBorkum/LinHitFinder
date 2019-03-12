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
  int threshold = 15;
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

    // Filling the vector of pedestal subtracted ADC values
    // Checking to see if there is a dissernable peak is the dataset
    int    PositivePeakAmplitude = 0    ;
    int    NegativePeakAmplitude = 0    ;
    int    PositivePeakPosition  = 0    ;
    int    NegativePeakPosition  = 0    ;
    int    PositivePeakArea      = 0    ;
    int    NegativePeakArea      = 0    ;
    double PositivePeakRMS       = 0    ;
    double NegativePeakRMS       = 0    ;
    int    BipolStartBin         = 0    ;
    int    BipolMiddleBin        = 0    ;
    int    BipolEndBin           = 0    ;
    bool   ThereIsAPeak          = false;
    bool   ThereIsAnotherPeak    = false;

    double RA = 0;
    for (size_t i=0; i<digit[event].ADCs().size(); ++i) {
      RA += ((double)digit[event].ADCs()[i] - RA)/((double)i + 1);
      std::cout << RA << " ";
    }
					    

    // Calculating the pedestal value
    for (size_t ADC=5; ADC<digit[event].ADCs().size(); ++ADC) 
      pedestal += digit[event].ADCs()[ADC];
    pedestal /= digit[event].ADCs().size() - 5;

    // Filling the dataset with the pedestal subtracted ADC values
    // Assessing whether there is a peak to look at in the dataset
    for (size_t ADC=5; ADC<digit[event].ADCs().size(); ++ADC) {
      int ADCValue = (int)digit[event].ADCs()[ADC] - (int)pedestal;
      data.push_back( ADCValue );
      if (ADCValue >=  (double)threshold) ThereIsAPeak       = true;
      if (ADCValue <= -(double)threshold) ThereIsAnotherPeak = true;
    }

    if (!ThereIsAPeak || !ThereIsAnotherPeak) continue;
    // Filling the output file with the pedistal subtracted ADC values if there
    // is a peak in the data. Threshold value is considered
    // This loop is also going to have the primitive calculations for the sake of
    // speed and efficiency
    for (size_t bin=0; bin<data.size(); ++bin) {

      if (data[bin] > PositivePeakAmplitude) { PositivePeakAmplitude = data[bin]; PositivePeakPosition = bin; }
      if (data[bin] < NegativePeakAmplitude) { NegativePeakAmplitude = data[bin]; NegativePeakPosition = bin; }
      if (abs(NegativePeakPosition - PositivePeakPosition) > 20) continue;
      
    }

    // Calculating the start, middle and end points of the bipole
    BipolStartBin  = PositivePeakPosition;
    BipolMiddleBin = PositivePeakPosition;
    BipolEndBin    = NegativePeakPosition;
    
    while (data[BipolStartBin ] > 0) { BipolStartBin --; }
    while (data[BipolMiddleBin] > 0) { BipolMiddleBin++; }
    while (data[BipolEndBin   ] < 0) { BipolEndBin   ++; }

    for (int adc=(int)BipolStartBin ; adc<(int)BipolMiddleBin; ++adc) PositivePeakArea += abs(data[adc]);
    for (int adc=(int)BipolMiddleBin; adc<(int)BipolEndBin   ; ++adc) NegativePeakArea += abs(data[adc]);

    PositivePeakRMS = (double)(BipolMiddleBin - BipolStartBin ) / 2.0; 
    NegativePeakRMS = (double)(BipolEndBin    - BipolMiddleBin) / 2.0;
    
    if (BipolStartBin != BipolMiddleBin &&
	abs(NegativePeakPosition - PositivePeakPosition) < 20) {

      // Starting with all of the positive peak stuff
      // Space line will seperate the negative peak stuff
      fOutputFilePrims << PositivePeakAmplitude << "\t"
		       << PositivePeakPosition  << "\t"
		       << PositivePeakArea      << "\t"
		       << PositivePeakRMS       << "\t"
		       << NegativePeakAmplitude << "\t"
		       << NegativePeakPosition  << "\t" 
		       << NegativePeakArea      << "\t"
		       << NegativePeakRMS       << "\t"
		       << BipolStartBin         << "\t"
		       << BipolMiddleBin        << "\t"
		       << BipolEndBin           << std::endl;

      for (auto const& i: data) 
	fOutputFileData << i << " ";
    }
    fOutputFileData << std::endl;

  }
  
}


DEFINE_ART_MODULE(HitDumper)
