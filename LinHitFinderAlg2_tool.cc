#include "dune/LinHitFinder/LinHitFinderAlg2.h"
#include "dune/LinHitFinder/Algorithms.h"

#include <algorithm>

std::vector<short> GetLinFitParameters(short posAmp   ,    // Maximum adc value of the peak of the hit
					short posAmpPos,    // Time of the maximum adc value in tick
					short threshold,    // Threshold value for the current scheme
					short peakStart,    // Time the hit crosses the threshold value on the left in ticks
					short peakEnd       // Time the hit crosses the threshold value on the right in ticks
					) {
    
  // Define the values that you want to get out of the function
  short risingGradient  = 0; short risingIntercept  = 0; short leftZeroPoint  = 0;
  short fallingGradient = 0; short fallingIntercept = 0; short rightZeroPoint = 0;

  // Calcualte the gradients
  risingGradient  = (posAmp -  threshold) / (posAmpPos - peakStart);
  fallingGradient = (posAmp -  threshold) / (posAmpPos - peakEnd  );

  // Calculate the intercepts
  risingIntercept  = posAmp - (posAmpPos * risingGradient );
  fallingIntercept = posAmp - (posAmpPos * fallingGradient);

  // Calculate the zero points
  leftZeroPoint  = (-1 * ( risingIntercept)) / ( risingGradient);
  rightZeroPoint = (-1 * (fallingIntercept)) / (fallingGradient);

  std::vector<short> OutputParameters = {risingGradient  ,
					 risingIntercept ,
					 leftZeroPoint   ,
					 fallingGradient ,
					 fallingIntercept,
					 rightZeroPoint  };
  return OutputParameters;
  
}

LinHitFinderAlg2::LinHitFinderAlg2(fhicl::ParameterSet const & p)
  : fThreshold          (p.get< unsigned int       >("Threshold"              , 20  )),
    fUseSignalKill      (p.get< bool               >("UseSignalKill"          , true)),
    fSignalKillLookahead(p.get< short              >("SignalKillLookahead"    , 5   )),
    fSignalKillThreshold(p.get< short              >("SignalKillThreshold"    , 15  )),
    fSignalKillNContig  (p.get< short              >("SignalKillNContig"      , 1   )),
    fCautiousNContig    (p.get< short              >("CautiousPedestalNContig", 10  )),
    fDoFiltering        (p.get< bool               >("DoFiltering"            , true)),
    fDownsampleFactor   (p.get< unsigned int       >("DownsampleFactor"       , 1   )),
    fFilterTaps         (p.get< std::vector<short> >("FilterCoeffs"           , {2, 9, 23, 31, 23, 9, 2})),
    fMultiplier         (std::accumulate(fFilterTaps.begin(), fFilterTaps.end  (), 0)) {
  std::cout << "Threshold is: " << fThreshold << std::endl;
}

/*
|   _____                         _____                       _ _              
|  |  __ \                       / ____|                     | (_)             
|  | |  | | _____      ___ __   | (___   __ _ _ __ ___  _ __ | |_ _ __   __ _  
|  | |  | |/ _ \ \ /\ / / '_ \   \___ \ / _` | '_ ` _ \| '_ \| | | '_ \ / _` | 
|  | |__| | (_) \ V  V /| | | |  ____) | (_| | | | | | | |_) | | | | | | (_| | 
|  |_____/ \___/ \_/\_/ |_| |_| |_____/ \__,_|_| |_| |_| .__/|_|_|_| |_|\__, | 
|                                                      | |               __/ | 
|                                                      |_|              |___/
*/

std::vector<short> LinHitFinderAlg2::downSample(const std::vector<short>& origional) {
  if (fDownsampleFactor == 1) {
    return origional;
  }
  else {
    std::vector<short> waveform;
    for (size_t i=0; i<origional.size(); i += fDownsampleFactor) {
      waveform.push_back(origional[i]);
    }
    return waveform;
  }
}

/*
|   ____  _             _   _    _ _ _     ______ _           _ _             
|  |  _ \(_)           | | | |  | (_) |   |  ____(_)         | (_)            
|  | |_) |_ _ __   ___ | | | |__| |_| |_  | |__   _ _ __   __| |_ _ __   __ _ 
|  |  _ <| | '_ \ / _ \| | |  __  | | __| |  __| | | '_ \ / _` | | '_ \ / _` |
|  | |_) | | |_) | (_) | | | |  | | | |_  | |    | | | | | (_| | | | | | (_| |
|  |____/|_| .__/ \___/|_| |_|  |_|_|\__| |_|    |_|_| |_|\__,_|_|_| |_|\__, |
|          | |                                                           __/ |
|          |_|                                                          |___/ 
*/

void LinHitFinderAlg2::hitFinding(const std::vector<short>& waveform,
				  std::vector<LinHitFinderAlgorithm::Hit>& hits,
				  int channel) {

  double pedestal = 0;
  
  bool  isPosHit  = false; bool  isNegHit  = false;
  bool  wasPosHit = false; bool  wasNegHit = false;

  int AmplitudePos = 0; int AmplitudeNeg = 0;
  int PeakTimePos  = 0; int PeakTimeNeg  = 0;
  
  LinHitFinderAlgorithm::Hit hit(channel, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  for (size_t iSample=0; iSample<waveform.size(); ++iSample) {

    int   sampleTime = iSample * fDownsampleFactor;
    short adc        = waveform[iSample];

    // Pedestal calculation
    if (iSample == 4) pedestal  = (double)adc;
    if (iSample  > 4) pedestal += ((double)adc - pedestal)/((double)iSample);

    // See if we're in a hit peak!
    // After 50 ticks the pedestal has pretty much converged
    // Recast the value of the individial ADCs to ADC - pedestal
    if (iSample > 50) {
      adc -= (short)pedestal;
      isPosHit = adc >  (short)fThreshold;
      isNegHit = adc < -(short)fThreshold;
    }

/*
|   _____          _ _   _           
|  |  __ \        (_) | (_)          
|  | |__) |__  ___ _| |_ ___   _____ 
|  |  ___/ _ \/ __| | __| \ \ / / _ \
|  | |  | (_) \__ \ | |_| |\ V /  __/
|  |_|   \___/|___/_|\__|_| \_/ \___|
|                                    
|                                    
*/
    // Find the positive hit
    // record all of the starting data about it
    if (isPosHit && !wasPosHit) {
      hit.startTimePos = sampleTime                ;
      hit.chargePos    = adc                       ;
      hit.timeOverThresholdPos += fDownsampleFactor;
    }

    // Now inside the Hit, find some primatives
    // Amplitude and the peak time
    if (isPosHit && wasPosHit) {
      if (adc > AmplitudePos) {AmplitudePos = adc; PeakTimePos = sampleTime;}
      hit.chargePos            += adc * fDownsampleFactor;
      hit.timeOverThresholdPos += fDownsampleFactor;
    }

    // Now at the end of the hit
    // record the end time of the hit
    if (!isPosHit && wasPosHit) {
      hit.endTimePos = sampleTime;
    }

    wasPosHit = isPosHit;

/*
|   _   _                  _   _           
|  | \ | |                | | (_)          
|  |  \| | ___  __ _  __ _| |_ ___   _____ 
|  | . ` |/ _ \/ _` |/ _` | __| \ \ / / _ \
|  | |\  |  __/ (_| | (_| | |_| |\ V /  __/
|  |_| \_|\___|\__, |\__,_|\__|_| \_/ \___|
|               __/ |                      
|              |___/                       
| 
*/
    // Find the negative hit
    // record all of the starting data about it
    if (isNegHit && !wasNegHit) {
      hit.startTimeNeg = sampleTime;
      hit.chargeNeg += adc;
      hit.timeOverThresholdNeg += fDownsampleFactor;
    }

    // Now inside the Hit, find some of the primatives
    // Amplitude and peak time
    if (isNegHit && wasNegHit) {
      if (adc < AmplitudeNeg) {AmplitudeNeg = adc; PeakTimeNeg = sampleTime;}
      hit.chargeNeg += abs(adc) * fDownsampleFactor;
      hit.timeOverThresholdNeg += fDownsampleFactor;
    }

    // Now at the end of the hit
    // record the end time of the hit
    if (!isNegHit && wasNegHit) {
      hit.endTimeNeg = sampleTime;
      if ((PeakTimeNeg - PeakTimePos) < 50) {hits.push_back(hit);}
    }

    wasNegHit = isNegHit;
    
  }
}
  

std::vector<LinHitFinderAlg2::Hit>
LinHitFinderAlg2::findHits(const std::vector<unsigned int>& ChannelNumbers, 
			   const std::vector<std::vector<short>>& inductionSamples) {

  // Make the storage for the hits you're going to make
  auto hits = std::vector<LinHitFinderAlgorithm::Hit>();
  std::cout << "findHits called with "
	    << inductionSamples.size()
	    << " channels. First channel has "
	    << inductionSamples[0].size()
	    << " samples"
	    << std::endl;

  for (size_t it=0; it<inductionSamples.size(); ++it) {
    const std::vector<short>& origionalWaveform = inductionSamples[it];

    std::vector<short> waveform = downSample(origionalWaveform);
    hitFinding(waveform, hits, ChannelNumbers[it]);

  }

  std::cout << "Returning "      << hits.size() << " hits" << std::endl;
  std::cout << "hits/channel = " << float(hits.size())/inductionSamples.size()    << std::endl;
  std::cout << "hits/tick = "    << float(hits.size())/inductionSamples[0].size() << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(LinHitFinderAlg2)
