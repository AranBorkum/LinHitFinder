#include "dune/LinHitFinder/LinHitFinderAlg2.h"
#include "dune/LinHitFinder/Algorithms.h"

#include <algorithm>

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
  short posAmp    = 0    ; short negAmp    = 0    ;
  short posAmpPos = 0    ; short negAmpPos = 0    ;

  LinHitFinderAlgorithm::Hit hit(channel, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
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

    // Inside the positive hit
    if (isPosHit && !wasPosHit) {
      hit.startTimePos          = sampleTime       ;
      hit.chargePos             = adc              ;
      hit.timeOverThresholdPos += fDownsampleFactor;
    }
    if (isPosHit && wasPosHit) {
      hit.chargePos            += adc * fDownsampleFactor;
      hit.timeOverThresholdPos += fDownsampleFactor      ;
      if (adc > posAmp) {
	posAmp = adc;
	posAmpPos = sampleTime;
      }
    }
    if (!isPosHit && wasPosHit) {
      hit .chargePos /= fMultiplier;
      hit .posAmplitude = posAmp;
      hit .posAmpPosition = posAmpPos;
    }

    // Inside the negative hit
    if (isNegHit && !wasNegHit) {
      hit.startTimeNeg          = sampleTime       ;
      hit.chargeNeg             = abs(adc)         ;
      hit.timeOverThresholdNeg += fDownsampleFactor;
    }
    if (isNegHit && wasNegHit) {
      hit.chargeNeg            += abs(adc) * fDownsampleFactor;
      hit.timeOverThresholdNeg += fDownsampleFactor           ;
      if (abs(adc) > negAmp) {
	negAmp = abs(adc);
	negAmpPos = sampleTime;
      }
    }
    if (!isNegHit && wasNegHit) {
      hit .chargeNeg /= fMultiplier;
      hit .negAmplitude = negAmp;
      hit .negAmpPosition = negAmpPos;

      if ((negAmpPos - posAmpPos) < 50)
	hits.push_back(hit);
    }
    
    wasPosHit = isPosHit;
    wasNegHit = isNegHit;

/*
| The hit finder works! Or at least it is actually finding and making hits.
| Things to do now:
|   - Now I need to actually impliment my hit finding method because I think it's better
|   - Update the hit class from LinHitAlgorithm to include all of the other primatives
|   - Neaten up this code because it looks pretty ugly
|   - See if the noise count goes down after implimenting my algorithm (see HitDumper_module.cc)
|   - Produce trees for positive peaks, negative peaks, bipole peaks (purely for completeness)
|   - Run over more MCC11 files to see what larger statistics will look like
|   - Do the clustering and see what the clustering efficiency will turn out like
|
|   - REALLY IMPORTANT: include a similar procedure for the collection hits. Should be really easy
*/    
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
