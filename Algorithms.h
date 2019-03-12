#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <vector>

void cautiousUpdate(short& median, int& runningDifference, const short sample, const int nContig) {
  
  if (sample > median) ++runningDifference;
  if (sample < median) --runningDifference;

  if (runningDifference > nContig   ) { ++median; runningDifference = 0; }
  if (runningDifference < -1*nContig) { --median; runningDifference = 0; }
}

std::vector<short> FIRFilering(const std::vector<short>& in,
				 const size_t ntaps,
				 const short* taps) {

  std::vector<short> filtered(in.size(), 0);
  for (size_t i=0; i<in.size(); ++i) {
    for (size_t it=0; it<ntaps; ++it) {
      const size_t index = i>it ? i-it: 0;
      filtered[i] += in[index] * taps[it];
    }
  }
  return filtered;
}

std::vector<short> cautiousPedestalSubtraction(const std::vector<short>& in,
					       const int nContig) {

  short median = in[0];
  std::vector<short> pedestal(in.size(), 0);

  int runningDifference = 0;
  for (size_t i=0; i<in.size(); ++i) {
    short sample = in[i];
    cautiousUpdate(median, runningDifference, sample, nContig);
    pedestal[i] = median;
  }
  return pedestal;
} 

std::vector<short> cautiousPedestalSigKill(const std::vector<short>& raw_in,
                                           const int lookahead, const int threshold,
                                           const int ncontig)
{
    short median=raw_in[0];
    std::vector<short> ped(raw_in.size(), 0);
    int runningDiff=0;
    // Are we updating the median (because we're not within a hit)?
    bool updating=true;
    for(size_t i=0; i<raw_in.size()-lookahead; ++i){
        short s=raw_in[i]; // The current sample
        short sig_cand=raw_in[i+lookahead]; // The sample a little way ahead
        // Do we later go over the threshold?
        bool cand_above=(sig_cand>median+threshold);
        // Are we currently below the threshold?
        bool current_below=(s<median+threshold);
        // Is there an upcoming transition over threshold? If so, freeze the pedestal...
        if(updating && cand_above){
            updating=false;
        }
        // ...until we fall below the pedestal again
        if( (!updating) && (current_below)){
            updating=true;
        }
        // Do the frugal streaming if we're not in a hit
        if(updating){
            cautiousUpdate(median, runningDiff, s, ncontig);
        }
        ped[i]=median;
    }
    // Get the last few samples, which we couldn't do in the main loop
    // because we'd go out-of-bounds
    for(size_t i=raw_in.size()-lookahead; i<raw_in.size(); ++i){
        ped[i]=median;
    }
    return ped;
}

#endif
