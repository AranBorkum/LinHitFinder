#ifndef LINHITFINDERALG2_H
#define LINHITFINDERALG2_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "dune/LinHitFinder/LinHitFinderAlgorithm.h"

class LinHitFinderAlg2 : public LinHitFinderAlgorithm {

 public:
  explicit LinHitFinderAlg2(fhicl::ParameterSet const & p);

  virtual std::vector<LinHitFinderAlgorithm::Hit>
    findHits(const std::vector<unsigned int>& channelNumbers,
	     const std::vector<std::vector<short>>& inductionSamples);

 protected:
  std::vector<short> downSample  (const std::vector<short>& origional);
  std::vector<short> findPedestal(const std::vector<short>& origional);
  std::vector<short> filter      (const std::vector<short>& origional);
  void hitFinding(const std::vector<short>& waveform           ,
		  std::vector<LinHitFinderAlgorithm::Hit>& hits,
		  int channel                                  );

  unsigned int       fThreshold          ;             
  bool               fUseSignalKill      ;
  short              fSignalKillLookahead;
  short              fSignalKillThreshold;
  short              fSignalKillNContig  ;
  short              fCautiousNContig    ;
  bool               fDoFiltering        ;
  unsigned int       fDownsampleFactor   ;
  std::vector<short> fFilterTaps         ;
  int                fMultiplier         ;

};

#endif
