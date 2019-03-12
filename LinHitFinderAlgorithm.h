#ifndef LINHITFINDERALGORITHM_H
#define LINHITFINDERALGORITHM_H

#include <vector>
#include <iostream>

class LinHitFinderAlgorithm {

 public:
  struct Hit {
    Hit(int _channel,
	int _startTimePos, int _chargePos, int _timeOverThresholdPos, int _posAmplitude,
	int _startTimeNeg, int _chargeNeg, int _timeOverThresholdNeg, int _negAmplitude
	)
    : channel             (_channel             ),
      startTimePos        (_startTimePos        ),
      chargePos           (_chargePos           ),
      timeOverThresholdPos(_timeOverThresholdPos),
      posAmplitude        (_posAmplitude        ),
      startTimeNeg        (_startTimeNeg        ),
      chargeNeg           (_chargeNeg           ),
      timeOverThresholdNeg(_timeOverThresholdNeg),
      negAmplitude        (_negAmplitude        )
    {}
    int channel             ;
    int startTimePos        ;
    int chargePos           ;
    int timeOverThresholdPos;
    int posAmplitude        ;
    int startTimeNeg        ;
    int chargeNeg           ;
    int timeOverThresholdNeg;
    int negAmplitude        ;
  };
  
  virtual ~LinHitFinderAlgorithm() = default;
  virtual std::vector<LinHitFinderAlgorithm::Hit>
  findHits(const std::vector<unsigned int>& channelNumbers,
	   const std::vector<std::vector<short>>& inductionSamples) = 0;

};
		
#endif
