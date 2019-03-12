<p align="center"><img width=75% src="https://github.com/AranBorkum/LinHitFinder/blob/master/InductionWireHitFit.png"></p>
![dunetpc_version](https://img.shields.io/badge/dunetpc_version-v08_12_00-brightgreen.svg) ![Dependencies](https://img.shields.io/badge/dependencies-standard-orange.svg)  ![issues](https://img.shields.io/badge/issues-none_known-orange.svg) ![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg) ![LArSoft](https://img.shields.io/badge/LArSoft-v08_12_00-brightgreen.svg) 
# Hit finder module for induction wires

# Basic Overview:
  - Specifically targets the induction wires
  - Extracts hit primatives from bipolar induction wire hits
  - Returns a hit object containing these primatives
 
# Features of the module
#### Down Sampling:
  - Takes in the origional stream of data and outputs every n-th entry
#### Pedestal Subtraction:
  - Calcualates the pedestal of the incoming data stream and subtracts it from every ADC value 
#### Filtering:
  - Does FIR filtering to reduce the noise levels of the incoming data stream
#### Note on current version:
  - This is a very basic version of what the module will be
  - Strongly based on the TriggerPrimativeFinder by PRogriguez
