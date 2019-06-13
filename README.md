# realtime_shimming

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Getting started](#getting-started)
- [Class definitions](#class-definitions)
- [Contributors](#contributors)
- [License](#license)

## Overview

This library consists of programs to perform shimming (static and real-time).
First designed for use with the 24-channel spine shim ([Topfer R, et al., MRM,
2018](https://doi.org/10.1002/mrm.27089)).

## Installation

Before running this software you will need to install the following dependencies:
- MATLAB (tested on R2015A, but more recent are expected to work)
  - Optimization toolbox
  - Image processing toolbox
- [SCT v 4.0.0](https://github.com/neuropoly/spinalcordtoolbox)

Download (or `git clone`) this repository and add this folder (with sub-folders) to the Matlab path. 
For certain features (e.g. recording from a respiratory sensor in "daemon mode"),
the Matlab path should be configured automatically at session start-up by adding the following lines to '~/startup.m' 
(create the file if it does not exist):
 
```
% Change ~/Code/realtime_shimming_repository/ to wherever the repository was downloaded/cloned:
PATH_TO_REALTIME_SHIMMING_REPOSITORY = '~/Code/realtime_shimming_repository/' ;
addpath( genpath( PATH_TO_REALTIME_SHIMMING_REPOSITORY ) ) ;
```

For "daemon mode" to work, the MATLAB session must be started from the command line, and from the home folder, e.g.:

```
user@polymtl ~ $ matlab &
```

For the command line start, MATLAB also needs to exist within the system path, e.g. add the following lines (adapted to refer to your version of MATLAB) to ~/.bash_profile

```
# add MATLAB path
export PATH=$PATH:/Applications/MATLAB_R2015a.app/bin/
```

*For phase unwrapping:*

To use the optional Abdul-Rahman 3D phase unwrapper, binaries must be compiled from the source code found in /external/source/
      
      
## Getting started



## Class definitions

Series of classes pertaining to shimming:

**ShimSpecs( )**

*System specifications Re: amplifier, DAC, etc.*

**ShimCom( )**

*Low-level communication with the amplifier/hardware.*

**ShimOpt( )**

*Optimization of shim currents. Subclass of FieldEval.*

**ShimUse( )**

*Shim control via high-level commands.*

**AuxTracking( )**

*Dynamic tracking of resonance offsets (for now: only subject respiration vis-a-vis the respiratory probe...).*

For the documentation, in the Matlab command prompt type:
	doc [class name]

If you would like to register a new shim system, create a new folder Shim_MyNewShim, and use the .m classes provided in the folder **Shim_Template/**.

## Contributors

[List of contributors.](https://github.com/neuropoly/realtime_shimming/graphs/contributors)

## License

The MIT License (MIT)

Copyright (c) 2018 Polytechnique Montreal, Université de Montréal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
