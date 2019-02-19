# realtime_shimming

- [Overview](#overview)


## Overview

This library consists of programs to perform shimming (static and real-time).
First designed for use with the 24-channel spine shim (Topfer R, et al., MRM,
2018. https://doi.org/10.1002/mrm.27089)

## Dependencies

- Matlab R2015A
  - Optimization toolbox
  - Image processing toolbox
- [SCT v 4.0.0](https://github.com/neuropoly/spinalcordtoolbox)

## Installation

- Download (or `git clone`) this repository


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

## Notes

When contributing to the library, for maintainability, please adhere to the
Matlab programming guidelines outlined by Richard Johnson:
https://www.mathworks.com/matlabcentral/fileexchange/2529-matlab-programming-style-guidelines

For general information on OO-programming in Matlab:
http://www.cs.ubc.ca/~murphyk/Software/matlabTutorial/html/objectOriented.html

All classes are * handle * as opposed to * value * classes.
For more info on this distinction:
https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html  


## Contributors

[List of contributors.](https://github.com/neuropoly/realtime_shimming/graphs/contributors)

## License

The MIT License (MIT)

Copyright (c) 2018 Polytechnique Montreal, Université de Montréal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
