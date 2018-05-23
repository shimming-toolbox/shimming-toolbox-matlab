classdef ShimSpecs_IUGM_Prisma_fit < ShimSpecs
%SHIMSPECS_IUGM_Prisma_fit
% 
% Shim System Specifications for Siemens shims on Prisma @UNF 
%    
%    ShimSpecs_IUGM_Prisma_fit is a ShimSpecs subclass
%
% =========================================================================
% Updated::20180312::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecs_IUGM_Prisma_fit(  )
%SHIMSPECS - Shim System Specifications 
    
Shim.Com = [] ;  

% 1st channel refers to RF transmit freq., 2-4 to the gradient offsets, 5-9 to the 2nd order shims
Shim.Amp.nChannels       = 9 ;  
Shim.Amp.nActiveChannels = 9 ;

Shim.Amp.staticChannels  = true( 1, Shim.Amp.nActiveChannels ) ;  
Shim.Amp.dynamicChannels = false( 1, Shim.Amp.nActiveChannels ) ;  


% NB: One can use the Siemens commandline AdjValidate tool to get all the values below:
Shim.Amp.minLarmorFrequency   = 123100100 ; % [units: Hz]
Shim.Amp.maxLarmorFrequency   = 123265000 ; % [units: Hz]

% the range of possible transmit frequencies (f0Bw) is also in Hz, not A.
f0Bw = (Shim.Amp.maxLarmorFrequency - Shim.Amp.minLarmorFrequency)/2 ; % 82450 Hz
% Shim.Amp.maxCurrentPerChannel
%   First 3 terms correspond to the linear gradient offsets: [units: micro-T/m]
%   Last 5 terms correspond to the 2nd order shims: [units: micro-T/m^2]
%   These last 5 correspond to  9.998 A 
%
%   Though I'd prefer to have all the shim terms in units of A, I don't (yet)
%   know what the corresponding values would be for the linears. Moreover,
%   these "multipole units" are those displayed in the Syngo/3D Shim card,
%   so, it may simply be more convenient to deal with the shim currents in
%   "display units".
Shim.Amp.maxCurrentPerChannel = [ f0Bw, 2300, 2300, 2300, 4959.01, 3551.29, 3503.299, 3551.29, 3487.302 ]' ;


Shim.Dac = [] ;
% Shim.Dac.resolution = 16 ; % [bits]
% Shim.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end


