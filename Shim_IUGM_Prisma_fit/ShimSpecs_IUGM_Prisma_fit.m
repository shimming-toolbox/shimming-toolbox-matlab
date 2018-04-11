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
function Shims = ShimSpecs_IUGM_Prisma_fit(  )
%SHIMSPECS - Shim System Specifications 

    
Shims.Com = [] ;  

Shims.Amp.minLarmorFrequency          = 123100100 ; % [units: Hz]
Shims.Amp.maxLarmorFrequency          = 123265000 ; % [units: Hz]
Shims.Amp.maxGradientOffsetPerChannel = 2300 ; % [units: micro-T/m]
Shims.Amp.maxCurrentPerChannel        = 9.998 ; % [units: A ]

% 1st channel refers to RF transmit freq., 2-4 to the gradient offsets, 5-9 to the 2nd order shims
Shims.Amp.nChannels       = 9 ;  
Shims.Amp.nActiveChannels = 9 ;

Shims.Dac = [] ;
% Shims.Dac.resolution = 16 ; % [bits]
% Shims.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end


