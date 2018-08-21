classdef ShimSpecs_IUGM_Prisma_fit < ShimSpecs
%SHIMSPECS_IUGM_Prisma_fit
% 
% Shim System Specifications for Siemens shims on Prisma @UNF 
%    
%    ShimSpecs_IUGM_Prisma_fit is a ShimSpecs subclass
%
% =========================================================================
% Updated::20180821::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecs_IUGM_Prisma_fit(  )
%SHIMSPECS - Shim System Specifications 

Shim.Id.systemName   = 'IUGM_Prisma_fit' ;
Shim.Id.channelNames = { 'X' ; 'Y' ; 'Z' ; 'A20' ; 'A21' ; 'B21' ; 'A22' ; 'B22' } ;

Shim.Com = [] ;  

% channels 1-3 refer to the gradient offsets, 4-8 to the 2nd order shims
Shim.Amp.nChannels       = 8 ;  
Shim.Amp.nActiveChannels = 8 ;

Shim.Amp.staticChannels  = true( Shim.Amp.nActiveChannels, 1 ) ;  
Shim.Amp.dynamicChannels = false( Shim.Amp.nActiveChannels, 1 ) ;  


% NB: One can use the Siemens commandline AdjValidate tool to get all the values below:
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
Shim.Amp.maxCurrentPerChannel = [ 2300, 2300, 2300, 4959.01, 3551.29, 3503.299, 3551.29, 3487.302 ]' ;


Shim.Dac = [] ;
% Shim.Dac.resolution = 16 ; % [bits]
% Shim.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end


