classdef ShimSpecsMartinos7TSiemens < ShimSpecs
%SHIMSPECSMARTINOS7TSIEMENS
% 
% Shim System Specifications for Siemens shims on 7T (Bay 5) at Martinos
%
% .......
%   
% Usage
%
% Specs = ShimSpecsMartinos7TSiemens(  )
%
%   Specs contains fields
%           
%       .Amp    
%           relating to amplifcation
%
%       .Com
%           relating to communication (e.g. RS-232)
%
%       .Dac 
%           relating to digital-to-analog conversion
%   
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%    
%    ShimSpecsRri is a ShimSpecs subclass
%
% =========================================================================
% Updated::20170213::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimSpecsMartinos7TSiemens(  )
%SHIMSPECS - Shim System Specifications 

    
Shims.Com = [] ;  

Shims.Amp.maxCurrentPerChannel = 5 ; % (absolute) [units: amps]

Shims.Amp.nChannels       = 8 ;  
Shims.Amp.nActiveChannels = 8 ;

Shims.Dac = [] ;
% Shims.Dac.resolution = 16 ; % [bits]
% Shims.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

