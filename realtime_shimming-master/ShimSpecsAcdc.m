classdef ShimSpecsAcdc < ShimSpecs
%SHIMSPECSACDC
% 
% Shim System Specifications for the AC/DC neck coil
%
% .......
%   
% Usage
%
% Specs = ShimSpecsAcdc(  )
%
%   Specs contains fields
%
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
%    ShimSpecsAcdc is a ShimSpecs subclass
%
% =========================================================================
% Updated::20161125::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimSpecsAcdc(  )
%SHIMSPECS - Shim System Specifications 
    
Shims.Com.baudRate    = 57600 ;  
Shims.Com.readTimeout = 500 ; %[units: ms] 

Shims.Com.dataBits    = 8 ;
Shims.Com.stopBits    = 1 ;
Shims.Com.flowControl = 'NONE' ;
Shims.Com.parity      = 'NONE' ;
Shims.Com.byteOrder   = 'bigEndian' ;

% min delay (in seconds) between transmission and reception of data is 0.001 s
Shims.Com.txRxDelay       = 0.001 ; % [units: s]

Shims.Amp.maxCurrentPerChannel = [] ; % (absolute) [units: amps]
Shims.Amp.maxCurrentPerBank    = [] ; % (absolute) [units: amps]

Shims.Amp.nChannels       = [] ;  
Shims.Amp.nActiveChannels = [] ;

Shims.Dac.resolution = 16 ; % [bits]
Shims.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
