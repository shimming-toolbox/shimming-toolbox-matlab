classdef ShimSpecsRri < ShimSpecs
%SHIMSPECSRRI
% 
% Shim System Specifications for the RRI 24-channel shim array
%
% .......
%   
% Usage
%
% Specs = ShimSpecsRri(  )
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
%    ShimSpecsRri is a ShimSpecs subclass
%
% =========================================================================
% Updated::20161125::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecsRri(  )
%SHIMSPECS - Shim System Specifications 

    
Shim.Com.baudRate    = 57600 ;  
Shim.Com.readTimeout = 500 ; %[units: ms] 

% As specified in MXD user manual (9700043-0000) pg.16
Shim.Com.dataBits    = 8 ;
Shim.Com.stopBits    = 1 ;
Shim.Com.flowControl = 'NONE' ;
Shim.Com.parity      = 'NONE' ;
Shim.Com.byteOrder   = 'bigEndian' ;

% As specified in the manual, min delay (in seconds) between transmission and
% reception of data is 0.001 s. However, resetting shims at this interval
% causes a fatal system error (84). 0.25 s appears to be about the min
% delay time that won't eventually result in an error.
Shim.Com.txRxDelay       = 0.001 ; % [units: s]
Shim.Com.mxdUpdatePeriod = 0.250 ; % [units: s]

Shim.Amp.nChannels       = 32 ;  
Shim.Amp.nActiveChannels = 24 ;

Shim.Amp.maxCurrentPerChannel = 5*ones(Shim.Amp.nActiveChannels, 1); % (absolute) [units: amps]
Shim.Amp.maxCurrentPerBank    = 20 ; % (absolute) [units: amps]
Shim.Amp.maxCurrentPerRail    = 10 ; % +/- [units: amps]

Shim.Amp.maxSlices = 1000 ; % specific to DSU. Consider alt class def. 

Shim.Dac.resolution = 16 ; % [bits]
Shim.Dac.maxCurrent = 5 ; % (absolute) [units: amps]
% Shim.Dac.maxTimeConstant = hex2dec('3FFF') ; % is this used anywhere?

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
