classdef (Abstract) ShimSpecs
%SHIMSPECS
% 
% Shim System Specifications
%
% .......
%   
% Usage
%
% Specs = ShimSpecs(  )
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
%    ShimSpecs is an Abstract class.
%
% =========================================================================
% Updated::20161129::ryan.topfer@polymtl.ca
% =========================================================================

properties   
Amp; % relating to amplification
Com; % relating to communication (e.g. RS-232)
Dac; % relating to digital-to-analog conversion 

end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimSpecs(  )
%SHIMSPECS - Shim System Specifications 
% 
% The following default Specs are informed by the 24-channel RRI system.

% ------- 
% COM 
Shims.Com.baudRate    = 57600 ;  
Shims.Com.readTimeout = 500 ; % [units: ms] 

Shims.Com.dataBits    = 8 ;
Shims.Com.stopBits    = 1 ;
Shims.Com.flowControl = 'NONE' ;
Shims.Com.parity      = 'NONE' ;
Shims.Com.byteOrder   = 'bigEndian' ;

% min delay (in seconds) between transmission and reception of data.
Shims.Com.txRxDelay   = 0.001 ; % [units: s]

% ------- 
% AMP 
Shims.Amp.maxCurrentPerChannel = 5 ; % (absolute) [units: amps]
Shims.Amp.maxCurrentPerBank    = 20 ; % (absolute) [units: amps]
% Shims.Amp.maxCurrentPerRail    = 10 ; % +/- [units: amps]

Shims.Amp.nChannels       = 32 ;  
Shims.Amp.nActiveChannels = 24 ;

% ------- 
% DAC
Shims.Dac.resolution = 16 ; % [bits]
Shims.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
