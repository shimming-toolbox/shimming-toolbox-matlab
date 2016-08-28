classdef ShimSpecs
%SHIMSPECS
% 
% Shim System Specifications
%
%
% Specs = ShimSpecs(  )
%
%   Specs contains fields
%
%           
%       .Amp    
%           pertaining to amplifcation
%
%       .Com
%           pertaining to communication (e.g. RS-232)
%
%       .Dac 
%           pertaining to digital-to-analog conversion
%             
% .......
%
%   Description
%   
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimUse
%    ShimTest 
%     
% =========================================================================
% Updated::20160824::ryan.topfer@polymtl.ca
% =========================================================================

properties   
Amp; % pertaining to amplification
Com; % pertaining to communication (e.g. RS-232)
Dac; % pertaining to digital-to-analog conversion 

end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimSpecs(  )
%SHIMSPECS - Shim System Specifications 

    
Shims.Com.baudRate    = 57600 ;  
Shims.Com.readTimeout = 500 ; %[units: ms] 

% As specified in MXD user manual (9700043-0000) pg.16
Shims.Com.dataBits    = 8 ;
Shims.Com.stopBits    = 1 ;
Shims.Com.flowControl = 'NONE' ;
Shims.Com.parity      = 'NONE' ;
Shims.Com.byteOrder   = 'bigEndian' ;

% As specified in the manual, min delay (in seconds) between transmission and
% reception of data is 0.001 s. However, resetting shims at this interval
% causes a fatal system error (84). 0.125 s appears to be about the min
% delay time that won't eventually result in an error.
Shims.Com.txRxDelay            = 0.001 ; % [units: s]

Shims.Amp.maxCurrentPerChannel = 5 ; % (absolute) [units: amps]
Shims.Amp.maxCurrentPerBank    = 20 ; % (absolute) [units: amps]
Shims.Amp.maxCurrentPerRail    = 10 ; % +/- [units: amps]

Shims.Amp.nChannels       = 32 ;  
Shims.Amp.nActiveChannels = 24 ;

Shims.Amp.maxSlices = 1000 ; % specific to DSU. Consider alt class def. 

Shims.Dac.resolution = 16 ; % [bits]
Shims.Dac.maxCurrent = 5 ; % (absolute) [units: amps]
Shims.Dac.maxTimeConstant = hex2dec('3FFF') ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
