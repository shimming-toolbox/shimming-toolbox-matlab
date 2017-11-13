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
%       .Adc 
%           relating to analog-to-digital conversion
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
% Updated::20171107::ryan.topfer@polymtl.ca
% =========================================================================

properties
    Adc;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimSpecsAcdc(  )
%SHIMSPECS - Shim System Specifications 
    
Shims.Com.baudRate    = 9600 ;  
% Shims.Com.readTimeout = 500 ; %[units: ms] 

Shims.Com.dataBits    = 8 ;
Shims.Com.stopBits    = 1 ;
Shims.Com.flowControl = 'NONE' ;
Shims.Com.parity      = 'NONE' ;
Shims.Com.byteOrder   = 'bigEndian' ;

% min delay (in seconds) between transmission and reception of data is 1 s
Shims.Com.txRxDelay       = 1 ; % [units: s]

Shims.Amp.maxCurrentPerChannel = 0.400 ; % (absolute) [units: A]
Shims.Amp.maxVoltagePerChannel = 200 ; % (absolute) [units: mV]

Shims.Amp.nChannels       = 8 ;  
Shims.Amp.nActiveChannels = 8 ;

Shims.Adc.mVPerAdcCount = 2 ;

Shims.Dac.resolution = 12 ; % [bits]
Shims.Dac.maxVoltage = 4096 ; % (absolute) [units: mV]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
