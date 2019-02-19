classdef ShimSpecs_Greg < ShimSpecs
%SHIMSPECS_GREG
% 
% Shim System Specifications for the AC/DC neck coil
%
% .......
%   
% Usage
%
% Specs = ShimSpecs_Greg(  )
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
%    ShimSpecs_Greg is a ShimSpecs subclass
%
% =========================================================================
% Updated::20181105::ryan.topfer@polymtl.ca
% =========================================================================

properties
    Adc;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecs_Greg(  )
%SHIMSPECS - Shim System Specifications 

Shim.Id.systemName   = 'Greg' ;
Shim.Id.channelNames = cell(8,1) ;
for iCh = 1 :8 
    Shim.Id.channelNames(iCh) = { ['Ch' num2str(iCh) ] } ; 
end
    
Shim.Com.baudRate    = 115200 ;  
% Shim.Com.readTimeout = 500 ; %[units: ms] 

Shim.Com.dataBits    = 8 ;
Shim.Com.stopBits    = 1 ;
Shim.Com.flowControl = 'NONE' ;
Shim.Com.parity      = 'NONE' ;
Shim.Com.byteOrder   = 'bigEndian' ;

% min delay (in seconds) between transmission and reception of data is 1 s
%
% UNTESTED
Shim.Com.txRxDelay       = 0.005 ; % [units: s]
Shim.Com.updatePeriod    = 0.1 ;

Shim.Amp.nChannels       = 8 ;  
Shim.Amp.nActiveChannels = 8 ;

Shim.Amp.maxCurrentPerChannel = 2.5 * ones( Shim.Amp.nActiveChannels, 1 ) ; ; % (absolute) [units: A]
Shim.Amp.maxVoltagePerChannel = 2500 ; % [units: mV]

Shim.Amp.staticChannels  = true( Shim.Amp.nActiveChannels, 1 ) ;  
Shim.Amp.dynamicChannels = true( Shim.Amp.nActiveChannels, 1 ) ;  

Shim.Adc.mVPerAdcCount = 2 ;

Shim.Dac.resolution       = 8 ; % [bits]
Shim.Dac.referenceVoltage = 1250 ; % [units: mV]
Shim.Dac.maximum          = 26214 ; 


end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
