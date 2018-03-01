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
    
Shims.Com.baudRate    = 115200 ;  
% Shims.Com.readTimeout = 500 ; %[units: ms] 

Shims.Com.dataBits    = 8 ;
Shims.Com.stopBits    = 1 ;
Shims.Com.flowControl = 'NONE' ;
Shims.Com.parity      = 'NONE' ;
Shims.Com.byteOrder   = 'bigEndian' ;

% min delay (in seconds) between transmission and reception of data is 1 s
Shims.Com.txRxDelay       = 1 ; % [units: s]

Shims.Amp.maxCurrentPerChannel = 2.2 ; % (absolute) [units: A]
Shims.Amp.maxVoltagePerChannel = 2500 ; % (absolute) [units: mV]

Shims.Amp.nChannels       = 8 ;  
Shims.Amp.nActiveChannels = 8 ;

Shims.Adc.mVPerAdcCount = 2 ;

Shims.Dac.resolution = 8 ; % [bits]
Shims.Dac.maxVoltage = 4096 ; % (absolute) [units: mV]

% =========================================================================

Shims.Dac.feedbackcalibrationcoeff1=[0.65909, 0.65, 0.64547, 0.6499, 0.6591, 0.654, 0.65, 0.65];         %Calibration coefficient to transform Voltages to Amps
Shims.Dac.feedbackcalibrationcoeff2 = [19.09, -10.908, -23.636, 20.91, 32.728, 11.308, -42.728, 34.546]; %Calibration coefficient to transform Voltages to Amps

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end