classdef (Abstract) ShimSpecs
%SHIMSPECS - Shim System Specifications
%
% .......
% 
% Description
%
%   ShimSpecs defines all the relevant hardware specifications of a shim system.
%   
% .......
%
% Usage
%
%   Note: ShimSpecs is an Abstract class, i.e. it is not 'used' in itself,
%   rather, its subclasses (e.g. ShimSpecs_Greg) are used when instantiated as
%   objects.
%
%   Specs = ShimSpecs(  )
%
%   Specs contains fields
%
%       .Id
%           system identifiers
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
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties   
Id;  % system identifiers
Amp; % relating to amplification
Com; % relating to communication (e.g. RS-232)
Dac; % relating to digital-to-analog conversion 
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecs(  )
%SHIMSPECS - Shim System Specifications 
% 
% The following default Specs are informed by the 24-channel RRI system.

Shim.Id.systemName   = 'Rriyan' ;
Shim.Id.channelNames = cell(24,1) ;
for iCh = 1 : 24
    Shim.Id.channelNames(iCh) = { ['Ch' num2str(iCh) ] } ; 
end

% ------- 
% COM 
Shim.Com.baudRate    = 57600 ;  
Shim.Com.readTimeout = 500 ; % [units: ms] 

Shim.Com.dataBits    = 8 ;
Shim.Com.stopBits    = 1 ;
Shim.Com.flowControl = 'NONE' ;
Shim.Com.parity      = 'NONE' ;
Shim.Com.byteOrder   = 'bigEndian' ;

% min delay (in seconds) between transmission and reception of data.
Shim.Com.txRxDelay   = 0.001 ; % [units: s]

% ------- 
% AMP 
Shim.Amp.maxCurrentPerChannel = 5 ; % (absolute) [units: amps]
Shim.Amp.maxCurrentPerBank    = 20 ; % (absolute) [units: amps]
% Shim.Amp.maxCurrentPerRail    = 10 ; % +/- [units: amps]

Shim.Amp.nChannels       = 32 ;  
Shim.Amp.nActiveChannels = 24 ;

% ------- 
% DAC
Shim.Dac.resolution = 16 ; % [bits]
Shim.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
