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
%   
%
%
% =========================================================================
% NB
%
% =========================================================================
% *** TODO 
% ..... 
% AMPSTOINT()
%
%   consider output 'clipFlag' for instance where input current exceeds 
%   max allowable current?
%
%   likewise for any set() function in ShimCom for shim currents.
%
% ..... 
% ()
%
%
% ..... 
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
%

    
    Shims.Com.baudRate    = 57600 ;  
    Shims.Com.readTimeout = 500 ; %[units: ms] 
    
    % As specified in MXD user manual (9700043-0000) pg.16
    Shims.Com.dataBits    = 8 ;
    Shims.Com.stopBits    = 1 ;
    Shims.Com.flowControl = 'NONE' ;
    Shims.Com.parity      = 'NONE' ;
    Shims.Com.byteOrder   = 'bigEndian' ;
    
    % as specified in the manual, min delay (in seconds) 
    % btw transmission and reception of data 
    Shims.Com.txRxDelay            = 0.001 ; % [units: s]
    
    Shims.Amp.maxCurrentPerChannel = 5 ; % (absolute) [units: amps]
    Shims.Amp.maxCurrentPerBank    = 20 ; % (absolute) [units: amps]
    Shims.Amp.maxCurrentPerRail    = 10 ; % +/- [units: amps]
    
    Shims.Amp.nChannels       = 32 ;  
    Shims.Amp.nActiveChannels = 24 ;

    Shims.Dac.resolution = 16 ; % [bits]
    Shims.Dac.maxCurrent = 5 ; % (absolute) [units: amps]


end
% =========================================================================
function bitValue = ampstoint( Shims, current )
%AMPSTOINT - wrap real number in Amperes to range of DAC
% 
% bitValue = AMPSTOINT( Shims, current ) 

    MAX_DIGI = (2^(Shims.Dac.resolution-1)) - 1 ; % Digital to Analog Converter max value

    bitValue = int16( current*( MAX_DIGI/Shims.Dac.maxCurrent ) ) ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
