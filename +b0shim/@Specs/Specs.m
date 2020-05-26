classdef Specs < dynamicprops
%b0shim.Specs Shim System Specifications
%    
%    Specs = b0shim.Specs( )
% 
% Returns a b0shim.Specs object.
% Specs defines all the relevant hardware specifications of a shim system.
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

properties  

    % Shim system ID: 1 string scalar. Used to tabulate `Opt`imization results
    systemName(1,1) {valid.mustBeA(systemName, "string")} = "MyShim" ;

    % Channel IDs: 1-string/channel. Used to tabulate `Opt`imization results
    channelNames(:,1) {valid.mustBeA(channelNames, "string")} = ["Ch1"; "Ch2"; "Ch3"] ;

    % Shim units: 1-string/channel. Used to tabulate `Opt`imization results 
    channelUnits(:,1) {valid.mustBeA(channelUnits, "string")} = ["[A]"; "[A]"; "[A]"] ;

    maxCurrentPerChannel(:,1) {mustBeNumeric} = [5; 5; 5]

end

properties( Dependent, SetAccess='private' )

    nActiveChannels {mustBeScalar, mustBeIntegerOrEmpty} ;

end

properties
Amp=[]; % relating to amplification
Com=[]; % relating to communication (e.g. RS-232)
Dac=[]; % relating to digital-to-analog conversion 
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function self = Specs( jsonPath )
%    
%    specs = b0shim.Specs.load_json( filename ) ;
%    specs = b0shim.Specs.load_json( url ) ;
%
% Loads and decodes a json file to initialize a struct then used to form `specs` a
% b0shim.Specs object 

% ------------
%% Check inputs
    if nargin == 0
        return; 
    end

    if [ nargin ~= 1 ] || ...
          ~[ ischar(jsonPath) | isstring(jsonPath) ]
        error('Specs constructor accepts a single input: a filepath or URL to a json file') ;
    end

    specs = load_json( jsonPath ) 

% % ------- 
% % COM 
% self.Com.baudRate    = 57600 ;  
% self.Com.readTimeout = 500 ; % [units: ms] 
%
% self.Com.dataBits    = 8 ;
% self.Com.stopBits    = 1 ;
% self.Com.flowControl = 'NONE' ;
% self.Com.parity      = 'NONE' ;
% self.Com.byteOrder   = 'bigEndian' ;
%
% % min delay (in seconds) between transmission and reception of data.
% self.Com.txRxDelay   = 0.001 ; % [units: s]
%
% % ------- 
% % AMP 
% self.Amp.maxCurrentPerChannel = 5 ; % (absolute) [units: amps]
% self.Amp.maxCurrentPerBank    = 20 ; % (absolute) [units: amps]
% % self.Amp.maxCurrentPerRail    = 10 ; % +/- [units: amps]
%
% self.Amp.nChannels       = 32 ;  
% self.Amp.nActiveChannels = 24 ;
%
% % ------- 
% % DAC
% self.Dac.resolution = 16 ; % [bits]
% self.Dac.maxCurrent = 5 ; % (absolute) [units: amps]

% ------------
%% Local functions 
    
    function [specs] = load_json( jsonPath )
    %LOAD_JSON Return `specs` struct from json file (system path or URL)

        if isfile( jsonPath )
            jsonTxt = fileread( jsonPath );
        else
            jsonTxt = webread( jsonPath, weboptions( 'ContentType', 'text') );
        end
            
        specs = jsondecode( jsonTxt );

    end

end
% =========================================================================
function [nActiveChannels] = get.nActiveChannels( self )

    nActiveChannels = length( self.channelNames ) ;

end
% =========================================================================
function [] = set.maxCurrentPerChannel( self, maxCurrentPerChannel )

    switch length(maxCurrentPerChannel) 
        case self.nActiveChannels
            self.maxCurrentPerChannel = maxCurrentPerChannel ;
        case 1 
            self.maxCurrentPerChannel = repmat( maxCurrentPerChannel, [self.nActiveChannels 1] ) ;
        otherwise
            error( 'b0shim:Specs:nCurrents_neq_nChannels', ...
                ['maxCurrentPerChannel can be assigned either as a scalar (constant across channels),\n'
                 ' or as a vector with an entry for each of the `nActiveChannels`.\n'] ) ;
    end
end
% =========================================================================
function [] = set.channelUnits( self, channelUnits )

    switch length( channelUnits ) 
        case self.nActiveChannels
            self.channelUnits = channelUnits ;
        case 1 
            self.channelUnits = repmat( channelUnits, [self.nActiveChannels 1] ) ;
        otherwise
            error( 'b0shim:Specs:nChannelUnits_neq_nChannels', ...
                ['channelUnits can be assigned either as a scalar (constant across channels),\n'
                 ' or as a vector with an entry for each of the `nActiveChannels`.\n'] ) ;
    end
end
% =========================================================================
function [] = save_json( self, filename )
%SAVE_JSON Writes a `Specs` instance to json file
%    
%    save_json( self, filename )
% 
% Encodes `Specs`-type object `self` as json and writes to the path string
% `filename`.
%
% __NOTE__
% The function uses 

    narginchk(2,2);
    filename = string(filename) ; % typecast in case input as char
    assert( length(filename)==1, ...
        'save_json requires 2 inputs: The b0shim.Specs object and the ouput filename string') ;

    json = jsonencode( self ) ;
    
    fid = fopen( filename ) ;
    fwrite( fid, json ) ; 
    fclose(fid) ;
    
end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
