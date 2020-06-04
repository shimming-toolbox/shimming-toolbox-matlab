classdef Specs < dynamicprops
%b0shim.Specs Shim system specifications (hardware description)
% 
% A `Specs` object, together with a set of reference maps (basis set), forms
% the basic representation of a shim system in the Shimming Toolbox. Notably,
% to optimize the current configuration of a given array for a target b0-field
% distribution (namely, to *shim it!*) these are the key prerequisities.
%
% The `Specs` object "in action" (when called upon by other Toolbox components)
% behaves much like a static struct: only its properties (`systemName`,
% `maxCurrentPerChannel`, etc.) are generally at play. As such, system-specific
% objects can be conveniently saved as .json files, to be loaded and validated
% by the generic interface—the `Specs()` constructor:
%
% __CONSTRUCTOR SYNTAX__
%    
%    self = b0shim.Specs( )
%    self = b0shim.Specs( pathToJson )
%     
% The void call to `b0shim.Specs()` returns the default object `self`. 
% This is effectively a "wizard" mode to help define a new shim system (an
% alternative to writing the json from scratch) and should only need to be
% performed once: Properties can be assigned prior to writing to disk via  
%      
%     self.save_json( filename )  
%
% When reinitializing a `Specs` object from file, `pathToJson` can be supplied
% as a string or char vector to a publicly accessible URL, or to a local file. 
%
% __ETC__
%
% 1. `Specs` inherits from `dynamicprops`: New properties can optionally be
% added via `addprop` method. 
%
% 2. `dynamicprops` is a `handle` class: Beware of the distinction when copying
% a variable instance!
%
% 3. Should further customizing be desired (e.g. definining additional methods)
% `Specs` can be subclassed. 
%
% __ETC__
%
% See also 
% dynamicprops https://www.mathworks.com/help/matlab/ref/dynamicprops-class.html  
% addprop https://www.mathworks.com/help/matlab/ref/dynamicprops.addprop.html  
% handle https://www.mathworks.com/help/matlab/ref/handle-class.html

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

    % Shim system/package name: A single string-scalar.
    % 
    % In general, all the code specific to a given shim system should be fully
    % contained in a MATLAB subpackage folder. For example, if
    % `specs.name = "my_shim"`, then the specs.json file should be saved
    % (e.g. via `specs.save_json`) to the subpackage: +b0shim/+coils/+my_shim/specs.json
    %
    % `systemName` defines the column header when tabulating optimization results 
    % in `b0shim.Opt.optimizeshimcurrents`
    %
    % e.g.
    %
    %|    my_shim    | Original | Optimal | Update | Realtime | Relative_Power |
    %| ------------- | -------- |   ---   |  ---   |    ---   | -------------- |
    %|   MC-Ch.1 [A] |    0     |  1.23   |  1.23  |   0.01   |     0.05       |
    %|   Gz [mT/m]   |    0     |  0.1    |  0.1   |   0.00   |     0.00       |
    % 
    % __ETC__
    %
    % See also  
    % channelNames  
    % channelUnits  
    % save_json  
    name(1,1) string {valid.mustBeA(name, ["string" "char"])} = "my_shim" ;

    % Channel IDs: 1 string/channel. 
    %
    % Used to tabulate `Opt`imization results
    %
    % `channelNames` defines the row names when tabulating optimization results 
    % in `b0shim.Opt.optimizeshimcurrents`
    % 
    % See also
    % channelUnits
    channelNames(:,1) {valid.mustBeA(channelNames, "string")} = ["Ch1"; "Ch2"; "Ch3"] ;

    % Shim units: 1 string/channel.  
    %
    % Along with `channelNames`, `channelUnits` is used when tabulating
    % optimization results in `b0shim.Opt.optimizeshimcurrents`
    % 
    channelUnits(:,1) {valid.mustBeA(channelUnits, "string")} = ["[A]"; "[A]"; "[A]"] ;

    maxCurrentPerChannel(:,1) {mustBeNumeric} = [5; 5; 5]

end

properties( Dependent, SetAccess='private' )

    nActiveChannels {mustBeScalar, mustBeIntegerOrEmpty} ;

end

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
    narginchk(0, 1);
    
    %% Check inputs
    if nargin == 0
        return; 
    end

    assert( [ischar(jsonPath) | isstring(jsonPath)] , ...
        'Specs constructor accepts a single input: a filepath or URL to a json file') ;

    specs = load_json( jsonPath ) 

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
            self.maxCurrentPerChannel = repmat( maxCurrentPerChannel, [self.nActiveChannels 1] );
        otherwise
            error( 'b0shim:Specs:nCurrents_neq_nChannels', ...
                ['maxCurrentPerChannel can be assigned either as a scalar (constant across channels),\n'
                 ' or as a vector with an entry for each of the `nActiveChannels`.\n'] );
    end
end
% =========================================================================
function [] = set.channelUnits( self, channelUnits )

    switch length( channelUnits ) 
        case self.nActiveChannels
            self.channelUnits = channelUnits;
        case 1 
            self.channelUnits = repmat( channelUnits, [self.nActiveChannels 1] );
        otherwise
            error( 'b0shim:Specs:nChannelUnits_neq_nChannels', ...
                ['channelUnits can be assigned either as a scalar (constant across channels),\n'
                 ' or as a vector with an entry for each of the `nActiveChannels`.\n'] );
    end
end
% =========================================================================
function [] = save_json( self, filename, Options )
%SAVE_JSON Writes a `Specs` instance to json file
%    
%    specs.save_json( filename )
%    specs.save_json( filename, "isOverwriting", true )
% 
% Encodes the `specs` object as json and writes to the path string `filename`.
%
% To overwrite an existing file, use the `(..., "isOverwriting", true)`
% argument-pair. 
%
% __ETC__
%
% - `save_json` wraps to Matlab's `jsonencode`, which outputs content in "compact"
% form. If desired, many convenient tools exist to prettify (expand) the .json
% file once written, such as:  
%   1. [webbrowser](https://jsonformatter.org/json-pretty-print)  
%   2. python CLI: `python -m json.tool filename`
%   3. [vscode](https://marketplace.visualstudio.com/items?itemName=vthiery.prettify-selected-json)
%
% See also  
% jsonencode    https://www.mathworks.com/help/matlab/ref/jsonencode.html
    arguments
        self b0shim.Specs;
        filename string;
        Options.isOverwriting {mustBeBoolean} = false;
    end
    
    assert( length(filename)==1, ['save_json requires 2 inputs: ' ...
        'A b0shim.Specs instance and a string specifying the ouput file'] );

    if ~contains(filename, '.') % Add extension
        filename = strcat( filename, '.json') ;
    end

    assert( [~isfile(filename) | Options.isOverwriting], ['File already exists.' ...
         'Use the name-value pair (...,"isOverwriting", true) to force overwrite.'])

    % [~, folderName] = fileparts( fileparts( filename ) );
    %
    % if ~startsWith(folderName, '+')
    %     warning('specs.json should typically be saved to a coil-specific package directory...');
    % end

    json = jsonencode( self ) ;
    
    fid = fopen( filename, 'w' ) ;
    fwrite( fid, json ) ; 
    fclose(fid) ;
    
end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
