classdef Specs < dynamicprops
%b0shim.Specs Shim system specifications (hardware description)
%
% A `Specs` object, together with a set of reference maps (basis set), forms
% the basic representation of a shim system in the Shimming Toolbox. Notably,
% to optimize the current configuration of a given array for a target b0-field
% distribution (namely, to *shim it!*) these are the key prerequisities.
%
% When called upon by other Toolbox components, the `Specs` object "in action"
% behaves much like a static struct: only its *properties* are generally at
% play. As such, system-specific objects can be conveniently saved as .json
% files, to be loaded and validated by the generic interface—the `Specs()`
% constructor:
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
% __EXAMPLE__
%
% A rudimentary config.json file:
% (Note that functionality to mix different elements from different arrays-fixed & iso-fixed does not exist yet)
% {
%     "name": "my_shim",
%     "channels": [
%         {
%             "name": "Ch. 1",
%             "units": "A",
%             "limits": [
%                 -5,
%                 5
%             ],
%             "positioning": "table-fixed"
%         },
%         {
%             "name": "Ch. 2: X-Gradient",
%             "units": "mT/m",
%             "limits": [
%                 -100,
%                 100
%             ],
%             "positioning": "iso-fixed"
%         },
%     ],
%     "com": [],
%     "pathToReferenceMaps": "",
%     "this": []
% }
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
% dynamicprops https://www.mathworks.com/help/matlab/ref/dynamicprops-class.html  
% addprop https://www.mathworks.com/help/matlab/ref/dynamicprops.addprop.html  
% handle https://www.mathworks.com/help/matlab/ref/handle-class.html

properties  

    % Shim system/package name as a string-scalar.
    % 
    % In general, all the code specific to a given shim system should be fully
    % contained in a MATLAB subpackage folder. For instance, if
    % `specs.name = "my_shim"`, the specs.json file should be saved
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
    % channels  
    % save_json    
    name(1,1) string = "my_shim" ;
    
    % URL or local path string to the json config file 
    filename(1,1) string = "";
    
    % URL or local file path string to the shim basis maps (NIfTI images?TODO/TBD) 
    pathToReferenceMaps(:,1) string = "";

    % b0shim.parts.Channel object-array 
    %
    % See also  
    % b0shim.parts.Channel
    channels(:,1) b0shim.parts.Channel ;
    
    % b0shim.parts.Port object (Optional)
    %  
    % `ports` only needs to be defined for shim systems that subclass `b0shim.Com`
    % to enable direct hardware control from MATLAB.
    %
    % For other systems (e.g. scanner shims or virtual arrays) the property can
    % be left empty.
    %
    % See also
    % b0shim.parts.Port
    ports(:,1) b0shim.parts.Port ;

end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function self = Specs( jsonPath )
    
    narginchk(0, 1);
    
    %% Check inputs
    if nargin == 0 % Return default/template object
        return;
    % else % Decode json file to initialize a struct; convert to a b0shim.Specs object 
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
        Options.isOverwriting {mustBeScalarOrEmpty,valid.mustBeBoolean} = false;
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
function [isValid] = validate( self )
%VALIDATE Return true if Specs configuration meets requirements
    isValid=false;
end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
