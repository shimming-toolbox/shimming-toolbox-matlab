classdef Config < dynamicprops 
%b0shim.Config Shim system configuration (hardware description)
% 
% `Config` object, together with a set of reference maps (basis set), forms
% the basic representation of a shim system in the Shimming Toolbox. Notably,
% to optimize the current configuration of a given array for a target b0-field
% distribution (namely, to *shim it!*) these are the key prerequisities.
%
% When called upon by other Toolbox components, the `Config` object "in action"
% behaves much like a static struct: only its *properties* are generally at
% play. As such, system-specific objects can be conveniently saved as .json
% config files to be loaded and validated by the generic interface—the `Config()`
% constructor:
%
% __CONSTRUCTOR SYNTAX__
%    
%    self = b0shim.Config( filename )
%
% Loads the shim configuration object `self` into memory using the config.json
% file `filename`, which can be a local system file or URL. (If the file is
% hosted on GitHub, note that the URL must be to the "Raw" file.)
% 
% **NOTE**: *Object properties are immutable* (fixed upon construction);
% nonetheless, system configurations can easily be edited or created, and
% optionally saved to .json file. This can be achieved by editing/creating
% the config.json file directly, or within MATLAB. e.g. by adapting an existing script an
% such as [ b0shim.coils.greg ](./coils/greg/write_config_file.md)).
%
% Briefly, call the constructor without arguments to return a default
% object and convert it to a struct (ignore the warning issued by MATLAB)
% ```
%    params = struct( b0shim.Config() );
% ```
% The fields of `params` correspond to Config properties and can be assigned
% as needed, e.g. to add and configure 2 channels:
% ```
%    params.channels(1:2)    = b0shim.parts.Channel;
%    
%    % Configure the 1st channel
%    params.channels(1).name        = "Ch. 1";
%    params.channels(1).units       = "A";
%    params.channels(1).limits      = [-5 5]; 
%    params.channels(1).positioning = "iso-fixed";
%
%    % Configure the 2nd channel
%    params.channels(2).name        = "Ch. 2";
%    params.channels(2).units       = "A";
%    params.channels(2).limits      = [-3 3]; 
%    params.channels(2).positioning = "iso-fixed";
%
%    % Assign a name to the entire system:
%    params.name = "my_shim_array";
%
%    % Etc.
% ```
%
% Once configured, the struct can be recast as a proper Config object
% by passing it to the constructor, and saved by calling `write_json` with
% the desired filename, e.g.
% ``` 
%    self = b0shim.Config( params );
%    mkdir( './+b0shim/+coils/+my_shim_array' );
%    write_json( self, './+b0shim/+coils/+my_shim_array/config.json' );
% ``` 
% 

% When reinitializing a `Config` object from an existing json file, `filename`
% can be supplied as a string or char vector to a publicly accessible URL, or
% to a local file. 
%
% __EXAMPLE__
%
% The rudimentary config.json file should resemble
%
% __ETC__
%
% 1. `Config` inherits from `dynamicprops`: New properties can optionally be
% added via `addprop` method. 
%
% 2. `dynamicprops` is a `handle` class: Beware of the distinction when copying
% a variable instance!
%
% 3. Should further customizing be desired (e.g. definining additional methods)
% `Config` can be subclassed. 
% 
% 4. Though the configuration allows a mix of different `channel.positioning` modes
% the functionality to deal with the mix does not yet exist—**TODO**. (This
% probably doesn't apply to any existing hardware anyway but could be interesting!)
%
% See also  
% [ dynamicprops ](https://www.mathworks.com/help/matlab/ref/dynamicprops-class.html)  
% [ addprop ](https://www.mathworks.com/help/matlab/ref/dynamicprops.addprop.html)  
% [ handle ](https://www.mathworks.com/help/matlab/ref/handle-class.html)  
% [ b0shim.Contents ]( ../Contents.md )  

properties (SetAccess=immutable)

    % Shim system/package name as a string-scalar.
    % 
    % In general, all the code specific to a given shim system should be fully
    % contained in a MATLAB subpackage folder. For instance, if
    % `Config.name = "my_shim"`, the config.json file should be saved
    % (e.g. via `write_json()` to the subpackage: +b0shim/+coils/+my_shim/config.json
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

    % Shim coil-element properties as a `b0shim.parts.Channel` object array
    % 
    % Each element has the 4 configurable properties `name, units, limits, positioning`
    % (e.g., `"Ch.1", "A", [-5 5], "table-fixed"). Property desciptions are provided
    % in the class documentation.
    %
    % See also  
    % [ b0shim.parts.Channel ]
    channels(:,1) b0shim.parts.Channel ;
    
    % Serial-port parameters as a b0shim.parts.Port object (applies to specific systems only)
    %  
    % `ports` only needs to be defined for shim systems that subclass
    % `b0shim.Com` to enable direct hardware control from MATLAB. For other
    % systems (e.g. scanner shims or virtual arrays) the property can be left
    % unassigned (i.e. the object will be `empty`).
    %
    % See also  
    % [ b0shim.parts.Port ]
    ports(:,1) b0shim.parts.Port ;

end

%TODO: Add version propertie(s) for saving to file 
% (would facilitate supporting 'legacy code' features)
% 2 version types?
% - version of b0shim.Config used to save the config.json file 
% a warning could be issued when a config file loaded and seen to be out of date? e.g. using git tags?
% properties (SetAccess=private)
% end
% - version of the system configuration? (e.g. our "greg" coil should have 2
% versions for the revB, revC boards) could just be expressed in the filename.

properties( Dependent, Hidden )
    % Version of b0shim.Config as a string-scalar 
    % TODO proper implementation (i.e. no hardcode, e.g. via git tags?)
    softwareVersion = '0000';
end

% =========================================================================
% =========================================================================
methods  
% =========================================================================
function self = Config( varargin )
    
    %% Check inputs
    narginchk(0, 1);
    if nargin == 0 % Return default/template object
        return;
    elseif ischar(varargin{1}) || isstring(varargin{1})
        params = read_json( varargin{1} );
    elseif isstruct( varargin{1} )
        params = varargin{1};
    else
        error( [mfilename ' constructor requires a single input:\n ' ...
            'The path to a (json) config file, or a struct of configuration parameters'] )
    end

    %% Assign properties from input struct 
        
    fieldNames = fieldnames( params );
    metaConfig  = meta.class.fromName( 'b0shim.Config' );
    propNames  = {metaConfig.PropertyList.Name} ;

    for iField = 1 : numel( fieldNames )
        
        % value to be assigned
        value = getfield( params, fieldNames{iField} );
       
        % look for match (case-insensitive)
        iProp = strcmpi( fieldNames{iField}, propNames );

        if any(iProp) % assign a standard property 
            self.(propNames{iProp}) = value ;
        else % assign a custom parameter
            addprop(self, fieldNames{iField});
            self.(fieldNames{iField}) = value; 
        end
    end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods (Hidden)
% =========================================================================
function [self] = loadobj( params )
%LOADOBJ Special load process for .mat files
%
% As the class constructor generally requires a struct/json input argument,
% the procedure for loading from/saving to *.mat* files needs to be customized.
% (Note, this is separate from the .json load/save procedure.) 
%  
% [ See ](https://www.mathworks.com/help/matlab/matlab_oop/passing-arguments-to-constructors-during-load.html)  
%
% See also
% saveobj
% [ loadobj ](https://www.mathworks.com/help/matlab/ref/loadobj.html)

    self = b0shim.Config( params );

end
% =========================================================================
function [params] = saveobj( self )
%SAVEOBJ Special save process for .mat files
%
% See also  
% loadobj

    warning( 'OFF', 'MATLAB:structOnObject' );
    params = struct( self ); 
    warning( 'ON', 'MATLAB:structOnObject' );

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
