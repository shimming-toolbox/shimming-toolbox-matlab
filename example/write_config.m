%EG_CONFIGURE_SYSTEM Create example system configuration; save as config.json 
%
% __TODO__ *This example is a WIP* 
% 
% __SUMMARY__
%
% Briefly, the basic procedure for configuring a new system is the following:  
%
% 1. Create an instance of the default configuration object:   
% `params=b0shim.Config;`
%
% 2. *As the object's properties are immutable*, convert it to a struct:  
% `params=struct(params);`
%
% 3. Assign fields according to the given system, e.g.  
% ```
%    % Create 3 coil-elements
%    params.channels(1:3)           = b0shim.parts.Channel; 
%    
%    % Configure the 1st channel
%    params.channels(1).name        = "Ch. 1";
%    params.channels(1).units       = "A";
%    params.channels(1).limits      = [-5 5]; 
%    params.channels(1).positioning = "iso-fixed";
%
%    % Configure remaining channels, etc.
% ```
%
% 4. Once all the assignments have been made, convert back to the original type:  
% `params=b0shim.Config(params);`
%
% 5. Optionally, save as .json:a
% `write_json( params, "config.json" );`
%
% __ETC__
% 
% Additional example scripts for existing coils should be found in their
% respective subpackage folders, e.g.: 
%
% See also  
% [ b0shim.coils.greg.write_config ]
