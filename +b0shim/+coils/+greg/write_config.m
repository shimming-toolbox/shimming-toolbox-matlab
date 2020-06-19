function [ filename ] = write_config( filename )
%WRITE_CONFIG Writes system configuration file for the 8 ch. AC/DC neck coil
%     
%     filename = b0shim.coils.greg.write_config();
%     b0shim.coils.greg.write_config( filename );
%
% Defines the system configuration and saves it to the coil's package folder under
% ```
%   ./+b0shim/+coils/+greg/config.json
% ```
% b0
%
% and saves the system configuration for
% the 3T 8-channel AC/DC neck-coil 
% 
% See also
% [ b0shim.coils.greg.Contents ]( ./Contents.md )
    arguments
        filename(1,:) char = [mfilename('fullpath') filesep 'config_' datestr(date,'yymmdd') ];
    end
warning('Untested commit - may contain bugs')

%% 1. Create an instance of the default configuration object
config = b0shim.Config;

%% 2. *As the object's properties are immutable*, convert to struct:  
% NOTE: Conversion elicits a warning. Not a concern. Mute it momentarily)
warning('OFF', msgId);
config = struct( config );
warning('ON', msgId);

%% 3. Assign fields according to system specs 

% System/package name
config.name = "greg"; 

%% Channel configuration

% Add 8 channels
config.channels(1:8) = b0shim.parts.Channel ; 

% Assign channel properties 
for iCh = 1 : numel(config.channels)
    config.channels(iCh).name  = strcat( "Ch. ", num2str(iCh) );
    config.channels(iCh).units = "A"; % Amperes
    config.limits(iCh)         = [-2.5 2.5]; % [min. max.] in Amperes 
    config.positioning         = "table-fixed";  % shim position is table-dependent
end

%% Configuration of serial port parameters
config.ports(1) = b0shim.parts.Port();

% NOTE: These parameters are actually set *in the source code of the Teensy
% microcontroller.*
% TODO: Rather than hardcode the values in two places, some consistency ought
% to be enforced between the C++ and MATLAB sources, e.g. by MATLAB parsing the
% .cpp/.hpp files, or the config.json values being transfered over serial, or
% some other alternative. (Same goes for the Arduino-based respiratory bellows...) 

config.ports(1).BaudRate    = 115200;
config.ports(1).DataBits    = 8;
config.ports(1).StopBits    = 1;
config.ports(1).Parity      = "none"
config.ports(1).FlowControl = "none";
config.ports(1).ByteOrder   = "big-endian";

%% Configure remaining properties related to the microcontroller

%% 4. With everything configured, convert back to object
config = b0shim.Config( config );

%% 5. Save as .json
write_json( config, filename )
 
end
