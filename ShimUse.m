classdef ShimUse < ShimCom
%SHIMUSE
% .......
%
% Shim = ShimUse( pathToCalibrationInfo )
%
%   Shim contains fields
%
%       .Opt
%           Object of type ShimOpt
%
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%     ProbeTracking
%     ShimCom
%     ShimOpt
%     ShimSpecs
%     ShimUse
%
% =========================================================================
% Updated::ryan.topfer@polymtl.ca::Tue 28 Jun 2016 18:23:55 EDT
% =========================================================================







realtimeshimming...

    currents = ShimOpt.pressuretocurrents( p ) % ??
    ShimCom.setandloadshim( currents )

% 
% Shim Use (and: Shim Use[r settings])
%   (Highest level)
% 
% Shims = ShimUse(  )
%
%   Shims contains fields
%
%           
%       .runMode    
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
%
% ..... 
% ()
%
%
% ..... 
% =========================================================================

properties   
    runMode; % '

end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

methods
% =========================================================================
function Shim = ShimUse(  )
%SHIMUSE   

    % Shims.Parameters.isRunningGui     = false ;
    
end
% =========================================================================
end
% =========================================================================
% =========================================================================
methods(Static) 
% =========================================================================
function [systemResponses, meanings] = definesystemresponses( )
%DEFINESYSTEMRESPONSES
%
% [systemResponses, meanings] = DEFINESYSTEMRESPONSES( )
%
%   systemResponses is a string array of HEX messages
% 
%   meanings is a cell array where each entry gives a brief definition of the
%   corresponding entry in systemResponses.
% 
% -------------------------------------------------------------------------
% According to RRI HEX specification protocol, system responses from the 
% MXD/DSU can be one of the following:  
%
% i.   0x00: ACK , command has been received and processed
% ii.  0x01: Unimplemented command
% iii. 0x02: A command parameter is out of range
% iv.  0x03: An error is present that prevents execution of the command
% v.   0x08: Invalid command
% vi.  0xFF: NACK, unrecognized command or bad CRC

systemResponse = ['0x00'; '0x01'; '0x02'; '0x03'; '0x08'; '0xFF'] ;

meanings    = cell{length(systemResponse), 1} ;
meanings{1} = 'ACK , command has been received and processed' ;
meanings{2} = 'Unimplemented command' ;
meanings{3} = 'A command parameters is out of range' ;
meanings{4} = 'An error is present that prevents command execution' ;
meanings{5} = 'Invalid command' ;
meanings{6} = 'NACK, unrecognized command or bad CRC' ;

end
% =========================================================================
function [] = listcommands( option )
%LISTCOMMANDS
%
% [] = LISTCOMMANDS( )
% [] = LISTCOMMANDS( 'usage' )
%
% option
%       if omited [default], command names alone are listed
% 
% TODO:
% if == 'usage', command names are listed along with their usage
%
%%

Cmd = ShimCom.definecommands( ) ;

fprintf(['\n MXD Commands:\n\n'])
commands = fieldnames( Cmd.Mxd ) ;

for iCmd = 1 : length(commands)
    fprintf([ commands{iCmd} '\n'])
end

fprintf(['\n\n'])

% fprintf(['\n DSU Commands:\n'])
% commands = fieldnames( Cmd.Dsu ) ;
%
% for iCmd = 1 : length(commands)
%     fprintf(['\n' commands{iCmd}])
% end
    

end
% =========================================================================
% =========================================================================
end
