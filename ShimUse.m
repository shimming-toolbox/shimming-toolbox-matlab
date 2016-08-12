classdef ShimUse < ShimCom
%SHIMUSE
%
% .......
%
% Shim = ShimUse( pathToCalibrationInfo )
%
%   Shim contains fields
%
%       .Opt
%           Object of type ShimOpt
%
%       .Specs
%           Object of type ShimSpecs
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimUse
%    ShimTest 
%
% =========================================================================
% Updated::20160802::ryan.topfer@polymtl.ca
% =========================================================================
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

% Real-time
%
% 1. Calibration: GRE-scans (inspired, expired). Record pressure logs with: 
%   ProbeTracking.recordandplotpressurelog 
%
% 2. Load pressure logs. User selects begin & end points of apnea. Median extracted.
%    
%
%


% realtimeshimming...
%
%     currents = ShimOpt.pressuretocurrents( p ) % ??
%     ShimCom.setandloadshim( currents )



properties   

Opt;

end


% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimUse(  )
%SHIMUSE   

Shim.Opt = ShimOpt;

Shim.Parameters.runMode = 'isCmdLine' ; % vs. 'isGui' 


end
% =========================================================================
function [] = display( Shim, msg )
%DISPLAY

if nargin < 2 || isempty(msg)
    SystemInfo = Shim.getsysteminformation( )
    msg = SystemInfo ;
end

assert( isstr(msg), 'Given message is not a string.' ) ;

switch Shim.Parameters.runMode 
    case 'isCmdLine'
        fprintf(['\n' msg '\n']) ;
    case 'isGui'
        fprintf(['\n' 'Error: GUI not yet supported!' '\n\n']) ;
end

end
% =========================================================================
function ChannelOutputs = getallchanneloutputs( Shim )
%GETALLCHANNELSOUTPUTS
%
% ChannelOutputs = GETALLCHANNELOUTPUTS( Shim ) 
% 
% Returns struct ChannelOutputs with fields
%
%   .current
%   .voltage
%   .power
%   .dissipatedPower


channelsToBankKey = Shim.getchanneltobankkey ;

ChannelOutputs.current = zeros( 1, Shim.Specs.nActiveChannels ) ;
ChannelOutputs.voltage = zeros( 1, Shim.Specs.nActiveChannels ) ;
ChannelOutputs.power   = zeros( 1, Shim.Specs.nActiveChannels ) ;
ChannelOutputs.dissipatedPower = zeros( 1, Shim.Specs.nActiveChannels ) ; 

for iChannel = 1 : Shim.Specs.nActiveChannels 
    
    ChannelOutput = Shim.getchanneloutput( channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) ) ;

    ChannelOutputs.current(iChannel)         = ChannelOutputs.current ;
    ChannelOutputs.voltage(iChannel)         = ChannelOutputs.voltage ;
    ChannelOutputs.power(iChannel)           = ChannelOutputs.power ;
    ChannelOutputs.dissipatedPower(iChannel) = ChannelOutputs.dissipatedPower ;

end

end
% =========================================================================
function [] = setandloadallshims( Shim )
%SETANDLOADALLSHIMS


% setallshims( Shim, currents )
%   Sets all shims based on [nChannel x 1] current vector (in amps)
% 
% (e.g. for our system, currents vector has 32 (and not 24) entries.
% Currents for inactive channels should be zero.
%
% 1. convert 24-currents into 32-current (~zero-padded~ vector)
% 2. set all channel buffers using 32-component vector
% 3. issue load all channels command


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

end
% =========================================================================
% =========================================================================

end
