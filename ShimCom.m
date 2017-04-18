classdef (Abstract) ShimCom < matlab.mixin.SetGet
%SHIMCOM - Shim Communication
%
% .......
% 
% Description
%
%   SHIMCOM is responsible for all direct communication with the shim system
%   microcontoller (e.g. setting/querying shim currents). 
%   Declaration of a ShimCom object immediately opens a serial (Com) port.
% 
% .......
%   
% Usage
%
%   Shims = ShimCom(  )
%
%   Shims contains fields
%
%       .Cmd
%           
%       .ComPort    
%
%       .Data 
%           .output 
%           .input  
%
%       .Params
%
%       .Specs
%
% =========================================================================
% Notes
%
%   Based off RRI Hex Protocol Specification guide:
%   9700052-0000 HexProtocolSpecification_REV-G
%
%   For a primer on RS-232 communication in Matlab see
%   http://www.mathworks.com/help/matlab/matlab_external/overview-of-the-serial-port.html
%
%   Part of series of classes pertaining to shimming:
%
%    Tracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%    
%    ShimCom is an Abstract class.
%
% =========================================================================
% Updated::20170204::ryan.topfer@polymtl.ca
% =========================================================================



% *** TODO 
% 
% ..... 
% SPLITINT( ) 
%  
%   Documentation for how two's complement works.
% ..... 
% When converting to DAC counts (e.g ampstodac())
%   consider output 'clipFlag' for instance where input current exceeds 
%   max allowable current?
%     
% =========================================================================

properties   
    Cmd;
    ComPort;
    Data;
    Params;
    Specs;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimCom(  )
%SHIMCOM - Shim Communication

% For all concrete ShimCom subclasses, these fields need to be initialized 
% with the appropriate (concrete) subclass 
%
% e.g. for ShimComRri(), 
Shim.Specs   = [] ; % Shim.Specs   = ShimSpecsRri( ) ;
Shim.ComPort = [] ; % Shim.ComPort = ShimCom.initializecomport( Shim.Specs ) ;
Shim.Cmd     = [] ; % Shim.Cmd     = ShimComRri.getcommands( ) ;


Shim.Data = [] ;
Shim.Data.output = uint8(0) ;
Shim.Data.input  = uint8(0) ;

Shim.Params = [] ; 
Shim.Params.nBytesToRead     = [] ; % depends on cmd sent to system
Shim.Params.nSendAttemptsMax = 5; % number of communication attempts before error


end
% =========================================================================
function [] = delete( Shim )
%DELETE  (custom helper function)
% 
% DELETE( Shim )
% 
% Destructor. Calls Shim.deletecomport( ) 

Shim.deletecomport();
clear Shim ;

end
% =========================================================================
function [Shim] = deletecomport( Shim )
%DELETECOMPORT  (custom helper function)
% 
% Shim = DELETECOMPORT( Shim )
% 
% Correct way to delete and clear the serial port object

% check existence?
if myisfield( Shim, 'ComPort' ) && ~isempty( Shim.ComPort ) 
    fclose( Shim.ComPort ) ;
    delete( Shim.ComPort ) ;
    clear Shim.ComPort  ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true)

function [ComPort] = initializecomport( Specs )
% Initialize (RS-232) Communication Port 
%
% if ismac
%     portName = '/dev/tty.usbserial' ;
% elseif isunix
%     portName = '/dev/ttyS0' ;
% elseif ispc
%     portName = 'COM4' ;
% else
%     error('What kind of computer is this!?')
% end

% -------------------------------------------------------------------------
% Serial Port 

if ismac
    portName = '/dev/tty.usbserial' ;
elseif isunix
    portName = '/dev/ttyS0' ;
elseif ispc
    portName = 'COM4' ;
else
    error('What kind of computer is this!?')
end

ExistingPorts = instrfind( ) ;

for iPort = 1 : length( ExistingPorts )
    if strcmp( ExistingPorts( iPort ).Port, 'portName' ) ...
    && strcmp( ExistingPorts( iPort ).Status, 'open' )

        fclose( ExistingPorts( iPort ).Port )
    end
end

ComPort = serial( portName,...
    'BaudRate', Specs.Com.baudRate,...
    'DataBits', Specs.Com.dataBits,...
    'StopBits', Specs.Com.stopBits,...
    'FlowControl', Specs.Com.flowControl,...
    'Parity', Specs.Com.parity,...
    'ByteOrder', Specs.Com.byteOrder ) ;   

end
% =========================================================================
function [lowByte, highByte] = splitint( z )
% SPLITINT
%
% [lowByte, highByte] = SPLITINT( int16 z ) 
%
% if z is positive
%   z = (2^8)*highByte + lowByte ;        
%
% else if z is in two's complement
%   ...
 
    if isa(z, 'uint16' ) || z > 0 
        lowByte  = uint8( bitand( uint16(z), uint16(255) ) ) ;
        % bitshift by -8 discards the lowest 8 bits
        highByte = uint8( bitshift( z, -8 ) ) ;
    elseif isa(z, 'int16')
        % invert bits
        % [0, 65535] is the range of uint16.
        % 2^16 = 65536 = 1 0000 0000 0000 0000 
        % Subtraction of 2^16 from z works to flip the bits
        % and add 1  
        z = abs(int32(abs(z)) - int32(65536)) ;
        
        lowByte  = uint8( bitand( z, int32(255) ))  ;

        % bitshift by -8 discards the lowest 8 bits
        highByte = uint8( bitshift( z, -8 ) ) ;
    else
        error('Input must be int16 or uint16') ;
    end 

end
% =========================================================================
function z = mergeints( highByte, lowByte, isSigned )
% MERGEINTS
%
% int16 z = MERGEINTS( uint8 highByte, uint8 lowByte, isSigned )


assert( nargin == 3 ) 
assert( isa( lowByte, 'uint8' ) && isa( highByte, 'uint8' ) ) 

% bitshift by 8 equivalent to multiplying by 2^8
z = bitshift( uint16(highByte), 8 ) ; 
z = z + uint16( lowByte ) ;

if ~isSigned
    return

% 127 in binary: 0111 1111
% leftmost 1 indicates negative in two's complement scheme
elseif isSigned && highByte > 127

    % [0, 65535] is the range of uint16.
    % 2^16 = 65536 = 1 0000 0000 0000 0000 
    % Subtraction of 2^16 from z works to flip the bits
    % and add 1 once cast as int16? 
    
    z = int16( int32(z) - 65536 ) ;    
end

end
% =========================================================================
function z = convertfromtwoscomplement( zUInt ) 
% CONVERTFROMTWOSCOMPLEMENT
%   
% int16 z = CONVERTFROMTWOSCOMPLEMENT( uint8 z ) 

assert( isa(zUInt, 'uint8') )

if zUInt < 127
    z = int16( zUInt ) ;
else
    z = -int16( bitxor( zUInt, uint8(255) ) + 1 ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Abstract)
% =========================================================================
[Cmd] = getcommands( )
%GETCOMMANDS     (STATIC METHOD!)
%
% Cmd = GETCOMMANDS( )
%
% Returns struct of all available string commands,
% 
% e.g., For the RRI HEX implementation...
%
%   Cmd.Mxd.getSystemHeartbeat          = '09' ;
%   Cmd.Mxd.getChannelOutput            = '47' ;
%   Cmd.Mxd.setAndLoadShimByBankChannel = '44' ;
%   Cmd.Mxd.setAndLoadShimByChannel     = '54' ;
%   etc...
% =========================================================================
[isAckReceived] = getsystemheartbeat( Shim ) ;
%GETSYSTEMHEARTBEAT
%
% Queries shim controller and returns true if responsive
% =========================================================================
[] = setandloadshim( Shim, varargin )  ;
%SETANDLOADSHIM     (MXD cmd 0x44 || 0x54)
%
% Set shim current (in amps) for single channel 
% 
% [] = SETANDLOADSHIM( Shim, channelIndexGlobal, current ) 
% [] = SETANDLOADSHIM( Shim, bankIndex, channelIndexByBank, current ) 
% =========================================================================
[] = setandloadallshims( Shim, currents )
%SETANDLOADALLSHIMS     (custom cmd)
% 
% Sets all shim buffers (MXD cmd 0x22) and loads the settings (MXD cmd 0x23).
% =========================================================================
[] = resetallshims( Shim ) 
%RESETALLSHIMS  
%
% Reset all shims to 0 A.
%
% [] = RESETALLSHIMS( Shim )
% =========================================================================
[ChannelOutput] = getchanneloutput( Shim, iBank, iChannel ) 
%GETCHANNELOUTPUT   (MXD cmd 0x47)
% 
% ChannelOutput = getchanneloutput( Shim, bankIndex, channelIndex ) 
%
% ChannelOutput contains fields
%   .current [in amps]
%   .voltage [in volts]
%   .power [in Watts]
%   .disspitatedPower [in Watts]
% =========================================================================
[ChannelOutputs] = getallchanneloutputs( Shim )
%GETALLCHANNELSOUTPUTS      (custom cmd)
%
% ChannelOutputs = GETALLCHANNELOUTPUTS( Shim ) 
% 
% ChannelOutputs has fields
%
%   .current [amperes]
%   .voltage [volts]
%   .power [watts]
%   .dissipatedPower [watts]

% =========================================================================
% =========================================================================
end
    
end
% =========================================================================
% =========================================================================
