classdef ShimComAcdc < ShimCom 
%SHIMCOMACDC - Shim Communication for AC/DC neck coil 
%
% .......
%   
% Usage
%
%   Shims = ShimComAcdc(  )
%
%   Shims contains fields
%
%       .Cmd
%           
%       .ComPort    
%
%       .Data 
%
%       .Params
%
% =========================================================================
% Notes
%
%   Part of series of classes pertaining to shimming:
%
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%     
%    ShimComAcdc is a ShimCom subclass.
%
% =========================================================================
% Updated::20171020::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimComAcdc( Specs )
%SHIMCOM - Shim Communication

if nargin < 2 || isempty( Specs ) 
    Shim.Specs = ShimSpecsAcdc( );
end

Shim.ComPort = ShimComAcdc.initializecomport( Shim.Specs ) ;
Shim.Cmd     = ShimComAcdc.getcommands( ) ;

Shim.Data.output = uint8(0) ;
Shim.Data.input  = uint8(0) ;

Shim.Params.nBytesToRead     = [] ; % depends on cmd sent to system
Shim.Params.nSendAttemptsMax = 5; % # communication attempts before error

end
% =========================================================================
function [isAckReceived] = getsystemheartbeat( Shim ) ;
%GETSYSTEMHEARTBEAT
%
% Queries shim controller and returns TRUE if responsive
end
% =========================================================================
function [] = setandloadshim( Shim, voltages )  ;
%SETANDLOADSHIM    
%
% Set shim voltage for single channel 
% 
% [] = SETANDLOADSHIM( Shim, channelIndexGlobal, voltage ) 

% output message 
Shim.Data.output    = Shim.Cmd.setAndLoadShim ; 

end
% =========================================================================
function [] = setandloadallshims( Shim, voltages )
%SETANDLOADALLSHIMS     
% 

for iCh = 1 : Shim.Specs.nActiveChannels 
    Shim.setandloadshim( iCh, voltages(iCh) ) ;
end

end
% =========================================================================
function [] = resetallshims( Shim ) 
%RESETALLSHIMS  
%
% Reset all shims to 0 A.
%
% [] = RESETALLSHIMS( Shim )

Shim.Params.nBytesToRead = 5 ; % (ACK only)

% output message 
Shim.Data.output = Shim.Cmd.resetAllShims ; 

% Shim.Data.output(2) = 5 ; 
% Shim.Data.output(3) = hex2dec( Shim.Cmd.Mxd.setLoadAllShim ) ;

[Shim, isSendOk] = Shim.sendcmd() ;

% if isSendOk 
%     [isAckReceived, systemResponse] = Shim.isackreceived( ) ;
%     if ~isAckReceived
%         disp(systemResponse)
%     end
% end

end
% =========================================================================
function [ChannelOutput] = getchanneloutput( Shim, iBank, iChannel ) 
%GETCHANNELOUTPUT   
% 
% ChannelOutput = getchanneloutput( Shim, bankIndex, channelIndex ) 
%
% ChannelOutput contains fields
%   .current [in amps]
%   .voltage [in volts]
%   .power [in Watts]
%   .disspitatedPower [in Watts]

Shim.Params.nBytesToRead = 5 ; % (ACK only)

% output message 
Shim.Data.output = Shim.Cmd.getChannelOutput ; 

ChannelOutput.voltage = 

end
% =========================================================================
function [ChannelOutputs] = getallchanneloutputs( Shim )
%GETALLCHANNELSOUTPUTS      
%
% ChannelOutputs = GETALLCHANNELOUTPUTS( Shim ) 
% 
% ChannelOutputs has fields
%
%   .current [amperes]
%   .voltage [volts]
%   .power [watts]
%   .dissipatedPower [watts]

end
% =========================================================================
function [Shims, isSendOk]= sendcmd( Shim )
%SENDCMD 
%   
%   Transmits a command from client to shim microcontroller 
%

% isSendOk     = false ;
% iSendAttempt = 0 ;
%
% while ~isSendOk && (iSendAttempt < Shims.Params.nSendAttemptsMax) 

isSendOk  = false ;
% isMsgRead = false ;

if strcmp( Shim.ComPort.Status, 'closed' ) ;
    fopen( Shim.ComPort ) ;
end 
    
fprintf( Shim.ComPort, Shim.Data.output ) ;
pause( Shim.Specs.Com.txRxDelay ) ;
    
    % [Shim.Data.input, bytesRead, msg] = fread( Shim.ComPort, ...
    %                           Shim.Params.nBytesToRead,  'uint8' ) ;
    % Shim.Data.input = uint8( Shim.Data.input ) ;
    %
    % if bytesRead > 0
    %     isMsgRead = true;
    % else
    %     disp('No bytes read')
    % end



end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [Cmd] = getcommands( )
%GETCOMMANDS
%
% [] = GETCOMMANDS( )
% -------------------------------------------------------------------------
% System commands (as strings) 

% Cmd.sync                        = '02' ;
% Cmd.getSystemHeartbeat          = '09' ;
Cmd.setAndLoadShim   = 's' ;

Cmd.getChannelOutput = 'p' ;
Cmd.resetAllShims    = 'q' ;

end
% =========================================================================
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
%     error('Wtf kind of computer is u running?')
% end


% -------------------------------------------------------------------------
% Serial Port 

if ismac
    % portName = '/dev/tty.usbserial' ; % USB to serial adapter
    portName = '/dev/tty.usbmodem1421' % Ryan's 2012 MacBook Pro -- right USB port
elseif isunix
    portName = '/dev/ttyS0' ;
elseif ispc
    portName = 'COM4' ;
else
    error('Wtf kind of computer is this?')
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
function adcVoltage = adctovoltage( Shim, adcCount )
%ADCTOVOLTAGE
%
% adcVoltage = ADCTOVOLTAGE( Shim, adcCount )
%
% Returns ADC voltage reading (in mV) based on adcCount (the raw ADC 
% digitization)

adcReadingVoltage = adcCount * Shim.Adc.mvPerAdcCount ;

end
% =========================================================================
% =========================================================================
function dacCount = voltstodac( Shims, voltage )
%VOLTSTODAC
%
% Wraps a voltage (in V) to a count within the DAC's range. 
%

MAX_VOLTAGE = 2.5 ; % mV 

MAX_DIGI = (2^(Shims.Specs.Dac.resolution-1)) - 1 ; % Digital to Analog Converter max value

dacCount = int16( current*( MAX_DIGI/Shims.Specs.Dac.maxCurrent  ) ) ;

end
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================

end
