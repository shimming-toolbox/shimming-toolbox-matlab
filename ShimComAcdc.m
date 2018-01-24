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
isAckReceived=true;
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
    
command =strcat('i',num2str(voltages(1)),',',num2str(voltages(2)),',',num2str(voltages(3)),',',num2str(voltages(4)),...
',',num2str(voltages(5)),',',num2str(voltages(6)),',',num2str(voltages(7)),',',num2str(voltages(8)),',');


fprintf(Shim.ComPort,'%s',command,'sync');    

for i=1:9
a = fscanf(Shim.ComPort,'%s');
disp(a);
end
end
% =========================================================================
function [] = resetallshims( Shim ) 
%RESETALLSHIMS  
%
% Reset all shims to 0 A.

    
% Send Command set all channels to 0V--------------------------------------
fprintf(Shim.ComPort,'%s','w','sync');

% Send Command to querry the channels input--------------------------------
fprintf(Shim.ComPort,'%s','q','sync');

%Read the Feedback from the arduino----------------------------------------
for i=1:9
a = fscanf(Shim.ComPort,'%s');
disp(a);
end

end

% =========================================================================
function [] = Open_ComPort( Shim ) 
fopen(Shim.ComPort);

for i=1:7
a = fscanf(Shim.ComPort,'%s');
disp(a);
end
end
% =========================================================================

function [] = close_ComPort( Shim ) 
fclose(Shim.ComPort);

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

%ChannelOutput.voltage = 

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
    portName = '/dev/cu.usbmodem1421' ;
elseif isunix
    portName = '/dev/ttyS0' ;
elseif ispc
    portName = 'COM4' ;
else
    error('Wtf kind of computer is this?')
end

delete (instrfindall) ;

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
