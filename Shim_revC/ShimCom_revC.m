classdef ShimCom_revC < ShimCom 
%SHIMCOM_revC - Shim Communication for the 8-channel AC/DC neck coil 
%
% .......
%   
% Usage
%
%   Shims = ShimCom_revC(  )
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
%    ShimCom_revC is a ShimCom subclass.
%
% =========================================================================
% Updated::20181028::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimCom_revC( Specs )
%SHIMCOM - Shim Communication

if nargin < 2 || isempty( Specs ) 
    Shim.Specs = ShimSpecs_revC( );
end

Shim.ComPort = ShimCom_revC.initializecomport( Shim.Specs ) ;
Shim.Cmd     = ShimCom_revC.getcommands( ) ;

Shim.Data.output = '' ;
Shim.Data.input  = '' ;

Shim.Params.nBytesToRead     = [] ; % depends on cmd sent to system
Shim.Params.nSendAttemptsMax = 5; % # communication attempts before error

Shim.opencomport() ;

isAckReceived = Shim.getsystemheartbeat() ;

if isAckReceived
    isCalibrationSuccessful = Shim.calibratedac() ;
    if ~isCalibrationSuccessful
        warning('DAC calibration unsuccessful. Check power supply is ON and redo Shim.calibratedac().' ) ; 
    else
        display('DAC calibration successful. Ready to issue shim currents');
    end
else
    warning('System unresponsive');
end

end
% =========================================================================
function [isAckReceived] = getsystemheartbeat( Shim ) ;
%GETSYSTEMHEARTBEAT

Shim.Params.nBytesToRead = 1 ;   
Shim.Data.output = Shim.Cmd.getSystemHeartbeat ;

Shim.sendcmd( ) ;

isAckReceived = logical( str2num( fgetl( Shim.ComPort ) ) ) 

% if isSendOk
%     Shim.isackreceived() ;
% end

end
% =========================================================================
function [] = setandloadshim( Shim, iCh, current )  ;
%SETANDLOADSHIM
%
% Set shim current (in units of Amps) for single channel 
% 
% [] = SETANDLOADSHIM( Shims, channelIndex, current ) 

Shim.Params.nBytesToRead = 1 ;   

assert( ( round(iCh) == iCh ) & ( iCh > 0 ) & ( iCh <= Shim.Specs.Amp.nChannels ), ...
    ['Channel index must be an integer between 1 and ' num2str( Shim.Specs.Amp.nChannels)] ) ;

Shim.Data.output        = Shim.Cmd.setAndLoadShimByChannel ;
Shim.Data.output(end+1) = num2str( iCh - 1 ) ; % Arduino uses 0-based indexing
Shim.Data.output = [Shim.Data.output Shim.currenttostring( current )] ;
Shim.Data.output
Shim.sendcmd() ;

isSet = logical( str2num( fgetl( Shim.ComPort ) ) ) 

end
% =========================================================================
function [] = setandloadallshims( Shim, currents )
%SETANDLOADALLSHIM
% 
% [] = SETANDLOADALLSHIMS( Shim, currents ) 
%
% Update all channels with currents (8-element vector w/units in A)

Shim.Params.nBytesToRead = 1 ;   

assert( length(currents) == Shim.Specs.Amp.nChannels, ...
    'Vector of currents must possess an entry for each shim channel.'  ) ;

Shim.Data.output = Shim.Cmd.setAndLoadAllShims ;

for iCh = 1 : Shim.Specs.Amp.nChannels
    Shim.Data.output = [ Shim.Data.output Shim.currenttostring( currents(iCh) ) ] ;
end

Shim.sendcmd() ;

isSet = logical( str2num( fgetl( Shim.ComPort ) ) ) ;

end
% =========================================================================
function [] = setandrampallshims( Shim, currents )
%SETANDRAMPALLSHIMS
% 
% [] = SETANDRAMPALLSHIMS( Shim, currents ) 
%
% Update all channels with currents (8-element vector w/units in A) by ramping 
% current up over 1.0 s

Shim.Params.nBytesToRead = 1 ;   

assert( length(currents) == Shim.Specs.Amp.nChannels, ...
    'Vector of currents must possess an entry for each shim channel.'  ) ;

Shim.Data.output = Shim.Cmd.setAndRampAllShims ;

for iCh = 1 : Shim.Specs.Amp.nChannels
    Shim.Data.output = [ Shim.Data.output Shim.currenttostring( currents(iCh) ) ] ;
end

Shim.sendcmd() ;

isSet = logical( str2num( fgetl( Shim.ComPort ) ) ) 

end
% =========================================================================
function [] = resetallshims( Shim )
%RESETALLSHIMS 
%
%   Shim currents reset to 0 A
%
Shim.Params.nBytesToRead = 1 ;   

Shim.Data.output = Shim.Cmd.resetAllShims ;

Shim.sendcmd( ) ;

isReset = logical( str2num( fgetl( Shim.ComPort ) ) ) 

assert( isReset, 'Msg send fail' )

% [Shims, isSendOk] = Shims.sendcmd( Shim.Cmd.resetAllShims ) ;
%
% if isSendOk
%     Shims.isackreceived() ;
% end

end
% =========================================================================
function [] = opencomport( Shim ) 
%OPENCOMPORT
% 
% Open serial communication port & reset Arduino Board 

if isempty( instrfind() )
    Shim.ComPort = ShimCom_revC.initializecomport( Shim.Specs ) ;
end

if strcmp( Shim.ComPort.Status, 'open' )
    return;
else
    fopen( Shim.ComPort );
    display( 'Opening ComPort & waiting for Arduino response byte...' )

    % Arduino should output a single char upon initialization
    inByte = logical( str2num( fgetl( Shim.ComPort ) ) ) ;

    if ~inByte
        warning('System unresponsive. Check connections.') ;
        fclose( Shim.ComPort ) ;
    end
end

end
% =========================================================================
function [] = closecomport( Shim ) 
%CLOSECOMPORT
% 
% Close serial communication port 

fclose(Shim.ComPort);

end
% =========================================================================
function [ChannelOutput] = getchanneloutput( Shim , ~, iChannel )
%GETCHANNELOUTPUT
%
% [ChannelOutput] = getchanneloutput( Shim , [], iChannel )
 
ChannelOutputs = Shim.getallchanneloutputs() ;
ChannelOutput.current = ChannelOutputs.current( iChannel ) ;
ChannelOutput.voltage = ChannelOutputs.voltage( iChannel ) ;

end
% =========================================================================
function [ChannelOutputs] = getallchanneloutputs( Shim )
%GETALLCHANNELSOUTPUTS      
%
% ChannelOutputs = GETALLCHANNELOUTPUTS( Shim ) 
% 
% ChannelOutputs has fields
%
%   .current [units: A]
%   .voltage [units: mV]

ChannelOutputs.current = zeros( 1, Shim.Specs.Amp.nChannels ) ;
%ChannelOutputs.voltage = zeros( 1, Shim.Specs.Amp.nChannels ) ;

Shim.Params.nBytesToRead = Shim.Specs.Amp.nChannels ;   

Shim.Data.output = Shim.Cmd.getAllChannelCurrents ;

Shim.sendcmd( ) ;

% for iCh = 1 : Shim.Specs.Amp.nChannels
%     ChannelOutputs.current(iCh) = Shim.Data.input(iCh) ; 
% end
%
% Shim.Data.input = '' ;
%
% nBytesToRead *ignoring \t and \n !
for iCh = 1 : Shim.Specs.Amp.nChannels
    ChannelOutputs.current(iCh) = str2double( fgetl( Shim.ComPort ) ) ;
end

assert( ~isempty( ChannelOutputs.current(1) ), 'Msg send fail' )

%Shim.Data.output = Shim.Cmd.getAllChannelVoltages ;
%Shim.sendcmd( ) ;

% for iCh = 1 : Shim.Specs.Amp.nChannels
%     ChannelOutputs.voltage(iCh) = str2num( fgetl( Shim.ComPort ) ) ; 
% end

end
% =========================================================================
function [isSendOk] = sendcmd( Shim, command )
%SENDCMD 
% 
%   Transmits command from client to shim microcontroller 

isSendOk  = false ;

if strcmp( Shim.ComPort.Status, 'closed' ) ;
    Shim.opencomport();
end 

fwrite( Shim.ComPort, Shim.Data.output, 'char' ) ;
pause( Shim.Specs.Com.txRxDelay ) ;

% assert( Shim.ComPort.BytesAvailable > 0, 'System unresponsive' ) ; %  Arduino must have issued at least 1 char in response
% Shim.Data.input = '' ;
%
% % nBytesToRead *ignoring \t and \n !
% for iByte = 1 : Shim.Params.nBytesToRead
%     Shim.Data.input(iByte) = str2num( fgetl( Shim.ComPort ) ) ;
% end
%
% if ~isempty( Shim.Data.input )
%     isSendOk = true;
% end
%
% assert( isSendOk, 'Msg send fail' ) ;    

end
% =========================================================================
function [isCalibrationSuccesful, isChannelCalibrationSuccesful] = calibratedac( Shim )
%CALIBRATEDAC 

% Arduino prints TRUE/FALSE for each channel according to the success of its calibration 
Shim.Params.nBytesToRead = Shim.Specs.Amp.nChannels ; 

isCalibrationSuccesful = zeros( Shim.Params.nBytesToRead, 1 ) ;

Shim.Data.output = Shim.Cmd.calibrateDac ;
Shim.sendcmd( ) ;

% read Arduino response
for iByte = 1 : Shim.Params.nBytesToRead 
    isCalibrationSuccesful(iByte) = str2num( fgetl( Shim.ComPort ) ) ; 
end

% result per channel
isChannelCalibrationSuccesful = isCalibrationSuccesful 

% overall result (true if each channel succesful)
isCalibrationSuccesful = all( isCalibrationSuccesful ) 

end
% =========================================================================
function [current] = currenttostring( Shim, current )
%CURRENTTOSTRING 
%
% Scale current (float in amperes) to uint16, convert to string, and if the resulting length is < 5,
% pad with leading '0':

assert( numel(current) == 1, '1 channel at a time...' ) ;

% shift current to be >=0... then multiply by scaling factor
current = num2str( uint16( ( current/2 + Shim.Specs.Amp.maxCurrentPerChannel(1) ) ...
    *65535/(2*Shim.Specs.Amp.maxCurrentPerChannel(1) ) ) ) ;

for i0 = numel(current)+1 : 5
    current = ['0' current] ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access = private)
% =========================================================================
end
% =========================================================================
% =========================================================================
    
methods(Static)
% =========================================================================
function [Cmd] = getcommands( )
% getcommands :
%
% - Get shim system commands 
%
%--------------------------------------------------------------------------


% System commands (as strings)---------------------------------------------

Cmd.setAndLoadAllShims    = 'a';

Cmd.calibrateDac          = 'c';

Cmd.setAndLoadShimByChannel = 'e';
% Cmd.setAndRampShimByChannel   = 'f';

Cmd.getSystemHeartbeat    = 'h'; 


Cmd.getAllChannelCurrents  = 'q'; 

Cmd.resetAllShims         = 'r' ;

Cmd.getAllChannelVoltages  = 'v'; 

Cmd.resetArduino           = 'z';

end
% =========================================================================
function [ComPort] = initializecomport( Specs )
% initializecomport : 
%
% -  Initialize (RS-232) Communication Port
% 
%--------------------------------------------------------------------------

% Serial Port 
%
warning( 'Serial port device name may change depending on the host computer.' )

if ismac
   % portName = '/dev/cu.usbserial' ;          % USB to serial adapter
    portName = '/dev/tty.usbserial' ;          % USB to serial adapter
    % portName = '/dev/tty.usbmodem14101' ; % direct USB to Arduino
    % portName = '/dev/tty.usbmodem48679001' ; % direct USB to Arduino
elseif isunix
    portName = '/dev/ttyUSB0' ;   
elseif ispc
    portName = 'COM4' ;
else
    error('What kind of computer is this!?')
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


end
% =========================================================================

end
