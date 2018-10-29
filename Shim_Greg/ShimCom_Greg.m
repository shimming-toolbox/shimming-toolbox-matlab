classdef ShimCom_Greg < ShimCom 
%SHIMCOM_GREG - Shim Communication for 8ch AC/DC neck coil 
%
% .......
%   
% Usage
%
%   Shims = ShimCom_Greg(  )
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
%    ShimCom_Greg is a ShimCom subclass.
%
% =========================================================================
% Updated::20181028::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimCom_Greg( Specs )
%SHIMCOM - Shim Communication

if nargin < 2 || isempty( Specs ) 
    Shim.Specs = ShimSpecs_Greg( );
end

Shim.ComPort = ShimCom_Greg.initializecomport( Shim.Specs ) ;
Shim.Cmd     = ShimCom_Greg.getcommands( ) ;

Shim.Data.output = uint8(0) ;
Shim.Data.input  = uint8(0) ;

Shim.Params.nBytesToRead     = [] ; % depends on cmd sent to system
Shim.Params.nSendAttemptsMax = 5; % # communication attempts before error

end
% =========================================================================
function [isAckReceived] = getsystemheartbeat( Shim ) ;
%GETSYSTEMHEARTBEAT

warning('Unimplemented funct. ShimCom_Greg.GETSYSTEMHEARTBEAT()')
isAckReceived=true;

end
% =========================================================================
function [] = setandloadshim( Shim, channel, current )  ;
%SETANDLOADSHIM
%
% Set shim current (in units of Amps) for single channel 
% 
% [] = SETANDLOADSHIM( Shims, channelIndex, current ) 

% % TODO
% % temp. fix: scaling to mA
% current = current*1000;
%
% calibrationVal = ( current - Shim.Specs.Com.feedbackcalibrationcoeffy(channel) )/ Shim.Specs.Com.feedbackcalibrationcoeffx(channel) ;
%
% % Variable used to convert currents into DAC value-------------------------
%
%   preampresistance = 0.22;
%   DACmaxvalue      = 26214;
%   
% %Conversion----------------------------------------------------------------
%
%   DACcurrent = num2str((( Shim.Specs.Dac.referenceVoltage - calibrationVal * 0.001 * preampresistance) * DACmaxvalue));
%   Channel=num2str(channel);
%   
%   
%   command=strcat(Shim.Cmd.updateOneChannel,Channel,'_',DACcurrent);
%
% fprintf(Shim.ComPort,'%s',command,'sync');  

% dbstop in setandloadshim at 99
[currentLB, currentHB] = Shims.splitint( Shims.ampstodac( currents ) ) ;    
    
Shims.Data.output(iByteOut) = currentHB ;
Shims.Data.output(iByteOut) = currentLB ; 


end
% =========================================================================
function [] = setandloadallshims( Shim, currents )
%SETANDLOADALLSHIM
% 
% [] = SETANDLOADALLSHIMS( Shim, currents ) 
%
% Update all channels with currents (8-element vector w/units in A)

% TODO
% temp. fix: scaling to mA


% command = strcat('o',num2str(currentsDac(1)),num2str(currentsDac(2)),num2str(currentsDac(3)),num2str(currentsDac(4)),...
%           num2str(currentsDac(5)),num2str(currentsDac(6)),num2str(currentsDac(7)),num2str(currentsDac(8)));
%
% Shim.sendcmd( command ) ;
%
% fscanf(Shim.ComPort,'%s');

end
% =========================================================================
function [] = resetallshims( Shim )
%RESETALLSHIMS 
%
%   Shim currents reset to 0 A
%
Shims.Params.nBytesToRead = 1 ;   


[Shims, isSendOk] = Shims.sendcmd( Shim.Cmd.resetAllShims ) ;

if isSendOk
    Shims.isackreceived() ;
end

end
% =========================================================================
function [] = opencomport( Shim ) 
%OPENCOMPORT
% 
% Open serial communication port & reset Arduino Board 

instrument = instrfind;

if isempty(instrument)
    Shim.ComPort = ShimCom_Greg.initializecomport( Shim.Specs ) ;
end

fopen(Shim.ComPort);

% fprintf(Shim.ComPort,'%s',Shim.Cmd.resetArduino,'sync');

nlinesFeedback = 8 ; % Lines of feedback after opening serial communication

% Read the Feedback from the arduino---------------------------------------
for i=1:nlinesFeedback
    a = fscanf(Shim.ComPort,'%s')
end

end
% =========================================================================
function [] = closecomport( Shim ) 
%CLOSECOMPORT
% 
% Close serial communication port 

fclose(Shim.ComPort);
delete(Shim.ComPort);
clear Shim.ComPort;

end
% =========================================================================
function [ChannelOutput] = getchanneloutput( Shim , ~, Channel )
%GETCHANNELOUTPUT
% 
%   Querry and display a current feedback for one channel
 
%Command to querry the channel feedback------------------------------------
command = strcat(Shim.Cmd.getChannelFeedback,Channel);

fprintf(Shim.ComPort,'%s',command,'sync');

% Read the Feedback from the arduino---------------------------------------
ChannelOutput=fscanf(Shim.ComPort,'%s');

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

ChannelOutputs.current  = zeros( 1, Shim.Specs.Amp.nChannels ) ;
ChannelOutputs.voltages = zeros( 1, Shim.Specs.Amp.nChannels ) ;
 
Shim.sendcmd( Shim.Cmd.getAllChannelCurrents ) ;

for iCh = 1 : Shim.Specs.Amp.nChannels
    ChannelOutputs.current(iCh) = str2double( fscanf( Shim.ComPort,'%s' ) ); 
end

Shim.sendcmd( Shim.Cmd.getAllChannelVoltages ) ;

for iCh = 1 : Shim.Specs.Amp.nChannels
    ChannelOutputs.voltage(iCh) = str2double( fscanf( Shim.ComPort,'%s' ) ); 
end

end
% =========================================================================
function [Shims, isSendOk] = sendcmd( Shim, command )
%SENDCMD 
% 
%   Transmits command from client to shim microcontroller 

isSendOk  = false ;

if strcmp( Shim.ComPort.Status, 'closed' ) ;
    fopen( Shim.ComPort ) ;
end 
    
fprintf( Shim.ComPort, '%s', command, 'sync' ) ;
    
end
% =========================================================================
function [dacCorrections]= calibratedac( Shim)
%CALIBRATEDAC 
%
dacCorrections = [] ;

% if strcmp( Shim.ComPort.Status, 'closed' ) ;
%     fopen( Shim.ComPort ) ;
% end 
%     
% fprintf( Shim.ComPort,'%s', Shim.Cmd.calibrateDacCurrent,'sync' ) ;
%
% calibrationvalues = [] ;
% ncalibrationpoint = 5 ; %Number of calibration points to calculate coefficient
% calibrationvalues = zeros(ncalibrationpoint,8);
% display('Calibration in process, please wait')
% pause(41);
%
% % Read Feedback from the arduino-------------------------------------------
% for j=1:(Shim.Specs.Amp.nChannels)
%     for i=1:(ncalibrationpoint)  
%         a = fscanf(Shim.ComPort,'%s');
%         display(a);
%         calibrationvalues(i,j) = str2double(a);
%     end
% end
%
% Shim.resetallshims() ;
%
% display(calibrationvalues);

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access = private)
% =========================================================================
function [current] = currenttostring( Shim, current )
%CURRENTTOSTRING 
%
% Scale current (float in amperes) to uint16, convert to string, and if the resulting length is < 5,
% pad with leading '0':
%
% e.g.
%
%

% shift current to be >=0... then multiply by scaling factor
current = string( uint16( ( current + Shim.Specs.Amp.maxCurrentPerChannel(1) ) ...
    *65535/( 2*Shim.Specs.Amp.maxCurrentPerChannel(1) ) ) ) ;

assert( numel(current)<=5 )

for i0 = numel(current) : 5
    current = ['0' current] ;
end

end
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

Cmd.setAndLoadByChannel   = 'i';

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
    % portName = '/dev/tty.usbserial' ;          % USB to serial adapter
    portName = '/dev/tty.usbmodem14101' ;
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
