classdef ShimCom_Rriyan < ShimCom 
%SHIMCOM_RRIYAN - Shim Communication for RRI system
%
% .......
% 
% Description
%
%   SHIMCOMRRI is responsible for all direct communication with the shim system
%   (MXD & DSU). Declaration of a ShimComRri object immediately opens a serial
%   (Com) port.
% 
% .......
%   
% Usage
%
%   Shims = ShimComRri(  )
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
%   MXD and DSU commands are listed in the RRI Hex Protocol Specification guide:
%   9700052-0000 HexProtocolSpecification_REV-G
%
%   Large portions of the following code are effectively a translation, into
%   MATLAB, of the VB source code for "Demo DSU Host Software ver. 2.00"
%   courtesy of Resonance Research, Inc.  
%
%   For a primer on RS-232 communication in Matlab see
% http://www.mathworks.com/help/matlab/matlab_external/overview-of-the-serial-port.html
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
%    ShimComRri is a ShimCom subclass.
%
% =========================================================================
% Updated::20161121::ryan.topfer@polymtl.ca
% =========================================================================


% *** TODO 
% 
% ..... 
%   Error codes (system responses) as enumerations? 
% ..... 
%   Params.isDebug 
%   possible parameter for outputting system responses to screen
%   
%   e.g. When true: Notify user if 'isSendOk' + what type of system responses
%   are returned. 
% ..... 
%   Shims.Params.INVALID -
%   Main.vb and HexProtocol (manual pdf) define this differently. Therefore,
%   bytes 3-5 of INVALID might be wrong.
% ..... 
%   SETCHANNELBUFFER()
%   does this work? clean up &/or fix.
% ..... 
%   Way of checking what the ComPort names are?
%   +Optional input params:
%       portName 
% ..... 
%   .Status as property?
% ..... 
%   checkcontrolresponse(  )
%
%   only control responses have a single data bit returned.
%   replace ISACKRECEIVED() with this.
%   Whenever the length of the received data differs from expectation
%   (i.e. when it is 5) call CHECKCONTROLRESPONSE()
% .....   
%   Lookup table:
%   definitions of MXD/DSU errors from MXD manual such that a returned error
%   msg to the serial port can indicate sth useful to the user
% ..... 
% SPLITINT( ) 
%  
%   Documentation for how two's complement works.
% ..... 
% When converting to DAC counts (e.g ampstodac())
%   consider output 'clipFlag' for instance where input current exceeds 
%   max allowable current?
%
% Likewise for any set() function in ShimComRri for shim currents.
%     
% =========================================================================

% properties   
%     Cmd;
%     ComPort;
%     Data;
%     Params;
%     Specs;
% end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimCom_Rriyan(  )
%SHIMCOM - Shim Communication

if nargin < 2
    Shims.Specs = ShimSpecs_Rriyan( );
end

Shims.ComPort = ShimCom_Rriyan.initializecomport( Shims.Specs ) ;
Shims.Cmd     = ShimCom_Rriyan.getcommands( ) ;

Shims.Data.output = uint8(0) ;
Shims.Data.input  = uint8(0) ;

Shims.Params.nBytesToRead     = [] ; % depends on cmd sent to system
Shims.Params.nSendAttemptsMax = 5; % # communication attempts before error
% -------
% Format for ACK/NACK/INVALID:
% { Sync; TotalMsgLength; CommandItself; C0; C1 } ; 
% where C0 & C1 are the CRC bytes 

Shims.Params.ACK             = {'02'; '05'; '00'; '00'; 'CB'} ;
Shims.Params.Unimplemented   = '01' ;
Shims.Params.ParamOutOfRange = '02' ;
Shims.Params.ErrorPresent    = '03' ;
Shims.Params.NACK            = {'02'; '05'; 'FF'; '78'; 'C4'} ;
Shims.Params.INVALID         = {'02'; '05'; '08'; '08'; '4F'} ;



end
% =========================================================================
function [offset, slope] = getbankadccalibrationdata( Shims, iBank ) 
%GETBANKADCCALIBRATIONDATA      (*INVALID? MXD cmd 0x02)
%
% [offset, slope] = GETBANKADCCALIBRATIONDATA( Shims, bankIndex ) 
%
%   bankIndex = 0, 1, 2, or 3

Shims.Params.nBytesToRead = 8 ;  

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 6 ; 
Shims.Data.output(3) = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(4) = uint8( iBank ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end

Shims.Data.input
offset = double( ShimCom_Rriyan.mergeints( Shims.Data.input(3), ...
                               Shims.Data.input(4), true ) ) ;

slope = double( ShimCom_Rriyan.mergeints( Shims.Data.input(5), ...
                               Shims.Data.input(6), false ) ) ;
end
% =========================================================================
function [isAckReceived] = getsystemheartbeat( Shims ) 
%GETSYSTEMHEARTBEAT     (MXD cmd 0x09)
%
% Queries MXD and returns 'ACK' if responsive

isAckReceived = false;

Shims.Params.nBytesToRead = 5 ; % (ACK only)

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemHeartbeat ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
    disp(systemResponse)
end

end
% =========================================================================
function [] = clearsystemerrors( Shims ) 
%CLEARSYSTEMERRORS      (MXD cmd 0x0A)
%
% [] = CLEARSYSTEMERRORS( Shims ) 

Shims.Params.nBytesToRead = 5 ; % (ACK only)

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.clearSystemErrors ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end

end
% =========================================================================
function [] = resetallshims( Shim ) 
%RESETALLSHIMS      (Custom cmd)  
%
% Reset all shims to 0 A.
%
% [] = RESETALLSHIMS( Shim )

Shim.setallshims( zeros(Shim.Specs.Amp.nChannels, 1) ) ;
Shim.setloadallshims ;

end
% =========================================================================
function ThresholdData = getchannelthresholddata( Shims, iBank, iChannel ) 
%GETCHANNELTHRESHOLDDATA   (NACK? MXD cmd 0x16)
% 
% ThresholdData = getchannelthresholddata( Shims, bankIndex, channelIndex ) 
%
% ThresholdData contains fields
%   .maxControlVoltage [in volts]
%   .maxOutputCurrent [in amps]
%   .maxOutputPower [in Watts]

% return byte: (1) sync ; (2) length; (3-10) data; (11-12) crc;  
Shims.Params.nBytesToRead = 10 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;
Shims.Data.output(2) = 7 ;
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getChannelThresholdData ) ;
Shims.Data.output(4) = uint8( iBank ) ;
Shims.Data.output(5) = uint8( iChannel ) ;

Shims.Data.output    = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd( ) ;

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input );
end

% in volts
ThresholdData.maxControlVoltage = ...
     (1/1000)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(3), ...
                                      Shims.Data.input(4), true ) ) ;
% in amperes
ThresholdData.maxOutputCurrent = ...
     (1/1000)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(5), ...
                                     Shims.Data.input(6), true ) ) ;
% in Watts
ThresholdData.maxOutputPower   = ...
     double( ShimCom_Rriyan.mergeints( Shims.Data.input(7), ...
                                 Shims.Data.input(8), false ) ) ;

end
% =========================================================================
function [] = setpoweron( Shims ) 
%SETPOWERON     (MXD cmd 0x20)
% 
% [] = SETPOWERON( Shims )
Shims.Params.nBytesToRead = 5 ; % (ACK only)

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setPowerOn ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end

end
% =========================================================================
function [] = setpoweroff( Shims ) 
%SETPOWEROFF    (MXD cmd 0x21)
%
% [] = SETPOWEROFF( Shims )

Shims.Params.nBytesToRead = 5 ; % (ACK only)

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setPowerOff ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end

end
% =========================================================================
function [] = setallshims( Shims, currents ) 
%SETALLSHIMS    (MXD cmd: 0x22)
%
% SETALLSHIMS( Shims, currents )
% 
% Sets all shims based on [nChannel x 1] current vector (in amps)
% 
% NB: for our system, currents vector has 32 entries (=total # amplifier
% channels, including inactive channels).  Currents for inactive channels
% should be zero.

assert( length(currents) == Shims.Specs.Amp.nChannels )

if any( currents > Shims.Specs.Amp.maxCurrentPerChannel )
    fprintf( [ 'WARNING: \n Requested current exceeds MAX channel current of ' ...
            num2str(Shims.Specs.Amp.maxCurrentPerChannel) ...
            '. Resulting current will be clipped.\n'] ) ;
end

Shims.Params.nBytesToRead = 5 ;

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;
Shims.Data.output(2) = 69 ; % sync,lng, cmd, 32*2 ch, crc*2; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setAllShims ) ;

iByteOut = 4 ;

for iChannel = 1 : Shims.Specs.Amp.nChannels
    [currentLB, currentHB] = Shims.splitint( ...
        Shims.ampstodac( currents(iChannel) ) ) ;    
    
    Shims.Data.output(iByteOut) = currentHB ;
    iByteOut = iByteOut + 1 ;
    Shims.Data.output(iByteOut) = currentLB ; 
    iByteOut = iByteOut + 1 ;
end

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd( ) ;

if isSendOk
    [isAckReceived, systemResponse] = isackreceived( Shims ) ;
    if ~isAckReceived
        disp(systemResponse)
    end
end

end
% =========================================================================
function [] = setloadallshims( Shims ) 
%SETLOADALLSHIMS     (MXD cmd 0x23)
%
%   Loads buffered shim settings.
% 
%   [] = SETLOADALLSHIMS( Shims ) 

Shims.Params.nBytesToRead = 5 ; % (ACK only)

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setLoadAllShims ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
    if ~isAckReceived
        disp(systemResponse)
    end
end

end
% =========================================================================
function SystemStatus = getsystemlongstatus( Shims ) 
%GETSYSTEMLONGSTATUS    (MXD cmd 0x25)
%
% SystemStatus = GETSYSTEMLONGSTATUS( Shims )
%
% SystemStatus has fields
%   .totalOutputPower [in W]
%   .onOff [1 - ON, 0 - OFF] 
%   .globlErrorRegister  
%       [0-7 bitwise: Error on a bank, Bank is missing, Bank communication
%       error, TBD, TBD, TBD, Interlock Error, Fatal Error Requiring Shutdown
%       Present)]
    
% return byte: (1) sync ; (2) length; (3-6) data; (7-8) crc;  
Shims.Params.nBytesToRead = 8 ; % 

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemLongStatus ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end
% in Watts
SystemStatus.totalOutputPower = ...
    double( ShimCom_Rriyan.mergeints( Shims.Data.input(3), ...
                               Shims.Data.input(4), false ) ) ;

SystemStatus.onOff = logical( Shims.Data.input(5) ) ;

SystemStatus.globalErrorRegister = Shims.Data.input(6) ;
    
end
% =========================================================================
function BankStatus = getbanklongstatus( Shims, iBank ) 
%GETBANKLONGSTATUS      (MXD cmd 0x26)
%
% BankStatus = GETBANKLONGSTATUS( Shims, bankIndex )
%   
%   bankIndex = 0, 1, 2, or 3
%
% BankStatus has fields
%   .highVoltagePositiveRail [in V]
%   .highVoltageNegativeRail [in V]
%   .lowVoltagePositiveRail  [in V]
%   .lowVoltageNegativeRail  [in V]
%   .heatsinkTemperature     [in degrees C]
%   .channelsWithFaults     [0-7 bitwise: ch0, ch1, ch2,... ch7] 
%   .faults           [See HEX protocol document]

% return byte: (1) sync ; (2) length; (3-9) data; (10-11) crc;  
Shims.Params.nBytesToRead = 11 ; % 

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 6 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getBankLongStatus ) ;
Shims.Data.output(4) = uint8( iBank ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    disp('OK')
    % in volts 
    BankStatus.highVoltagePositiveRail = ...
        ShimCom_Rriyan.convertfromtwoscomplement( Shims.Data.input(3) ) ;
    BankStatus.highVoltageNegativeRail = ...
        ShimCom_Rriyan.convertfromtwoscomplement( Shims.Data.input(4) ) ;
    BankStatus.lowVoltagePositiveRail  = ...
        ShimCom_Rriyan.convertfromtwoscomplement( Shims.Data.input(5) ) ;
    BankStatus.lowVoltageNegativeRail  = ...
        ShimCom_Rriyan.convertfromtwoscomplement( Shims.Data.input(6) ) ;
    
    % in degrees Celsius
    BankStatus.heatsinkTemperature = Shims.Data.input(7) ;
    
    BankStatus.channelsWithFaults = Shims.Data.input(8) ;
    BankStatus.faults = Shims.Data.input(9) ;
else
    error('?')
end

end
% =========================================================================
function SystemInfo = getsysteminformation( Shims ) 
%GETSYSTEMINFORMATION (MXD cmd 0x29)
%
% SystemInfo = GETSYSTEMINFORMATION( Shims )
% 
% SystemInfo has fields
%
%   .SoftwareVersion
%       .majorRevision
%       .minorRevision
%   .SystemModel
%       .name
%       .numberOfChannels
%   .serialNumber
    
% return byte: (1) sync ; (2) length; (3-8) data; (9-10) crc;  
Shims.Params.nBytesToRead = 10 ; % 

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemInformation ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end

SystemInfo.SoftwareVersion.majorRevision = Shims.Data.input(3) ; 
SystemInfo.SoftwareVersion.minorRevision = Shims.Data.input(4) ; 

SystemInfo.SystemModel.name = Shims.Data.input(5) ; 
SystemInfo.SystemModel.numberOfChannels = Shims.Data.input(6) ; 

SystemInfo.serialNumber = ...
    double( ShimCom_Rriyan.mergeints( Shims.Data.input(7), ...
                               Shims.Data.input(8), false ) ) ;
    
end
% =========================================================================
function [] = setsystemcurrenttime( Shims, SystemTime ) 
%SETSYSTEMCURRENTTIME (MXD cmd 0x2B)
%
% [] = SETSYSTEMCURRENTTIME( Shims, SystemTime )
%
% SystemTime has fields
%   
%   .year
%   .month
%   .day
%   .hour
%   .minute
%   .second
    
Shims.Params.nBytesToRead = 5 ; % ACK only

Shims.Data.output  = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;

if nargin < 2
    error('Not enough arguments')
end


Shims.Data.output(2) = 12 ;
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setSystemCurrentTime ) ;

Shims.Data.output(4) = uint8( SystemTime.month ) ;
Shims.Data.output(5) = uint8( SystemTime.day ) ;

[yearLB, yearHB] = Shims.splitint( uint16( SystemTim.year ) ) ;    
Shims.Data.output(6) = yearHB ;
Shims.Data.output(7) = yearLB ; 

Shims.Data.output(8)  = uint8( SystemTime.hours ) ;
Shims.Data.output(9)  = uint8( SystemTime.minutes ) ;
Shims.Data.output(10) = uint8( SystemTime.seconds ) ;


Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk
    Shims.isackreceived() ;
end
    
end
% =========================================================================
function [] = setandloadallshims( Shim, currents )
%SETANDLOADALLSHIMS     (custom cmd)
% 
% Sets shim buffers (MXD cmd 0x22) and loads the settings (MXD cmd 0x23).
%
% [] = SETANDLOADALLSHIMS( Shim, currents )
%
% numel(currents) == Shims.Specs.nChannels || Shims.Specs.nActiveChannels
% 
% i.e. currents vector is either length 24 or 32
 

if numel(currents) == Shim.Specs.Amp.nChannels
    Shim.setallshims( currents ) ;
    
elseif numel(currents) == Shim.Specs.Amp.nActiveChannels
    Shim.setallshims( Shim.mapcurrentstomxd( currents ) ) ;
end

Shim.setloadallshims ;

end
% =========================================================================
function ChannelOutputs = getallchanneloutputs( Shim )
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

channelsToBankKey = Shim.getchanneltobankkey ;

ChannelOutputs.current = zeros( 1, Shim.Specs.Amp.nActiveChannels ) ;
ChannelOutputs.voltage = zeros( 1, Shim.Specs.Amp.nActiveChannels ) ;
ChannelOutputs.power   = zeros( 1, Shim.Specs.Amp.nActiveChannels ) ;
ChannelOutputs.dissipatedPower = zeros( 1, Shim.Specs.Amp.nActiveChannels ) ; 

for iChannel = 1 : Shim.Specs.Amp.nActiveChannels 
    
    ChannelOutput = Shim.getchanneloutput( channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) ) ;
    ChannelOutputs.current(iChannel)         = ChannelOutput.current ;
    ChannelOutputs.voltage(iChannel)         = ChannelOutput.voltage ;
    ChannelOutputs.power(iChannel)           = ChannelOutput.power ;
    ChannelOutputs.dissipatedPower(iChannel) = ChannelOutput.dissipatedPower ;

end

end
% =========================================================================
function SystemTime = getsystemcurrenttime( Shims ) 
%GETSYSTEMCURRENTTIME   (MXD cmd 0x2C)
%
% SystemTime = GETSYSTEMCURRENTTIME( Shims )
%
% SystemTime has fields
%   
%   .year
%   .month
%   .day
%   .hour
%   .minute
%   .second
    
% return byte: (1) sync ; (2) length; (3-10) data; (10-11) crc;  
Shims.Params.nBytesToRead = 11 ; % 

% output message to read in shim values
Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
Shims.Data.output(2) = 5 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemCurrentTime ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;
if isSendOk
    Shims.isackreceived() ;
end

Shims.Data.input(:) 
SystemTime.month = double( Shims.Data.input(3) ) ; 
SystemTime.day   = double( Shims.Data.input(4) ) ; 

SystemTime.year = double( ShimCom_Rriyan.mergeints( Shims.Data.input(5), ...
                    Shims.Data.input(6), false ) ) ;

SystemTime.hours   = double( Shims.Data.input(7) ) ; 
SystemTime.minutes = double( Shims.Data.input(8) ) ; 
SystemTime.seconds = double( Shims.Data.input(9) ) ; 
    
end
% =========================================================================
function [] = setandloadshim( Shims, varargin ) 
%SETANDLOADSHIM     (MXD cmd 0x44 || 0x54)
%
% Set shim current (in amps) for single channel 
% 
% [] = SETANDLOADSHIM( Shims, channelIndexGlobal, current ) 
% [] = SETANDLOADSHIM( Shims, bankIndex, channelIndexByBank, current ) 


Shims.Params.nBytesToRead = 5 ; 

Shims.Data.output  = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;

if nargin < 3
    error('Not enough arguments')

elseif nargin == 3

    iChannel = varargin{1} ;
    current  = varargin{2} ;

    [currentLB, currentHB] = Shims.splitint( Shims.ampstodac( current ) ) ;    

    Shims.Data.output(2) = 8 ;
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setAndLoadShimByChannel ) ;
    Shims.Data.output(4) = uint8( iChannel ) ;
    Shims.Data.output(5) = currentHB ;
    Shims.Data.output(6) = currentLB ; 

elseif nargin == 4
    iBank    = varargin{1} ;
    iChannel = varargin{2} ;
    current  = varargin{3} ;

    [currentLB, currentHB] = Shims.splitint( Shims.ampstodac( current ) ) ;    
    Shims.Data.output(2) = 9 ;
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setAndLoadShimByBankChannel ) ;
    Shims.Data.output(4) = uint8( iBank ) ;
    Shims.Data.output(5) = uint8( iChannel ) ;
    Shims.Data.output(6) = currentHB ; 
    Shims.Data.output(7) = currentLB ;
    
end
    
Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd( ) ;

if isSendOk
    [isAckReceived, systemResponse] = isackreceived( Shims ) ;
    if isAckReceived
        disp('OK')
    else
        disp(systemResponse)
    end
end

end
% =========================================================================
function [] = setchannelbuffer( Shims, iBank, iChannel, current ) 
%SETCHANNELBUFFER *INVALID?
if current > Shims.Specs.Amp.maxCurrentPerChannel
    fprintf( [ 'WARNING: \n Requested current of ' num2str(current) ...
            ' exceeds MAX channel current of ' ...
            num2str(Shims.Specs.Amp.maxCurrentPerChannel) 
            '\n Resulting current will be clipped.\n'] ) ;
end

[currentLB, currentHB] = Shims.splitint( Shims.ampstodac( current ) ) 

Shims.Params.nBytesToRead = 5 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;

% % apparently this command is INVALID
% Shims.Data.output(2) = 8 ;
% Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setShimBufferByChannel ) ;
% Shims.Data.output(4) = uint8( iChannel ) ;
% Shims.Data.output(5) = currentHB ;
% Shims.Data.output(6) = currentLB ; 

% apparently this command is also INVALID
Shims.Data.output(2) = 9 ;
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setShimBufferByBankChannel ) ;
Shims.Data.output(4) = uint8( iBank ) ;
Shims.Data.output(5) = uint8( iChannel ) ;
Shims.Data.output(6) = currentHB ;
Shims.Data.output(7) = currentLB ; 


Shims.Data.output    = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd( ) ;

if isSendOk
    [isAckReceived, systemResponse] = isackreceived( Shims ) ;
    if isAckReceived
        disp('OK')
    else
        disp(systemResponse)
    end
end

end
% =========================================================================
function ChannelStatus = getchannellongstatus( Shims, iBank, iChannel ) 
%GETCHANNELLONGSTATUS   (MXD cmd 0x45)
% 
% ChannelStatus = getchannellongstatus( Shims, bankIndex, channelIndex ) 
%
% ChannelStatus contains fields
%   .current [in amps]
%   .voltage [in volts]
%   .power [in Watts]
%   .dissipatedPower [in Watts]
%   .dacInputVoltage [in volts]
%   .auxiliaryInputVoltage [in volts]
%   .onOff [0=Off, 1=On]
%   .auxiliaryEnable [0=Off, 1=On]
%   .faults [see HEX protocol]

% return byte: (1) sync ; (2) length; (3-10) data; (11-12) crc;  
Shims.Params.nBytesToRead = 19 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;
Shims.Data.output(2) = 7 ;
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getChannelLongStatus ) ;
Shims.Data.output(4) = uint8( iBank ) ;
Shims.Data.output(5) = uint8( iChannel ) ;

Shims.Data.output    = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd( ) ;

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input );
end

% in amperes
ChannelStatus.current = ...
     (1/1000)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(3), ...
                                      Shims.Data.input(4), true ) ) ;
% in volts
ChannelStatus.voltage = ...
     (1/100)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(5), ...
                                     Shims.Data.input(6), true ) ) ;
% in Watts
ChannelStatus.power   = ...
     double( ShimCom_Rriyan.mergeints( Shims.Data.input(7), ...
                                 Shims.Data.input(8), false ) ) ;

ChannelStatus.dissipatedPower = ...
     double( ShimCom_Rriyan.mergeints( Shims.Data.input(9), ...
                                 Shims.Data.input(10), false ) ) ;

% in volts
ChannelStatus.dacInputVoltage = ...
     (1/100)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(11), ...
                                     Shims.Data.input(12), true ) ) ;

ChannelStatus.auxiliaryInputVoltage = ...
     (1/100)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(13), ...
                                     Shims.Data.input(14), true ) ) ;

ChannelStatus.onOff = logical( Shims.Data.input(15) ) ;

ChannelStatus.auxiliaryEnable = logical( Shims.Data.input(16) ) ;

ChannelStatus.faults = Shims.Data.input(17)  ;

end
% =========================================================================
function ChannelOutput = getchanneloutput( Shims, iBank, iChannel ) 
%GETCHANNELOUTPUT   (MXD cmd 0x47)
% 
% ChannelOutput = getchanneloutput( Shims, bankIndex, channelIndex ) 
%
% ChannelOutput contains fields
%   .current [in amps]
%   .voltage [in volts]
%   .power [in Watts]
%   .disspitatedPower [in Watts]

% return byte: (1) sync ; (2) length; (3-10) data; (11-12) crc;  
Shims.Params.nBytesToRead = 12 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;
Shims.Data.output(2) = 7 ;
Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getChannelOutput ) ;
Shims.Data.output(4) = uint8( iBank ) ;
Shims.Data.output(5) = uint8( iChannel ) ;

Shims.Data.output    = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd( ) ;

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input );
end

% in amperes
ChannelOutput.current = ...
     (1/1000)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(3), ...
                                      Shims.Data.input(4), true ) ) ;
% in volts
ChannelOutput.voltage = ...
     (1/100)*double( ShimCom_Rriyan.mergeints( Shims.Data.input(5), ...
                                     Shims.Data.input(6), true ) ) ;
% in Watts
ChannelOutput.power   = ...
     double( ShimCom_Rriyan.mergeints( Shims.Data.input(7), ...
                                 Shims.Data.input(8), false ) ) ;

ChannelOutput.dissipatedPower = ...
     double( ShimCom_Rriyan.mergeints( Shims.Data.input(9), ...
                                 Shims.Data.input(10), false ) ) ;

end
% =========================================================================
function percentSliceIntensity = getslice( Shims, iChannel, iSlice )
%GETSLICE    (DSU cmd 0xE0)
%
% percentSliceIntensity = GETSLICE( Shims, channelIndex, sliceIndex ) 
%
%   channelIndex : number between 1 to 32 (see ShimCom_Rriyan.getchanneltobankkey)
%   sliceIndex   : number between 1 and 1000 
%
%   percentSliceIntensity : number between -100 and 100  

Shims.Params.nBytesToRead = 6 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Dsu.sync) ) ; 
Shims.Data.output(2) = 8 ; 
Shims.Data.output(3) = uint8( hex2dec( Shims.Cmd.Dsu.getSlice ) ) ;
Shims.Data.output(4) = uint8( iChannel ) ;

[iSliceLowByte, iSliceHighByte] = Shims.splitint( uint16( iSlice ) ) ;

Shims.Data.output(5) = iSliceHighByte ;
Shims.Data.output(6) = iSliceLowByte ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
    if ~isAckReceived
        disp(systemResponse)
    end
end

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input ); ??
end

percentSliceIntensity = dactopercentsliceintensity( ...
    ShimCom_Rriyan.mergeints( Shims.Data.input(3), Shims.Data.input(4), true ) ) ;

function ps = dactopercentsliceintensity( dc )
% Convert dac counts (dc, between 0 and 65535)) to percent slice intensity (ps,
% between -100 to 100)  
    ps = double( dc ) - 32768 ; % offset
    ps = (100/32768)*ps ; % scale
end

end
% =========================================================================
function [] = setslice( Shims, iChannel, iSlice, percentSliceIntensity )
%SETSLICE    (DSU cmd 0xD0)
%
% [] = SETSLICE( Shims, channelIndex, sliceIndex, sliceIntensity ) 
%
%   channelIndex : number between 1 to 32 (see ShimCom_Rriyan.getchanneltobankkey)
%   sliceIndex   : number between 1 and 1000 
%   percentSliceIntensity : number between -100 and 100  

assert( abs(percentSliceIntensity) <= 100, ...
        'Input slice intensity (%) must be between -100 and 100' )

Shims.Params.nBytesToRead = 5 ; % ACK only

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Dsu.sync) ) ; 
Shims.Data.output(2) = 10 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Dsu.setSlice ) ;
Shims.Data.output(4) = uint8( iChannel ) ;

[iSliceLowByte, iSliceHighByte] = Shims.splitint( uint16( iSlice ) ) ;

Shims.Data.output(5) = iSliceHighByte ;
Shims.Data.output(6) = iSliceLowByte ;

[percentSliceIntensityLowByte, percentSliceIntensityHighByte] = ...
    Shims.splitint( percentsliceintensitytodac( percentSliceIntensity)  ) ;

Shims.Data.output(7) = percentSliceIntensityHighByte ;
Shims.Data.output(8) = percentSliceIntensityLowByte ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
    if ~isAckReceived
        disp(systemResponse)
    end
end

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input ); ??
end

function dc = percentsliceintensitytodac( ps )
% convert percent slice intensity (ps, between -100 to 100) to DAC counts 
% (dc, between 0 and 65535)
    assert( abs(percentSliceIntensity) <= 100, ...
        'Input slice intensity (%) must be between -100 and 100' )
    dc = ps + 100 ; % floor set to 0x0000
    dc = dc*32767.5/100 ;
    dc = uint16( dc ) ;
end

end
% =========================================================================
function [] = setmatrixinterconnectiontoself( Shims, iChannel, iSlot )
%SETMATRIXINTERCONNECTIONTOSELF    (DSU cmd 0x20)
%
% [] = SETMATRIXINTERCONNECTIONTOSELF( Shims, channelIndex, slotIndex ) 
%
%   channelIndex : number between 1 to 32 (see ShimCom_Rriyan.getchanneltobankkey)
%   slotIndex : [time constant module 1-4]: 
%       0x00 (no connection), 
%       0x01 (first module), 
%       0x02 (second), 
%       0x04 (third), 
%       0x08 (fourth), 
%       0x09 (first and fourth) etc. 
%
%       (See corresponding entry in Hex Communication manual)

Shims.Params.nBytesToRead = 5 ; % ACK only. 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Dsu.sync) ) ; 
Shims.Data.output(2) = 7 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Dsu.setMatrixInterconnectionToSelf ) ;
Shims.Data.output(4) = uint8( iChannel ) ;
Shims.Data.output(5) = uint8( hex2dec( iSlot ) ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
    if ~isAckReceived
        disp(systemResponse)
    end
end

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input ); ??
end

end
% =========================================================================
function interconnections = getchannelmatrixinterconnection( Shims, iChannel )
%GETCHANNELMATRIXINTERCONNECTION (DSU cmd 0xFB)
%
% interconnections = GETMATRIXINTERCONNECTION( Shims, channelIndex ) 
%
%   channelIndex : number between 1 to 32 (see ShimCom_Rriyan.getchanneltobankkey)


Shims.Params.nBytesToRead = 68 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Dsu.sync) ) ; 
Shims.Data.output(2) = 6 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Dsu.getMatrixInterconnections ) ;
Shims.Data.output(4) = uint8( iChannel ) ;

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
    if ~isAckReceived
        disp(systemResponse)
    end
end

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input ); ??
end

interconnections = zeros(32,2) ;

for iChannel = 3 : 2 : 65
    interconnections( iChannel, 1 ) = Shims.Data.input( iChannel ) ;
    interconnections( iChannel, 2 ) = Shims.Data.input( iChannel + 1) ;
end

end
% =========================================================================
function [] = testwave( Shims, iChannel, isOn )
%TESTWAVE    (DSU cmd 0x21)
%
% [] = TESTWAVE( Shims, channelIndex, isOn ) 
%
%   channelIndex : number between 1 to 32 (see ShimCom_Rriyan.getchanneltobankkey)
%   isOn : true (1) OR false (0) 

Shims.Params.nBytesToRead = 5 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Dsu.sync) ) ; 
Shims.Data.output(2) = 7 ; 
Shims.Data.output(3) = hex2dec( Shims.Cmd.Dsu.testWave ) ;
Shims.Data.output(4) = uint8( iChannel ) ;

if isOn == 1
    Shims.Data.output(5) = uint8( 1 ) ;
else
    Shims.Data.output(5) = uint8( 2 ) ;
end

Shims.Data.output = ShimCom_Rriyan.appendcrc( Shims.Data.output, ...
                        ShimCom_Rriyan.calculatecrc( Shims.Data.output ) ) ; 

[Shims, isSendOk] = Shims.sendcmd() ;

if isSendOk 
    [isAckReceived, systemResponse] = Shims.isackreceived( ) 
    if ~isAckReceived
        disp(systemResponse)
    end
end

if( Shims.Data.input(2) ~= Shims.Params.nBytesToRead ) 
    disp('error...')
    % Shims.checkcontrolresponse( Shims.Data.input ); ??
end

end
%=========================================================================
function [] = delete( Shim  )
%DELETE  (custom helper function)
% 
% DELETE( Shim )
% 
% Destructor. Calls Shim.deletecomport( ) 

Shim.deletecomport();
clear Shim ;

end
% =========================================================================
function [Shims] = deletecomport( Shims )
%DELETECOMPORT  (custom helper function)
% 
% Shims = DELETECOMPORT( Shims )
% 
% Correct way to delete and clear the serial port object

% check existence?
if myisfield( Shims, 'ComPort' ) && ~isempty( Shims.ComPort ) 
    fclose( Shims.ComPort ) ;
    delete( Shims.ComPort ) ;
    % clear Shims.ComPort  ;
end

end
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access = private)
% =========================================================================
function [isAckReceived, systemResponse] = isackreceived( Shims )
%ISACKRECEIVED
%
% Checks data returned from master board for 
%   'ACK','NACK','INVALID', or 'Wrong' 

isAckReceived = false ;

% '3' because 1st is SYNC, 2nd byte is LNG,...
systemResponse = Shims.Data.input(3) ; 

if systemResponse == hex2dec( Shims.Params.ACK(3) )
    systemResponse = 'ACK';
    isAckReceived  = true;

elseif systemResponse == hex2dec( Shims.Params.NACK(3) )
    systemResponse = 'NACK'; 

elseif systemResponse == hex2dec( Shims.Params.INVALID(3) )
    systemResponse = 'INVALID'; 

else
    systemResponse = 'Wrong' ;
end

end
% =========================================================================
function [Shims, isSendOk]= sendcmd( Shims )
%SENDCMD transmits a command from client to MXD or DSU systems
%   order of output: 
%
%   sync -> lng  -> cmd  -> outputData -> C0 -> C1
%
% .......
%
%   Description
%
%   sync 
%
%       sync command for MXD (=='0x02') or DSU (=='0xC0') 
%
%   lng
%
%       length in bytes of entire message to be transmitted 
%       (including sync, lng, cmd, ,outputData, crc)
%
%   cmd 
%
%       system command in hex format
%
%   outputData
%
%       formated according to command (see RRI manual)
%
%
% .......

isSendOk     = false ;
iSendAttempt = 0 ;

while ~isSendOk && (iSendAttempt < Shims.Params.nSendAttemptsMax) 

    [Shims, isMsgRead] = Shims.writetomachine() ;
    if isMsgRead % check message was read properly:
        crc = ShimCom_Rriyan.calculatecrc( Shims.Data.input( 1 :end-2) ) ;
    
        if ( Shims.Data.input( end - 1) == bitand(crc, uint16(255) ) ) ...
            && ( Shims.Data.input( end ) == bitshift( crc, -8 ) ) 
            isSendOk = true ;
            return 
        end 
    end
    iSendAttempt = iSendAttempt + 1 ;
end

error('System unresponsive.')

end
% =========================================================================
function [Shims, isMsgRead] = writetomachine( Shims )
%WRITETOMACHINE
%
% Write to machine (and read back the system response)

isMsgRead = false ;

if strcmp( Shims.ComPort.Status, 'closed' ) ;
    fopen( Shims.ComPort ) ;
end 



fwrite( Shims.ComPort, Shims.Data.output, 'uint8' ) ;
pause( Shims.Specs.Com.txRxDelay ) ;

[Shims.Data.input, bytesRead, msg] = fread( Shims.ComPort, ...
                          Shims.Params.nBytesToRead,  'uint8' ) ;
Shims.Data.input = uint8( Shims.Data.input ) ;



if bytesRead > 0
    isMsgRead = true;
else
    disp('No bytes read')
end

end
% =========================================================================
function dacCount = ampstodac( Shims, current )
%AMPSTODAC
%
% Wraps a current value in Amperes to a count within the DAC's range. 
%
% dacCount = AMPSTODAC( Shims, current ) 
%
% Note: Beware of clipping:
%
%   Shims.Specs.Dac.maxCurrent = 5 A % [By default]
%
%   if current >= Shims.Specs.Dac.maxCurrent, dacCount = 32767 ;
%   elseif current <= Shims.Specs.Dac.maxCurrent, dacCount = -32768 ;

MAX_DIGI = (2^(Shims.Specs.Dac.resolution-1)) - 1 ; % Digital to Analog Converter max value

dacCount = int16( current*( MAX_DIGI/Shims.Specs.Dac.maxCurrent  ) ) ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [X0,X1,X2,X3] = getchanneltobankmatrices( )
%GETCHANNELTOBANKMATRICES   (custom helper function)
%
% [X0,X1,X2,X3] = GETCHANNELTOBANKMATRICES()
%
% Returns four [nChannelsPerBank x 24] matrices that when applied to 
% the [24 x 1] vector of shim currents will each return a truncated version
% of the vector corresponding strictly to the channels belonging to amplifier
% banks 0 through 3.

channelsToBankKey = ShimCom_Rriyan.getchanneltobankkey() ;

nActiveChannels   = size( channelsToBankKey, 1 ) ;
banks             = channelsToBankKey(:,2) ;
nBanks            = max( banks ) + 1 ;

X                 = cell( nBanks, 1 ) ;

for iBank = 0 : nBanks - 1
    
    nActiveChannelsPerBank = sum( banks == iBank ) ;
    
    X{iBank + 1} = full( sparse( ...
                    [1:nActiveChannelsPerBank], find( banks == iBank ), ...
                    ones(nActiveChannelsPerBank,1), ...
                    nActiveChannelsPerBank, nActiveChannels ) )  ;

end 

X0 = X{1} ;
X1 = X{2} ;
X2 = X{3} ;
X3 = X{4} ;

end
% =========================================================================
function [currentsOut] = mapcurrentstomxd( currents )
%MAPCURRENTSTOMXD   (custom helper function)
% 
% There are only 24 active channels in the 32 channel amplifier. MAPCURRENTSTOMXD
% takes a 24-element vector and rearranges the elements into a 32-element vector
% where indices corresponding to inactive channels are 0 and the original 1-24
% channel values are arranged according to their global MXD channel index.
%
% currentsOut = MAPCURRENTSTOMXD( currentsIn )
%
% e.g. Insert active channel indices into corresponding cells of full currents
%   vector (i.e. 32 channels, including inactive channels)
%
% >> currentsOut = ShimCom_Rriyan.mapcurrentstomxd( [1 : 24] ) 
%
% currentsOut =
%
%      2
%      5
%      6
%      9
%     10
%      1
%      3
%      0
%      4
%      7
%      8
%     11
%     12
%      0
%      0
%      0
%     14
%     17
%     18
%     21
%     22
%     13
%     15
%      0
%     16
%     19
%     20
%     23
%     24
%      0
%      0
%      0

channelToBankKey = ShimCom_Rriyan.getchanneltobankkey() ;
currentsOut = zeros( 32, 1 ) ;
currentsOut( channelToBankKey( :, 4 ) ) = currents ;

end    
% =========================================================================
function [channelToBankKey] = getchanneltobankkey( )
%GETCHANNELTOBANKKEY    (custom helper function)
%
% Maps the spine coil channels (1-24) to the corresponding MXD banks (0-3) 
% & channels (1-32)
% 
% .......   
% 
% Syntax
%
%   X = GETCHANNELTOBANKKEY();
%
% .......   
% 
% Description
% 
%   E.g. (top 6 rows of output)
%   
%   >> getchanneltobankkey()
% 
% Shim Coil || Bank || Channel || Channel || 
%                     (by bank)  (global)
%     1         0        5          6           
%     2         0        0          1           
%     3         0        6          7           
%     4         1        0          9           
%     5         0        1          2 
%     6         0        2          3 
%
%                   ...
%
% .......

nChannels = 24 ;

channelToBankKey       = zeros( nChannels, 4 ) ;
channelToBankKey(:, 1) = 1 : nChannels ;

% correspondence
channelToBankKey(1, 2:4)  = [0 5 6] ;
channelToBankKey(2, 2:4)  = [0 0 1] ;
channelToBankKey(3, 2:4)  = [0 6 7] ;
channelToBankKey(4, 2:4)  = [1 0 9] ;
channelToBankKey(5, 2:4)  = [0 1 2] ;
channelToBankKey(6, 2:4)  = [0 2 3] ;

channelToBankKey(7, 2:4)  = [1 1 10] ;
channelToBankKey(8, 2:4)  = [1 2 11] ;
channelToBankKey(9, 2:4)  = [0 3 4]  ;
channelToBankKey(10, 2:4) = [0 4 5]  ;
channelToBankKey(11, 2:4) = [1 3 12] ;
channelToBankKey(12, 2:4) = [1 4 13] ;

channelToBankKey(13, 2:4) = [2 5 22] ;
channelToBankKey(14, 2:4) = [2 0 17] ;
channelToBankKey(15, 2:4) = [2 6 23] ;
channelToBankKey(16, 2:4) = [3 0 25] ;
channelToBankKey(17, 2:4) = [2 1 18] ;
channelToBankKey(18, 2:4) = [2 2 19] ;

channelToBankKey(19, 2:4) = [3 1 26] ;
channelToBankKey(20, 2:4) = [3 2 27] ;
channelToBankKey(21, 2:4) = [2 3 20] ;
channelToBankKey(22, 2:4) = [2 4 21] ;
channelToBankKey(23, 2:4) = [3 3 28] ;
channelToBankKey(24, 2:4) = [3 4 29] ;

end
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static, Hidden=true)
% =========================================================================
function [dataMsg] = appendcrc( dataMsg, crc )
%APPENDCRC
% 
% uint8 msg = APPENDCRC( uint8 msg, uint16 crc ) 
%
% appendedMsg = APPENDCRC( msg, crc ) adds the 2 (low followed by high) CRC
% bytes to the end of 'msg'

assert( isa( crc, 'uint16' ) ) ;
assert( isa( dataMsg, 'uint8' ) ) ;

[lowByte, highByte] = ShimCom_Rriyan.splitint( crc ) ; 

dataMsg(end+1) = lowByte ;    
dataMsg(end+1) = highByte ;

end
% =========================================================================
function [crc] = calculatecrc( dataMsg )
%CALCULATECRC
%   
%   uint16 crc = CALCULATECRC( uint8 dataMsg ) ;
%
%   dataMsg
%       integer array
%
%   crc
%       integer
    
assert( strcmp( class( dataMsg ), 'uint8' ) ) ;

dataMsg = uint16( dataMsg ) ;

poly = uint16( hex2dec( '8408' ) ) ;
crc  =  uint16( 0 ) ;  

for iMsg = 1 : length( dataMsg )

    crc = bitxor( crc, dataMsg(iMsg) ) ;

    for iBitShift = 1 : 8 

        if (bitand(crc,1)~=0)
            crc = bitshift(crc,-1);
            crc = bitxor(crc,poly);
        else
            crc = bitshift(crc,-1);
        end
    end
end

end
% =========================================================================
function [Cmd] = getcommands( )
%GETCOMMANDS
%
% [] = GETCOMMANDS( )
% -------------------------------------------------------------------------
% System commands (as Hexadecimal strings) 
% I.  for the MXD 
% II. for DSU

Cmd.Mxd.sync                        = '02' ;
Cmd.Mxd.getSystemHeartbeat          = '09' ;
Cmd.Mxd.clearSystemErrors           = '0A' ;
Cmd.Mxd.getChannelThresholdData     = '16' ;
Cmd.Mxd.setPowerOn                  = '20' ;
Cmd.Mxd.setPowerOff                 = '21' ;
Cmd.Mxd.setAllShims                 = '22' ;
Cmd.Mxd.setLoadAllShims             = '23' ;
% Cmd.Mxd.getSystemStatShort        = '24' ;
Cmd.Mxd.getSystemLongStatus         = '25' ;
Cmd.Mxd.getBankLongStatus           = '26' ;
% Cmd.Mxd.getBankShortStatus        = '27' ;
% Cmd.Mxd.getBankChannelsStat       = '28' ;
Cmd.Mxd.getSystemInformation        = '29' ;

% Cmd.Mxd.getSystemHistoryFile  = '2A' ;
Cmd.Mxd.setSystemCurrentTime       = '2B' ;
Cmd.Mxd.getSystemCurrentTime       = '2C' ;
% Cmd.Mxd.setAllShimsOn       = '2D' ;
% Cmd.Mxd.setAllShimsOff       = '2E' ;


Cmd.Mxd.setAndLoadShimByBankChannel = '44' ;
Cmd.Mxd.getChannelLongStatus        = '45' ;
Cmd.Mxd.getChannelOutput            = '47' ;

% Cmd.Mxd.setChannelPowerON         = '50' ;
% Cmd.Mxd.setChannelPowerOFF        = '51' ;
% Cmd.Mxd.setChannelAuxEnable       = '52' ;
% Cmd.Mxd.setChannelAuxDisable      = '53' ;
Cmd.Mxd.setAndLoadShimByChannel     = '54' ;
% Cmd.Mxd.getChannelStatLong        = '55' ;
% Cmd.Mxd.getChannelStatShort       = '56' ;
% Cmd.Mxd.getChannelOutput          = '57' ;
% Cmd.Mxd.getChannelInput           = '58' ;
% Cmd.Mxd.getChannelControl         = '59' ;
Cmd.Mxd.setShimBufferByBankChannel  = '4A' ;
Cmd.Mxd.loadShimByBankChannel       = '4B' ;

Cmd.Mxd.setShimBufferByChannel      = '5A' ;
Cmd.Mxd.loadShimByChannel           = '5B' ;

% -------------------------------------------------------------------------
% II. Hex commands for DSU


% Cmd.Dsu.heartbeat                     = '09';
% Cmd.Dsu.intTrigger                    = '10';
%

Cmd.Dsu.setMatrixInterconnectionToSelf  = '20';
Cmd.Dsu.testWave                        = '21';
% Cmd.Dsu.getEcc                        = 'A0';
% Cmd.Dsu.setEcc                        = 'B0';
Cmd.Dsu.sync                          = 'C0';

Cmd.Dsu.setSlice                      = 'D0';
Cmd.Dsu.getSlice                      = 'E0';
%
% Cmd.Dsu.diagnosticsOn                 = 'D2';
% Cmd.Dsu.triggerOn                     = 'D3';
% Cmd.Dsu.triggerOff                    = 'D4';
% Cmd.Dsu.triggerStat                   = 'D5';
% Cmd.Dsu.stopPulse                     = 'D6';
%
% Cmd.Dsu.errorReport                   = 'F0';
% Cmd.Dsu.writeRampRate                 = 'F1';
% Cmd.Dsu.readRampRate                  = 'F2';
%
% Cmd.Dsu.setNewSliceNumber             = 'F4';
% Cmd.Dsu.resetToSlice01                = 'F5';
% Cmd.Dsu.resetToSlice02                = 'F6';
% Cmd.Dsu.setLoopLastSliceLoopReady     = 'F7';
% Cmd.Dsu.stopLoop                      = 'F8';
% Cmd.Dsu.pulseShapeSequence            = 'F9';
% Cmd.Dsu.setMatrixIntercon             = 'FA';
Cmd.Dsu.getMatrixInterconnections       = 'FB';
% Cmd.Dsu.getTcData                     = 'FD';
% Cmd.Dsu.getTcConfig                   = 'FE';
% Cmd.Dsu.triggerAck                    = 'FF';

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
    portName = '/dev/tty.usbserial' ;
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
end
% =========================================================================
% =========================================================================

end

