classdef ShimCom
%SHIMCOM - Shim Communication
%
% Shims = ShimCom(  )
%
%   Shims contains fields
%
%       .Cmd
%           
%       .ComPort    
%
%       .Data 
%
%       .Parameters
% .......
%
%   Description
%
%   Declaration of ShimCom object immediately opens a serial (aka Com) port
%   (see ShimCom.opencomport() for the expected port names).
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%     ShimCom
%     ShimOpt
%     ShimSpecs
% 
% To Write:
%     ShimUse (should be high-level commands for user + private methods)
%     ShimGui (GUI to call ShimUse functions)
%
% =========================================================================
% NB
%
%   Large portions of the following code are effectively a translation, into
%   MATLAB, of the VB source code for "Demo DSU Host Software ver. 2.00"
%   courtesy of Resonance Research, Inc.  
%
% For primer on RS-232 communication in Matlab see 
%
% http://www.mathworks.com/help/matlab/matlab_external/overview-of-the-serial-
% port.html
%
%   Careful with implicit casting in matlab -
%   Hex Communication Protocol for the shims deals exclusively with 8-bit 
%   words. By default, however, Matlab sets numerical variables to type double.
% =========================================================================
% *** TODO 
% 
%     .debug parameter or something? When true: Notify user if 'isSendOk' + what type
%     of system responses are returned. 
% ..... 
% Shims.Parameters.INVALID -
%     Main.vb and HexProtocol (manual pdf) define this differently.
%     Therefore, bytes 3-5 of INVALID might be wrong.
%
% ..... 
%     Way of checking what the ComPort names are?
%
% ..... 
%     .Status as property?
%
% ..... 
%     -Optional input params:
%           portName 
%    
% ..... 
%     -listcommands( 'usage' )
%          isImplemented('cmdName')? 
% ..... 
%     -checkcontrolresponse(  )
%
%   only control responses have a single data bit returned.
%   replace ISACKRECEIVED() with this.
%   Whenever the length of the received data differs from expectation
%   (i.e. when it is 5) call CHECKCONTROLRESPONSE()
%
% ..... 
% CALCULATECRC() 
%
%   shouldn't be static. 
%     
% ..... 
%     
%     -Lookup table:
%      definitions of MXD/DSU errors from MXD manual such that a 
%      returned error msg to the serial port can indicate sth useful to the
%      user
% ..... 
% OPENCOMPORT( ) 
%
%   should check if the port is already in use - if it is, clear the buffer.
%   Likewise CLOSECOMPORT should make sure the serial object exists?
%
% ..... 
% SPLITINT( ) 
%  
%   Documentation for how two's complement works.
%
% ..... 
% MAPCOILSTOMXD()
%
%   currently performs 2 functions: the mapping of coil index (1-24) to bank/channel
%   and also the display function, showing current for each channel. The latter 
%   should be a separate function entirely, and should not be within this class.
%
%   also,
%   skewed 4th column output corresponding to current for MXD channels > 9 ??
%
%     
% =========================================================================

properties   

    Cmd;
    ComPort;
    Data;
    Parameters;
    Specs;

end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shims = ShimCom(  )
%SHIMCOM - Shim Communication

    % if nargin < 1
    %     ShimUse('debug') ; 
    if nargin < 2
        Shims.Specs = ShimSpecs();
    end

    Shims.Cmd     = ShimCom.definecommands( ) ;
    Shims.ComPort = ShimCom.initialisecomport( Shims.Specs ) ;

    Shims.Data.output = uint8(0) ;
    Shims.Data.input  = uint8(0) ;

    % -------
    % Format for ACK/NACK/INVALID:
    % { Sync; TotalMsgLength; CommandItself; C0; C1 } ; 
    % where C0 & C1 are the CRC bytes 
    %
    Shims.Parameters.ACK     = {'02'; '05'; '00'; '00'; 'CB'} ;
    Shims.Parameters.Unimplemented = '01';
    Shims.Parameters.ParamOutOfRange = '02' ;
    Shims.Parameters.ErrorPresent = '03' ;
    Shims.Parameters.NACK    = {'02'; '05'; 'FF'; '78'; 'C4'} ;
    Shims.Parameters.INVALID = {'02'; '05'; '08'; '08'; '4F'} ;
    

    Shims.Parameters.nBytesToRead     = [] ; % depends on cmd sent to system
    Shims.Parameters.nSendAttemptsMax = 5; % # communication attempts before error

    % Shims.Parameters.isConnectedMxd   = true ;
    % Shims.Parameters.isConnectedDsu   = true ;

    
    % % -------
    % % ?? could absorb isSendOk into this
    % Shims.status = '';
end
% =========================================================================
function [] = clearsystemerrors( Shims ) 
    
    Shims.Parameters.nBytesToRead = 5 ; % (ACK only)

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.clearSystemErrors ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk
        Shims.isackreceived() ;
    end
end
% =========================================================================
function SystemStatus = getsystemlongstatus( Shims ) 
% SystemStatus = GETSYSTEMLONGSTATUS( Shims )
    
    % return byte: (1) sync ; (2) length; (3-6) data; (7-8) crc;  
    Shims.Parameters.nBytesToRead = 8 ; % 

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemLongStatus ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk
        Shims.isackreceived() ;
    end

    % in Watts
    SystemStatus.totalOutputPower = ...
        double( ShimCom.mergeints( Shims.Data.input(3), ...
                                   Shims.Data.input(4), false ) ) ;
    
    SystemStatus.onOff = Shims.Data.input(5) ;

    SystemStatus.globalErrorRegister = Shims.Data.input(6) ;
    
end
% =========================================================================
function BankStatus = getbanklongstatus( Shims, iBank ) 
% BankStatus = GETBANKLONGSTATUS( Shims, bankIndex )
%   where bankIndex = 0, 1, 2, or 3
%
% BankStatus has fields
%   .highVoltagePositiveRail [in V]
%   .highVoltageNegativeRail [in V]
%   .lowVoltagePositiveRail  [in V]
%   .lowVoltageNegativeRail  [in V]
%   .heatsinkTemperature     [in degrees C]
%   .channelsWithFaults      [    ]
%   .faultsPresent           [    ]

    % return byte: (1) sync ; (2) length; (3-9) data; (10-11) crc;  
    Shims.Parameters.nBytesToRead = 11 ; % 

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 6 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getBankLongStatus ) ;
    Shims.Data.output(4) = uint8( iBank ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk
        disp('OK')
        % in volts 
        BankStatus.highVoltagePositiveRail = ...
            ShimCom.convertfromtwoscomplement( Shims.Data.input(3) ) ;
        BankStatus.highVoltageNegativeRail = ...
            ShimCom.convertfromtwoscomplement( Shims.Data.input(4) ) ;
        BankStatus.lowVoltagePositiveRail  = ...
            ShimCom.convertfromtwoscomplement( Shims.Data.input(5) ) ;
        BankStatus.lowVoltageNegativeRail  = ...
            ShimCom.convertfromtwoscomplement( Shims.Data.input(6) ) ;
        
        % in degrees Celsius
        BankStatus.heatsinkTemperature = Shims.Data.input(7) ;
        
        BankStatus.channelsWithFaults = Shims.Data.input(8) ;
        BankStatus.faultsPresent = Shims.Data.input(9) ;
    else
        error('?')
    end
end
% =========================================================================
function SystemInfo = getsysteminformation( Shims ) 
% SystemInfo = GETSYSTEMINFORMATION( Shims )
    
    % return byte: (1) sync ; (2) length; (3-6) data; (7-8) crc;  
    Shims.Parameters.nBytesToRead = 10 ; % 

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemInformation ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk
        Shims.isackreceived() ;
    end

    SystemInfo.SoftwareVersion.majorRevision = Shims.Data.input(3) ; 
    SystemInfo.SoftwareVersion.minorRevision = Shims.Data.input(4) ; 
    
    SystemInfo.SystemModel.name = Shims.Data.input(5) ; 
    SystemInfo.SystemModel.numberOfChannels = Shims.Data.input(6) ; 
    
    SystemInfo.serialNumber = ...
        double( ShimCom.mergeints( Shims.Data.input(7), ...
                                   Shims.Data.input(8), false ) ) ;
    
    
end
% =========================================================================
function [] = resetallshims( Shims ) 
%RESETALLSHIMS
%
% RESETALLSHIMS( Shims )
% 
% Sets all shims to 0 amps.

    Shims.setallshims( zeros(Shims.Specs.Amp.nChannels, 1) ) ;
    Shims.setloadallshims ;

end
% =========================================================================
function [] = setallshims( Shims, currents ) 
%SETALLSHIMS
%
% SETALLSHIMS( Shims, currents )
% 
% Sets all shims based on [nChannel x 1] current vector (in amps)
% 
% (e.g. for our system, currents vector has 32 (and not 24) entries.
% Currents for inactive channels should be zero.

    assert( length(currents) == Shims.Specs.Amp.nChannels )

    if any( currents > Shims.Specs.Amp.maxCurrentPerChannel )
        fprintf( [ 'WARNING: \n Requested current ' ...
                ' exceeds MAX channel current of ' ...
                num2str(Shims.Specs.Amp.maxCurrentPerChannel) 
                '\n Resulting current will be clipped.\n'] ) ;
    end

    Shims.Parameters.nBytesToRead = 5 ;

    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;
    Shims.Data.output(2) = 69 ; % sync,lng, cmd, 32*2 ch, crc*2; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setAllShims ) ;
    
    iByteOut = 4 ;

    for iChannel = 1 : Shims.Specs.Amp.nChannels
        [currentLB, currentHB] = Shims.splitint( ...
            Shims.Specs.ampstoint( currents(iChannel) ) ) ;    
        
        Shims.Data.output(iByteOut) = currentHB ;
        iByteOut = iByteOut + 1 ;
        Shims.Data.output(iByteOut) = currentLB ; 
        iByteOut = iByteOut + 1 ;
    end

    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 

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
function [] = setandloadshim( Shims, varargin ) 
%SETANDLOADSHIM
%
% Set shim current for single channel 
% 
% systemResponse = setandloadshim( Shims, channelIndexGlobal, current ) 
% systemResponse = setandloadshim( Shims, bankIndex, channelIndexByBank, current ) 
%
% current is in amps
%
% systemResponse is a HEX code. For the list of responses, type 
%   help definesystemresponses 
%
% -------
Shims.Parameters.nBytesToRead = 5 ; 

Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;

if nargin < 3
    error('Not enough arguments')

elseif nargin == 3

    iChannel = varargin{1} ;
    current  = varargin{2} ;

    [currentLB, currentHB] = Shims.splitint( Shims.Specs.ampstoint( current ) ) ;    

    Shims.Data.output(2) = 8 ;
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setAndLoadShimByChannel ) ;
    Shims.Data.output(4) = uint8( iChannel ) ;
    Shims.Data.output(5) = currentHB ;
    Shims.Data.output(6) = currentLB ; 

elseif nargin == 4
    iBank    = varargin{1} ;
    iChannel = varargin{2} ;
    current  = varargin{3} ;

    [currentLB, currentHB] = Shims.splitint( Shims.Specs.ampstoint( current ) ) ;    
    Shims.Data.output(2) = 9 ;
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setAndLoadShimByBankChannel ) ;
    Shims.Data.output(4) = uint8( iBank ) ;
    Shims.Data.output(5) = uint8( iChannel ) ;
    Shims.Data.output(6) = currentHB ; 
    Shims.Data.output(7) = currentLB ;
    
end
    
    Shims.Data.output    = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 

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
%SETCHANNELBUFFER
    if current > Shims.Specs.Amp.maxCurrentPerChannel
        fprintf( [ 'WARNING: \n Requested current of ' num2str(current) ...
                ' exceeds MAX channel current of ' ...
                num2str(Shims.Specs.Amp.maxCurrentPerChannel) 
                '\n Resulting current will be clipped.\n'] ) ;
    end

    [currentLB, currentHB] = Shims.splitint( Shims.Specs.ampstoint( current ) ) 

    Shims.Parameters.nBytesToRead = 5 ; 

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
    

    Shims.Data.output    = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 

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
function [] = setpoweroff( Shims ) 
% [] = SETPOWEROFF( Shims )
    Shims.Parameters.nBytesToRead = 5 ; % (ACK only)

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setPowerOff ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk
        Shims.isackreceived() ;
    end
end
% =========================================================================
function [] = setpoweron( Shims ) 
% [] = SETPOWERON( Shims )
    Shims.Parameters.nBytesToRead = 5 ; % (ACK only)

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setPowerOn ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk
        Shims.isackreceived() ;
    end
end
% =========================================================================
function ChannelOutput = getchanneloutput( Shims, iBank, iChannel ) 
%GETCHANNELOUTPUT 
%
% Return channel values 
% 
% ChannelOutput = getchanneloutput( Shims, bankIndex, channelIndex ) 
%
% ChannelOutput contains fields
%   .current [in amps]
%   .voltage [in volts]
%   .power [in Watts]
%   .disspitatedPower [in Watts]

    % return byte: (1) sync ; (2) length; (3-10) data; (11-12) crc;  
    
    Shims.Parameters.nBytesToRead = 12 ; 
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ;
    Shims.Data.output(2) = 7 ;
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getChannelOutput ) ;
    Shims.Data.output(4) = uint8( iBank ) ;
    Shims.Data.output(5) = uint8( iChannel ) ;

    Shims.Data.output    = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 

    [Shims, isSendOk] = Shims.sendcmd( ) ;
    
    if( Shims.Data.input(2) ~= Shims.Parameters.nBytesToRead ) 
        disp('error...')
        % Shims.checkcontrolresponse( Shims.Data.input );
    end
    
    % in amperes
    ChannelOutput.current = ...
         (1/1000)*double( ShimCom.mergeints( Shims.Data.input(3), ...
                                          Shims.Data.input(4), true ) ) ;
    % in volts
    ChannelOutput.voltage = ...
         (1/100)*double( ShimCom.mergeints( Shims.Data.input(5), ...
                                         Shims.Data.input(6), true ) ) ;
    % in Watts
    ChannelOutput.power   = ...
         double( ShimCom.mergeints( Shims.Data.input(7), ...
                                     Shims.Data.input(8), false ) ) ;
    
    ChannelOutput.dissipatedPower = ...
         double( ShimCom.mergeints( Shims.Data.input(9), ...
                                     Shims.Data.input(10), false ) ) ;

end
% =========================================================================
function [] = getsystemheartbeat( Shims ) 
% %SENDACK
% % write ACK or NACK to the ComPort
% %
    
    Shims.Parameters.nBytesToRead = 5 ; % (ACK only)

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync ) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.getSystemHeartbeat ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
    if isSendOk 
        [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
        if isAckReceived
            disp('OK')
        else
            disp(systemResponse)
        end
    end

end
% =========================================================================
function [] = setloadallshims( Shims ) 
    
    Shims.Parameters.nBytesToRead = 5 ; % (ACK only)

    % output message to read in shim values
    Shims.Data.output    = uint8( hex2dec( Shims.Cmd.Mxd.sync) ) ; 
    Shims.Data.output(2) = 5 ; 
    Shims.Data.output(3) = hex2dec( Shims.Cmd.Mxd.setLoadAllShims ) ;
    
    Shims.Data.output = ShimCom.appendcrc( Shims.Data.output, ...
                            ShimCom.calculatecrc( Shims.Data.output ) ) ; 
    
    [Shims, isSendOk] = Shims.sendcmd() ;
   
    if isSendOk 
        [isAckReceived, systemResponse] = Shims.isackreceived( ) ;
        if isAckReceived
            disp('OK')
        else
            disp(systemResponse)
        end
    end

end
% % =========================================================================
% function [] = sendack( Shims )
% %SENDACK
% % write ACK or NACK to the ComPort
% %
%
%     if crc == 0
%         Shims.Data.output = uint8( hex2dec( Shims.Parameters.ACK ) ) ;
%         Shims.writetomachine( );
%         disp('ACK sent')
%     else
%         Shims.Data.output = uint8( hex2dec( Shims.Parameters.NACK ) ) ;
%         Shims.writetomachine( );
%         disp('NACK sent')
%     end   
%
% end
% =========================================================================
function [Shims] = deletecomport( Shims )
%DELETECOMPORT  
% Shims = DELETECOMPORT( Shims )
% 
% Proper way to delete the serial port object

% check existence?
if myisfield( Shims, 'ComPort' ) && ~isempty( Shims.ComPort ) 
    fclose( Shims.ComPort ) ;
    delete( Shims.ComPort ) ;
    clear Shims.ComPort  ;
else
    disp('ComPort is not open.');
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
 
    if systemResponse == hex2dec( Shims.Parameters.ACK(3) )
        systemResponse = 'ACK';
        isAckReceived  = true;
    
    elseif systemResponse == hex2dec( Shims.Parameters.NACK(3) )
        systemResponse = 'NACK'; 
    
    elseif systemResponse == hex2dec( Shims.Parameters.INVALID(3) )
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
   
    while ~isSendOk && (iSendAttempt < Shims.Parameters.nSendAttemptsMax) 
    
        [Shims, isMsgRead] = Shims.writetomachine() ;
        
        if isMsgRead % now check message was read properly:
            crc = ShimCom.calculatecrc( Shims.Data.input( 1 :end-2) ) ;
        
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
% (also reads back the system response)

    isMsgRead = false ;
    
    if strcmp( Shims.ComPort.Status, 'closed' ) ;
        fopen( Shims.ComPort ) ;
    else 
        %-----
        % write
        fwrite( Shims.ComPort, Shims.Data.output, 'uint8' ) ;
        pause( Shims.Specs.Com.txRxDelay ) ;
        
        % -----
        % read
        [Shims.Data.input, bytesRead, msg] = fread( Shims.ComPort, ...
                                  Shims.Parameters.nBytesToRead,  'uint8' ) ;
        Shims.Data.input = uint8( Shims.Data.input ) ;
        
        if bytesRead > 0
            isMsgRead = true;
        else
            disp('no bytes read')
        end

        % fclose( Shims.ComPort ) ;
    end

end
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
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

    [lowByte, highByte] = ShimCom.splitint( crc ) ; 

    dataMsg(end+1) = lowByte ;    
    dataMsg(end+1) = highByte ;

end
% =========================================================================
function [crc] = calculatecrc( dataMsg )
%CALCULATECRC
%
%   
%   uint16 crc = CALCULATECRC( uint8 dataMsg ) ;
%
%
%       dataMsg
%           integer array
%
%       crc
%           integer
% .......
%
%   Description
    
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
function [ComPort] = initialisecomport( Specs )
% Initialise (RS-232) Communication Port 
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

    % doesn't make sense to presuppose these names in this way!
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

    ComPort = serial(portName,...
        'BaudRate', Specs.Com.baudRate,...
        'DataBits',Specs.Com.dataBits,...
        'StopBits',Specs.Com.stopBits,...
        'FlowControl',Specs.Com.flowControl,...
        'Parity',Specs.Com.parity,...
        'ByteOrder',Specs.Com.byteOrder ) ;   

end
% =========================================================================
function [X0,X1,X2,X3] = getchanneltobankmatrices( )
%GETCHANNELTOBANKMATRICES
%
% [X0,X1,X2,X3] = GETCHANNELTOBANKMATRICES()
%
% Returns four [nChannelsPerBank x 24] matrices that when applied to 
% the [24 x 1] vector of shim currents will each return a truncated version
% of the vector corresponding strictly to the channels belonging to amplifier
% banks 0 through 3.


channelsToBankKey = ShimCom.getchanneltobankkey() ;

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

    channelToBankKey = ShimCom.getchanneltobankkey() ;

    nActiveChannels = size( channelToBankKey,1 ) ;
    currentsOut = zeros( 32, 1) ;
   
    for iChannel = 1 : nActiveChannels 
        currentsOut( channelToBankKey(iChannel,4) ) = currents(iChannel);
    end

end    
% =========================================================================
function [channelToBankKey] = getchanneltobankkey( )
%GETCHANNELTOBANKKEY
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

    % Indices: 1st column 'spine coil' | 2nd 'MXD' bank # | 3rd MXD ch #
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
%
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

    assert( isa(zUInt, 'uint8') )
    
    if zUInt < 127
        z = int16( zUInt ) ;
    else
        z = -int16(bitxor( zUInt, uint8(255) ) + 1) ;
    end

end
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
function [Cmd] = definecommands( )
%DEFINECOMMANDS
%
% [] = DEFINECOMMANDS( )
% -------------------------------------------------------------------------
% System commands (as Hexadecimal strings) 
% I.  for the MXD 
% II. for DSU

    Cmd.Mxd.sync                        = '02' ;
    
    Cmd.Mxd.getSystemHeartbeat          = '09' ;

    Cmd.Mxd.clearSystemErrors         = '0A' ;
    %
    Cmd.Mxd.setPowerOn                  = '20' ;
    Cmd.Mxd.setPowerOff                 = '21' ;
    Cmd.Mxd.setAllShims                 = '22' ;
    Cmd.Mxd.setLoadAllShims             = '23' ;
    % Cmd.Mxd.getSystemStatShort          = '24' ;
    Cmd.Mxd.getSystemLongStatus           = '25' ;
    Cmd.Mxd.getBankLongStatus             = '26' ;
    % Cmd.Mxd.getBankShortStatus          = '27' ;
    % Cmd.Mxd.getBankChannelsStat         = '28' ;
    Cmd.Mxd.getSystemInformation          = '29' ;
    



    Cmd.Mxd.setAndLoadShimByBankChannel = '44' ;
    Cmd.Mxd.getChannelOutput            = '47' ;
    
    % Cmd.Mxd.setChannelPowerON           = '50' ;
    % Cmd.Mxd.setChannelPowerOFF          = '51' ;
    % Cmd.Mxd.setChannelAuxEnable         = '52' ;
    % Cmd.Mxd.setChannelAuxDisable        = '53' ;
    Cmd.Mxd.setAndLoadShimByChannel       = '54' ;
    % Cmd.Mxd.getChannelStatLong          = '55' ;
    % Cmd.Mxd.getChannelStatShort         = '56' ;
    % Cmd.Mxd.getChannelOutput            = '57' ;
    % Cmd.Mxd.getChannelInput             = '58' ;
    % Cmd.Mxd.getChannelControl           = '59' ;
    Cmd.Mxd.setShimBufferByBankChannel    = '4A' ;
    Cmd.Mxd.loadShimByBankChannel         = '4B' ;
    
    Cmd.Mxd.setShimBufferByChannel        = '5A' ;
    Cmd.Mxd.loadShimByChannel             = '5B' ;
    
    % -------------------------------------------------------------------------
    % II. Hex commands for DSU


    % Cmd.Dsu.heartbeat                     = '09';
    % Cmd.Dsu.intTrigger                    = '10';
    %
    % Cmd.Dsu.tcConnect                     = '20';
    % Cmd.Dsu.testWave                      = '21';
    % Cmd.Dsu.getEcc                        = 'A0';
    % Cmd.Dsu.setEcc                        = 'B0';
    % Cmd.Dsu.sync                          = 'C0';
    % Cmd.Dsu.setSlice                      = 'D0';
    %
    % Cmd.Dsu.diagnosticsOn                 = 'D2';
    % Cmd.Dsu.triggerOn                     = 'D3';
    % Cmd.Dsu.triggerOff                    = 'D4';
    % Cmd.Dsu.triggerStat                   = 'D5';
    % Cmd.Dsu.stopPulse                     = 'D6';
    %
    % Cmd.Dsu.getSlice                      = 'E0';
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
    % Cmd.Dsu.getMatrixIntercon             = 'FB';
    % Cmd.Dsu.getTcData                     = 'FD';
    % Cmd.Dsu.getTcConfig                   = 'FE';
    % Cmd.Dsu.triggerAck                    = 'FF';

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





