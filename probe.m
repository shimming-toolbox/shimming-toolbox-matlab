function [] = probe( Params )
%PROBE
%
%   Description
%   
%   Reads pressure data from microcontroller 
%
%
%   Syntax
%
%   PROBE( Parameters )
%
%    .......................
%   
%   The following Parameter.fields are supported
%   
%   .isPlottingInRealtime
%       default = false
%
%       Plots pressure vs. time during acquisition. 
%
%   .isSavingPressureLog
%       default = true
%
%   .pressureLogFilename
%       default = [./pressureLog.bin]
%
%   .pressureTimesFilename
%       default = [PATH_TO_PRESSURE_LOG /pressureTimes.bin]
%
%   .isForcingOverwrite
%       default = false
%
%       Will overwrite the log files if they already exist. 
%
%   .runTime 
%       default = 60
%   
%       Total sampling time in seconds.
%
%   .skipSampleFactor  
%       default = 1 
% 
%       The sampling rate of the microcontroller itself is fixed at 100 Hz.
%       100 Hz is therefore the maximum sampling rate.
%       Should computational efficiency become an issue (MATLAB not keeping up
%       with the microcontroller), lesser rates are possible by skipping samples.
%       E.g. Should one wish to skip every third sample, Options.skipSampleFactor
%       (which must always be an integer > 0) would be set to 3, corresponding 
%       to an effective sample rate of 100/3 Hz.
%     
%    To be supported:
%
%    .portName
%
%
%
%
% Based on sonde.m written by Imanne Al Maachi.
% 
% Updated 08/2015 by Ryan
% topfer@ualberta.ca
% =========================================================================

DEFAULT_ARDUINOPERIOD    = 10 ; % [units: ms] ;
DEFAULT_SKIPSAMPLEFACTOR = 1;

DEFAULT_ISPLOTTINGINREALTIME  = false ;
DEFAULT_ISPLOTTINGPRESSURELOG = true ;

DEFAULT_ISSAVINGPRESSURELOG   = true ;
DEFAULT_PRESSURELOGFILENAME   = './pressureLog.bin' ;
DEFAULT_PRESSURETIMESFILENAME = '/pressureTimes.bin' ;
DEFAULT_ISFORCINGOVERWRITE    = false ;

DEFAULT_RUNTIME       = 60 ; % [units: s] ;

% =========================================================================

% =========================================================================
if nargin < 1 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'isPlottingInRealtime' ) || isempty(Params.isPlottingInRealtime)
    Params.isPlottingInRealtime = DEFAULT_ISPLOTTINGINREALTIME ;
end

if  ~myisfield( Params, 'isPlottingPressureLog' ) || isempty(Params.isPlottingPressureLog)
    Params.isPlottingPressureLog = DEFAULT_ISPLOTTINGPRESSURELOG ;
end

if  ~myisfield( Params, 'isSavingPressureLog' ) || isempty(Params.isSavingPressureLog)
    Params.isSavingPressureLog = DEFAULT_ISSAVINGPRESSURELOG ;
end

if  ~myisfield( Params, 'pressureLogFilename' ) || isempty(Params.pressureLogFilename)
    Params.pressureLogFilename = DEFAULT_PRESSURELOGFILENAME ;
end

if  ~myisfield( Params, 'pressureTimesFilename' ) || isempty(Params.pressureTimesFilename)
    [pathToSamplesFilename,~,~] = fileparts(Params.pressureLogFilename) ;
    Params.pressureTimesFilename  = [pathToSamplesFilename DEFAULT_PRESSURETIMESFILENAME] ; 
end

if  ~myisfield( Params, 'isForcingOverwrite' ) || isempty(Params.isForcingOverwrite)
    Params.isForcingOverwrite = DEFAULT_ISFORCINGOVERWRITE ;
end

if Params.isSavingPressureLog && ~Params.isForcingOverwrite 
        if exist( Params.pressureLogFilename )
            disp( [Params.pressureLogFilename ' already exists.'] ) 
            Params.pressureLogFilename = [Params.pressureLogFilename '_' datestr(now)];
            disp( ['Writing to ' Params.pressureLogFilename ' instead.'] ) 
        end
        if exist( Params.pressureTimesFilename )
            disp( [Params.pressureTimesFilename ' already exists.'] ) 
            Params.pressureTimesFilename = [Params.pressureTimesFilename '_' datestr(now)];
            disp( ['Writing to ' Params.pressureTimesFilename ' instead.'] ) 
            disp( 'To force overwrite, set input Params.isForcingOverwrite to 1.' )
        end
end

if  ~myisfield( Params, 'runTime' ) || isempty(Params.runTime)
    Params.runTime = DEFAULT_RUNTIME ;
end

if  ~myisfield( Params, 'skipSampleFactor' ) || isempty(Params.skipSampleFactor)
    Params.skipSampleFactor = DEFAULT_SKIPSAMPLEFACTOR ;
end


Params.arduinoPeriod = Params.skipSampleFactor * DEFAULT_ARDUINOPERIOD ;


% =========================================================================

% =========================================================================
%% Create the serial object
if ismac
    portName = '/dev/tty.usbmodem1411' ; % my left USB port
    % portName = '/dev/tty.usbmodem1431' ; % my right USB port
    % flag = exist( portName, 'file' )
    % if flag == 0
    %     error(['Device file ( ' portName ...
    %         ' ) not found. Is the microcontroller connected?']) 
    % end
% elseif isunix
%     portName = '/dev/ttyS100';
%     flag = exist( portName, 'file' )
%     if flag == 0
%         errorMsg = ['Device file ( ' portName ' ) not found.'  ...
%             'Is the microcontroller connected?' ...
%             'Ensure symbolic link exists between actual Arduino device' ...
%             'file and the phantom-device file discoverable by MATLAB:' ...
%             'In a terminal, type: ' ... 
%             'ln /dev/ttyACM0 /dev/ttyS100' ] ; 
%     end
end

dbstop in probe at 163
ComPort = serial(portName);
fopen(ComPort);
% =========================================================================

% =========================================================================


%% send refresh time 
isSamplingFrequencyReceived = false;

iAttempt    = 1 ;
maxCommunicationAttempts = 3; 
disp('Assigning sample frequency via com port')

while(~isSamplingFrequencyReceived && iAttempt <= maxCommunicationAttempts )
    
    disp(['Attempt #' num2str(iAttempt)]);    
    
    fprintf(ComPort, '%d', Params.arduinoPeriod);
    
    pause(0.1)
    
    start_word = fscanf(ComPort,'%f')
    
    if(~isempty(start_word) && isnumeric(start_word) && start_word == Params.arduinoPeriod/10)  
        isSamplingFrequencyReceived = true;
        disp('Communication successful')
    end
    
    iAttempt = iAttempt + 1;

end




%% Collect data


if ~Params.isPlottingInRealtime 

    nSamples = floor(Params.runTime /(Params.arduinoPeriod*0.001)) ;
    
    
    fprintf('Debut enregistrement,\n');
    
    pressureTimes = zeros(nSamples,1) ;
    pressureLog   = zeros(nSamples,1) ;
    

    tic ;
    for iSample = 1 : nSamples
        pressureLog(iSample) = fscanf(ComPort,'%f');
        pressureTimes(iSample) = toc;
    end
    
    if Params.isSavingPressureLog

       pressureLogFid   = fopen( Params.pressureLogFilename, 'w+' )
       fwrite( pressureLogFid, pressureLog, 'double' ) ;
       fclose( pressureLogFid );
   
       pressureTimesFid = fopen( Params.pressureTimesFilename, 'w+' )
       fwrite( pressureTimesFid, pressureTimes, 'double' ) ;
       fclose( pressureTimesFid );
    
    end

    %% Clean up 
    fclose(ComPort);    delete(ComPort);       clear ComPort;
    
    if Params.isPlottingPressureLog 
        figure ;
        plot( pressureTimes,pressureLog);
        xlabel('Time (s)');
        ylabel('Pressure (0.01 mbar)');
    end

else
    %% Set up the figure window
    time = now;
    pressureLog = 0;
    timeStep = 0;
    
    figureHandle = figure('NumberTitle','off',...
        'Name','Sonde de respiration',...
        'Color',[0 0 0],'Visible','off');
    
    % Set axes
    axesHandle = axes('Parent',figureHandle,...
        'YGrid','on',...
        'YColor',[0.9725 0.9725 0.9725],...
        'XGrid','on',...
        'XColor',[0.9725 0.9725 0.9725],...
        'Color',[0 0 0]);
    
    hold on;
    
    plotHandle = plot(axesHandle,time,pressureLog,'Marker','.','LineWidth',1,'Color',[1 0 0]);
    
    xlim(axesHandle,[min(time) max(time+0.001)]);
    
    % Create xlabel
    xlabel('Time','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    
    % Create ylabel
    ylabel('Pression (0.01mbar)','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    
    % Create title
    title('Cycles de respiration ','FontSize',15,'Color',[1 1 0]);
    drawnow; 
    stopTime = datestr(now+Params.runTime/(24*3600),'mm/DD HH:MM:SS');
    count = 1;
    while ~isequal(datestr(now,'mm/DD HH:MM:SS'),stopTime)
        time(count) = datenum(clock);
        feed = fscanf(ComPort,'%f');
        if(   (~isempty(feed)) && (isnumeric(feed))   )
            pressureLog(count)= feed;
            set(plotHandle,'YData',pressureLog,'XData',time);
            set(figureHandle,'Visible','on');
            datetick('x','mm/DD HH:MM:SS');
            
            count = count +1;
        end
        
       drawnow;
    end
end

% =========================================================================

% =========================================================================

end
