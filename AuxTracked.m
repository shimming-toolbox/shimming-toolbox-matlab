classdef (Abstract) AuxTracked < Aux 

properties   
    % Data ;
    % Specs ;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Aux = AuxTracked( Specs )
%TRACKING  

if nargin < 1
    Specs = [] ;
end

Aux.Specs  = Specs ;
Aux.Specs.state = 'inert' ;
Aux.Data   = [] ; 
Aux.Data.p = [] ; 

end    
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Abstract)
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [measurementLog, sampleTimes] = loadmeasurementlog( measurementLogFilename, sampleTimesFilename )
%LOADMEASUREMENTLOG
% 
% Reads binary file of data measurements (e.g. pressure recording) to return
% vector(s) of doubles.
%
% measurementLog                = LOADMEASUREMENTLOG( measurementLogFilename ) ;
% [measurementLog, sampleTimes] = LOADMEASUREMENTLOG( measurementLogFilename, sampleTimesFilename )

if nargin < 1
    error( 'Insufficient arguments. Must provide full path to measurement log .bin file.' ) ;

else
    if nargin >= 1
        measurementLogFid = fopen( measurementLogFilename, 'r' ) ;
        measurementLog    = fread( measurementLogFid, inf, 'double' ) ;
        fclose( measurementLogFid );
    end

    if nargin == 2 
        sampleTimesFid = fopen( sampleTimesFilename, 'r' ) ;
        sampleTimes    = fread( sampleTimesFid, inf, 'double' ) ;
        fclose( sampleTimesFid );
    end
end

end
% =========================================================================
function [] = plotmeasurementlog( measurementLog, Params )
%PLOTMEASUREMENTLOG
%
% PLOTMEASUREMENTLOG( measurementLog ) ;
% PLOTMEASUREMENTLOG( measurementLog, Params )
%
% Supported fields to Params struct
%
%   .figureTitle
%       [default: 'Pressure log']
%
%   .sampleTimes
%       vector (length == length(measurementLog)) of sample times in seconds
%
%   .yLabel
%       [default: 'Pressure (kPa)']

DEFAULT_FIGURETITLE = 'Pressure log' ;
DEFAULT_YLABEL      = 'Pressure (kPa)' ;

if nargin < 1
    error( 'Insufficient arguments. Must provide measurement log vector.' ) ;
end

if nargin == 1 || isempty( Params ) 
    Params.dummy = [] ;
end

if ~myisfield( Params, 'figureTitle' ) || isempty( Params.figureTitle ) 
    Params.figureTitle = DEFAULT_FIGURETITLE ;
end

if ~myisfield( Params, 'yLabel' ) || isempty( Params.yLabel ) 
    Params.yLabel = DEFAULT_YLABEL ;
end

% ------- 
figure 

if myisfield( Params, 'sampleTimes' ) && ~isempty( Params.sampleTimes ) 
    plot( Params.sampleTimes, measurementLog, '+' ) ;
    xlabel('Time (s)');
else
    plot( measurementLog, '+' ) ;
    xlabel('Sample index');
end
    
title( Params.figureTitle ) ;
ylabel( Params.yLabel ) ;

end
% =========================================================================
function [medianMeasure] = userselectmedianmeasurement( measurementLog )
% USERSELECTMEDIANMEASUREMENT
%
%   medianMeasure = USERSELECTMEDIANMEASUREMENT( measurementLog ) 
%
%   Plots measurementLog and the user selects START and END (apnea) indices
%   over which to calculate the median. The median measurement is superposed
%   over the measurementLog graph and the user is asked if the result is 
%   satisfactory (or redo).

isUserSatisfied = false ;

while ~isUserSatisfied

    gcf ;
    plot( measurementLog, '+' ) ;
    title( 'Measure Log' ) ;
    
    xlabel('Sample index');
    ylabel('Amplitude');
    
    trainingFrameStartIndex = ...
        input( ['Identify sample index corresponding to beginning of training frame ' ...
            '([Enter] selects sample 1): '] ) ;
    
    if isempty(trainingFrameStartIndex)
        trainingFrameStartIndex = 1;
    end

    trainingFrameEndIndex = ...
        input( ['Identify sample index corresponding to end of training frame ' ...
            '([Enter] selects the last recorded sample): '] ) ;

    if isempty(trainingFrameEndIndex)
       medianMeasure = ...
           median( measurementLog( trainingFrameStartIndex : end ) ) ;
    else
       medianMeasure = ...
           median( measurementLog( trainingFrameStartIndex : trainingFrameEndIndex ) ) ;
    end

    gcf; 
    plot( measurementLog, '+' );
    hold on;
    plot( medianMeasure*ones( size( measurementLog ) ) ) ;
    title( 'Measure Log' ) ;
    xlabel('Sample index');
    ylabel('Amplitude');
    legend('Measure log','Median measurement over given interval');    
    hold off;

    response = input(['Is the current median estimate satisfactory?' ...
        '0 to re-enter the data range; 1 (or enter) to continue: ']) ;

     if ~isempty(response)
        isUserSatisfied = logical(response) ;
     else
         isUserSatisfied = true;
     end

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

