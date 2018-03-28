classdef ProbeTracked < AuxTracked

properties   
    % ;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Aux = ProbeTracked( Specs  )
%PROBE - ProbeTracked  

if nargin < 1
    Specs = [] ;
end
Specs.state = 'inert';

Aux.Data.p = [];

end    
% =========================================================================
function [AuxCopy] = copy( Aux )
%COPY  
% 
% Aux = COPY( Aux )

AuxCopy = ProbeTracked( Aux.Specs ) ; % or just allow that constructor to accept ProbeTracking as an input and copy the fields...

AuxCopy.Data = Aux.Data ;

end
% =========================================================================
function [] = delete( Aux )
%DELETE  
% DELETE( Aux )
% 

clear Aux ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [pressureLog, sampleTimes] = loadpressurelog( pressureLogFilename, sampleTimesFilename )
%LOADPRESSURELOG
% 
% Wraps to AuxTracking.loadmeasurementlog()
%
% pressureLog                = LOADPRESSURELOG( pressureLogFilename ) ;
% [pressureLog, sampleTimes] = LOADPRESSURELOG( pressureLogFilename, sampleTimesFilename )

if nargin < 1
    error( 'Insufficient arguments. Must provide full path to pressure log .bin file.' ) ;
else
    if nargin == 1
        [pressureLog] = AuxTracking.loadmeasurementlog( pressureLogFilename ) ;
    elseif nargin == 2
        [pressureLog, sampleTimes] = ...
            AuxTracking.loadmeasurementlog( pressureLogFilename, sampleTimesFilename ) ;
    end
end

end
% =========================================================================
function [] = plotpressurelog( pressureLog, Params )
%PLOTPRESSURELOG
% 
% Wraps to AuxTracking.plotmeasurementlog()
%
% PLOTPRESSURELOG( pressureLog ) ;
% PLOTPRESSURELOG( pressureLog, Params )
%
% Supported fields to Params struct
%
%   .figureTitle
%       [default: 'Pressure log']
%
%   .sampleTimes
%       vector (length == length(pressureLog)) of sample times in seconds

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

AuxTracking.plotmeasurementlog( pressureLog, Params ) ;

end
% =========================================================================
function [medianPressure] = userselectmedianpressure( pressureLog )
% USERSELECTMEDIANPRESSURE 
%
% Wraps to AuxTracking.userselectmedianmeasurement()
%
% medianPressure = USERSELECTMEDIANPRESSURE( pressureLog ) 
%
% Plots pressureLog and the user selects START and END (apnea) indices over
% which to calculate the median. The median pressure is superposed over the
% pressureLog graph and the user is asked if the result is satisfactory (or
% redo).

medianPressure = AuxTracking.userselectmedianmeasurement( pressureLog )

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

