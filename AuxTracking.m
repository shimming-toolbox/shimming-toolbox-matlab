classdef (Abstract) AuxTracking < AuxTracked 
% AUXTRACKING - Respiration/Field tracking for real-time shimming 
%
% Aux = TRACKING(  )
%
%   Aux contains fields
%
%       .Data 
%           .p 
%               data measurements themselves (e.g. pressure or navigator phase)
%           .t
%               (not implemented, but could be the sample times of .p)
%
%       .Specs
%           .dt 
%               sampling interval in units of ms?
%
% .......
%
%   Description
%
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%    AuxTracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%
% =========================================================================
% Updated::20180302::ryan.topfer@polymtl.ca
% =========================================================================

% *** TODO 
%
% ..... 
%
% =========================================================================

properties   
    % Data ;
    % Specs ;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Aux = AuxTracking( Specs )
%TRACKING  

if nargin < 1
    Specs = [] ;
end

Aux.Specs  = Specs ;
Aux.Data   = [] ; 
Aux.Data.p = [] ; 

end    
% =========================================================================
function [predictedMeasurement] = predictmeasurement( Aux, delay, order )
% PREDICTMEASUREMENT
% 
% pPredicted = PREDICTMEASUREMENT( Aux, delay, p0, dpdt ) 
% pPredicted = PREDICTMEASUREMENT( Aux, delay, p0, dpdt, d2pdt2 ) 
%
% delay
%   [units: ms]


if nargin < 3
    predictedMeasurement = Aux.Data.p( end ) ; % 0th order prediction 
    return;
elseif nargin == 3
    order = 0;
elseif nargin == 4
    order = 1;
elseif nargin == 5
    order = 2;
end

predictedMeasurement = p0 ; % 0th order term

if order >= 1  
    
    % 1st order Taylor expansion 
    predictedMeasurement = predictedMeasurement + delay*dpdt ;

    if order == 2
        predictedMeasurement = predictedMeasurement + ((delay^2)/2)*d2dp2 ;
    end 

end

end
% =========================================================================

% =========================================================================
% =========================================================================
end
% =========================================================================
% =========================================================================
methods(Abstract)
% =========================================================================
% =========================================================================
[AuxInert] = copyinert( Aux )
%COPYINERT  
% 
% Takes an AuxTracking object (i.e. possibly associated with an open serial port
% or something, and returns an inert AuxTracked object with copies of
%   Aux.Data
%   Aux.Specs
% etc.

% =========================================================================
[isAuxTracking] = begintracking( Aux )
%BEGINTRACKING 
% 
% (e.g. open com port)
%
% Returns true if successful. 

% =========================================================================
[] = stoptracking( Aux )
%STOPTRACKING 
%
% (e.g. close com port)

% =========================================================================
[p] = getupdate( Aux )
%GETUPDATE 
%
% Read in a single new tracking measurement (p) (e.g. from open com port)
%
% p = GETUPDATE( Aux )


% =========================================================================
end


methods(Static)
% =========================================================================

end
% =========================================================================
% =========================================================================

end
