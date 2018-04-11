classdef Aux < matlab.mixin.SetGet
% AUX - Auxiliary data/measurement(s) associated with a MaRdI instance
% 
% for now, the only valid "Aux" objects are ProbeTracking & ProbeTracked
% (the home-built respiratory probe)
% however, other possiblilities exist for future implementation:
% e.g.
% -TTL watchdog for slicewise shimming
% -resp bellows from the PMU (which may have a very different implementation to the homebuilt system)
% -some sort of image-file (e.g. socket) or measurement-data (e.g. nav) watchdog?
% -perhaps even other sorts of images might be linked here?
%
% AuxObj = Aux(  )
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
% Updated::20170407::ryan.topfer@polymtl.ca
% =========================================================================

% *** TODO 
% 
% Not sure how general the methods would be so there aren't any for now.
% ..... 
% 
% How is Aux.Specs.state managed?
%
% =========================================================================

properties   
    Data ;
    Specs ; % state = {active, inactive, inert, void}
    Params ;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Aux = Aux( Specs )
%TRACKING  

if nargin < 1
    Specs = [] ;
    Specs.state = 'void' ;
end

Aux.Specs  = Specs ;
Aux.Params = [] ;
Aux.Data   = [] ; 
Aux.Data.p = [] ; 

end    
% =========================================================================

end
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
methods(Abstract)
% =========================================================================

end
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================

end

