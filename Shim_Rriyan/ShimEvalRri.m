classdef ShimEvalRri < ShimEval & ShimOptRri 
%SHIMEVAL RRI - Shim Evaluation for RRI shim
%
% .......
% 
% Usage
%
% Shim = ShimEvalRri( )
% Shim = ShimEvalRri( Params )
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
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
% ShimEvalRri is a ShimEval & ShimOptRri subclass 
%     
% =========================================================================
% Updated::20170411::ryan.topfer@polymtl.ca
% =========================================================================


% properties
%     ShimmedField;
% end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimEvalRri( Params )
%SHIMEVAL - Shim Evaluation 

if nargin < 1
    Params = [];
end

Shim = Shim@ShimOptRri( Params ) ;

Shim.ShimmedField = [] ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

