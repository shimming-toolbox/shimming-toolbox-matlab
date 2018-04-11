classdef ShimOptMartinos7TSiemens < ShimOpt
%SHIMOPTMARTINOS7TSIEMENS - Shim Optimization for Siemens shims on 7T (Bay 5) at Martinos
% 
% .......
% 
% Usage
%
% Shim = ShimOptMartinos7TSiemens( Params )
% 
% Defaults 
% 
% Params.pathToShimReferenceMaps =
%   '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007.mat'
%
% Params.ProbeSpecs = [] ;
%   .ProbeSpecs is a parameters struct for ProbeTracking(). See HELP
%   ProbeTracking() for more information.
%
%   Shim contains fields
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info Re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%       .Field
%           Object of type MaRdI pertaining to field distribution to be shimmed
%
%       .Model
%           .currents  
%               Optimal shim current vector (i)
%               [units A]
%           .field     
%               Optimal shim field from projection of i onto reference maps (Ai)
%               [units Hz]
%           .couplingCoefficients
%               For realtime shimming, relates vector relating field to pressure
%               [units Hz/Pa]
%           .dcCurrentsOffsets
%               For realtime shimming, vector of "y-intercept" currents 
%               (i.e. currents for pressure = 0 Pa)
%               [units A]
%
%       .Probe
%           Object of type ProbeTracking
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
% ShimOpt is a MaRdI subclass [ShimOpt < MaRdI]
%     
% =========================================================================
% Updated::20170213::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
% .....
% OPTIMIZESHIMCURRENTS()
%   Run (unconstrained) Cg solver;
%   If solution achievable given system constraints, return solution;
%   Else, run fmincon given constraints & return that solution instead;
%
% =========================================================================

% properties % defined in parent class ShimOpt 
    % Field ; % object of type MaRdI
    % Model ;
    % Probe ; % object of type ProbeTracking
% end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOptMartinos7TSiemens( Params )
%SHIMOPT - Shim Optimization

% dbstop in ShimOptMartinos7TSiemens at 93;

Shim.Field = [ ] ;       
Shim.Model = [ ] ; 
% Shim.Probe = ProbeTracking( Params.ProbeSpecs )  ; 
Shim.Probe = [ ]  ; 

DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/Martinos7TReferenceMaps20170216_fit.mat';
DEFAULT_PROBESPECS = [] ;

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'ProbeSpecs' ) || isempty(Params.ProbeSpecs)
   Params.ProbeSpecs = DEFAULT_PROBESPECS ;
end

ShimUse.display(['\n Preparing for shim ...  \n\n'...
        'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

Tmp = load( Params.pathToShimReferenceMaps ) ;

%%-----
% dB/dI linear 'Current-to-Field' operator
Shim.img              = Tmp.Shim.img ;
Shim.Hdr              = Tmp.Shim.Hdr ;



end
% =========================================================================
function Shim = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Shim = OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]

Specs = ShimSpecsMartinos7TSiemens();

if nargin < 1
    error('Function requires at least 1 argument of type ShimOpt')
elseif nargin == 1
    Params.dummy = [];
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end

% Params for conjugate-gradient optimization
CgParams.tolerance     = 1E-6 ;
CgParams.maxIterations = 10000 ;    

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

b = M*(-Shim.Field.img(:)) ;

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*b, ... % effective solution vector
                            zeros( [Specs.Amp.nActiveChannels 1] ), ... % initial model (currents) guess
                            CgParams ) 

pseudoInvCurrents = Shim.Model.currents
% ------- 
% TODO: CHECK SOL. VECTOR FOR CONDITIONS 
% + allow more optional params for optimization 
isCurrentSolutionOk = false ;

if ~isCurrentSolutionOk

    Options = optimset(...
        'DerivativeCheck','off',...
        'GradObj','on',...
        'Display', 'off',... %'iter-detailed',...
        'MaxFunEvals',36000,...
        'TolX',1e-11,...
        'TolCon',1E-8);

    A=M*A;

    tic
    [Shim.Model.currents] = fmincon( ...
        @shimcost,...
        Shim.Model.currents,...
        [],...
        [],...
        [],...
        [],...
        -Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
        Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
        [],...
        Options);
    toc

end

function [f, df] = shimcost( currents )

     y=A*currents - b;
     f=y'*y;
     df=2*A'*y;

end

% function [C, Ceq] = first_order_norm( currents )
% % C(x) <= 0
% % (e.g. x = currents)
%
%     C   = 0; 
%     Ceq = [];
%     waterLevel = 1E-8; % small number (relative to |currents|) for stability
%
%     % split currents up into banks
%     % -------
%     % bank 0
%     i0 = X0*currents;
%
%     % -------
%     % bank 1
%     i1 = X1*currents;
%
%     % -------
%     % bank 2
%     i2 = X2*currents;
%
%     % -------
%     % bank 3
%     i3 = X3*currents;
%
%     % Overall abs current cannot exceed Specs.Amp.maxCurrentPerBank (e.g. 20 A)
%     % This condition is redundant given the following 2 on pos/neg currents
%     C(1) = sum( abs(i0) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
%     C(2) = sum( abs(i1) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
%     C(3) = sum( abs(i2) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
%     C(4) = sum( abs(i3) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
%
%     % pos. current cannot exceed Specs.Amp.maxCurrentPerRail (e.g. + 10 A) 
%     C(5) = abs(sum( ((i0>0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%     C(6) = abs(sum( ((i1>0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%     C(7) = abs(sum( ((i2>0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%     C(8) = abs(sum( ((i3>0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%
%     % neg. current cannot be below Specs.Amp.maxCurrentPerRail (e.g. - 10 A) 
%     C(9)  = abs(sum( ((i0<0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%     C(10) = abs(sum( ((i1<0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%     C(11) = abs(sum( ((i2<0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%     C(12) = abs(sum( ((i3<0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
%
% end
    
Shim = Shim.setforwardmodelfield ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Access=protected)
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

end
% =========================================================================
% =========================================================================

end



