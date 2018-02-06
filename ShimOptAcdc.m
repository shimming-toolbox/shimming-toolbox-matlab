classdef ShimOptAcdc < ShimOpt
%SHIMOPTACDC - Shim Optimization for Ac/Dc 8 channel array (cervical spine shim)
%
% .......
% 
% Usage
%
% Shim = ShimOptAcdc( Params )
% 
% Defaults 
% 
% Params.pathToShimReferenceMaps =
%
% Params.TrackerSpecs = [] ;
%   .TrackerSpecs is a parameters struct for ProbeTracking(). See HELP
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
%       .Tracker
%           Object of type TrackerTracking
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    FieldEval
%    Tracking
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
% Updated::20171025::ryan.topfer@polymtl.ca
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
    % Tracker ; % object of type ProbeTracking
% end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOptAcdc( Params )
%SHIMOPTRRI - Shim Optimization

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOptAcdc.assigndefaultparameters( Params )


ShimUse.display(['\n Preparing for shim ...  \n\n'...
        'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

% Loads .mat containing Shim struct
% which has fields
% Shim.img              - linear dB/dI 'current-to-field' opterator
% Shim.Hdr              - defines info like voxel locations 
% Shim.Hdr.MaskingImage - defines spatial support of reference maps
RefMaps = load( Params.pathToShimReferenceMaps ) ; % load shim ref maps

%%-----
% dB/dI linear 'Current-to-Field' operator
Shim.img              = RefMaps.Shim.img ;
Shim.Hdr              = RefMaps.Shim.Hdr ;

% if myisfield( SpineShim, 'mask' )
%     Shim.Hdr.MaskingImage = SpineShim.mask ;
% end
%
% load( Params.pathToShimReferenceMaps ) ;

Shim.Field = [ ] ;       
Shim.Model = [ ] ; 
Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 

if (nargin == 2) && (~isempty(Field))
    
    Shim = Shim.setoriginalfield( Field ) ;
    
    if Params.isInterpolatingReferenceMaps
        
        Shim = interpolatetoimggrid( Shim, Field ) ;
        Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

    end

else
    Shim.Field = [ ] ; % user can assign later via Shim.setoriginalfield() 

end

end
% =========================================================================
function [currents, currentsExpired] = optimizeshimcurrents( Shim, Params, FieldExpired )
%OPTIMIZESHIMCURRENTS 
%
% currents = OPTIMIZESHIMCURRENTS( Shim, Params )
% [currentsInspired, currentsExpired] = OPTIMIZESHIMCURRENTS( Shim, Params, FieldExpired )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: determined by class ShimSpecsAcdc.Amp.maxCurrentPerChannel]
% 
% deprecated???
%   .maxVoltagePerChannel
%       [default: determined by class ShimSpecsAcdc.Amp.maxVoltagePerChannel]

Specs = ShimSpecsAcdc();

DEFAULT_REGULARIZATIONPARAMETER     = 0;
DEFAULT_ISRETURNINGPSEUDOINVERSE    = 0;

if nargin < 1
    error('Function requires at least 1 argument of type ShimOpt')
elseif nargin == 1
    Params.dummy = [];
    DEFAULT_ISSOLVINGAUGMENTEDSYSTEM    = false;
    DEFAULT_ISPENALIZINGFIELDDIFFERENCE = false;
elseif nargin == 2
    DEFAULT_ISSOLVINGAUGMENTEDSYSTEM    = false;
    DEFAULT_ISPENALIZINGFIELDDIFFERENCE = false;
elseif nargin == 3
    DEFAULT_ISSOLVINGAUGMENTEDSYSTEM    = true;
    DEFAULT_ISPENALIZINGFIELDDIFFERENCE = true;
end

if ~myisfield(Params, 'isReturningPseudoInverse') || isempty( Params.isReturningPseudoInverse ) 
    Params.isReturningPseudoInverse = DEFAULT_ISRETURNINGPSEUDOINVERSE ; 
end

if ~myisfield(Params, 'maxVoltagePerChannel') || isempty( Params.maxVoltagePerChannel ) % deprecated?
    Params.maxVoltagePerChannel = Specs.Amp.maxVoltagePerChannel ; 
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end

if ~myisfield(Params, 'isSolvingAugmentedSystem') || isempty( Params.isSolvingAugmentedSystem )
    Params.isSolvingAugmentedSystem = DEFAULT_ISSOLVINGAUGMENTEDSYSTEM ;
end

if ~myisfield(Params, 'isPenalizingFieldDifference') || isempty( Params.isPenalizingFieldDifference ) 
    Params.isPenalizingFieldDifference = DEFAULT_ISPENALIZINGFIELDDIFFERENCE ;
end

if ~myisfield(Params, 'regularizationParameter') || isempty( Params.regularizationParameter ) 
    Params.regularizationParameter = DEFAULT_REGULARIZATIONPARAMETER ;
end
 

% Params for conjugate-gradient optimization
CgParams.tolerance     = 1E-10 ;
CgParams.maxIterations = 100000 ;    

nImg = numel( Shim.Field.img(:) ) ; % number of voxels

% -------
% define matrix of data-weighting coefficients : W
if ~myisfield( Params, 'dataWeights' ) || isempty( Params.dataWeights ) 

    W = speye( nImg, nImg ) ;

else

    assert( numel( Params.dataWeights ) == nImg ) 

    if ( size( Params.dataWeights, 1 ) ~= nImg ) || ( size( Params.dataWeights, 2) ~= nImg )
        
        W = spdiags( Params.dataWeights(:), 0, nImg, nImg ) ;

    end

end

M  = Shim.gettruncationoperator*W ;


A  = M*Shim.getshimoperator ; % masked current-to-field operator
MA = A;

b = M*(-Shim.Field.img(:)) ;

solutionVectorLength = Specs.Amp.nActiveChannels ;

if Params.isSolvingAugmentedSystem
   
   
    % stacked/augmented solution vector length 
    solutionVectorLength = 2*Specs.Amp.nActiveChannels ;

    % augmented data vector : inspired field, expired field, vertically concatenated
    b = [b ; M*(-FieldExpired.img(:)) ] ;
    
    % augmented matrix operator
    AA = [A zeros( size(A) ) ;  zeros( size(A) ) A ] ;
    
    if Params.isPenalizingFieldDifference
        
        % field-difference penalizer
        P  = sqrt(Params.regularizationParameter) * [MA -MA] ;
        % % current-difference penalizer
        % P =  sqrt(Params.regularizationParameter) *[eye(24) -eye(24)] ;

        % augmented again :     
        % the added vector is the target inspired-expired difference-field
        % ------
        % changed 20171005 by RT :
        if ~myisfield(Params, 'targetFieldDifference') || isempty( Params.targetFieldDifference ) 
            % Previously was targetDifference was a vector of zeros :
            Params.targetDifference = zeros(size(M*FieldExpired.img(:))) ; 
            % riro = Shim.Field.img - FieldExpired.img ;
            % Params.targetDifference = M*riro(:) ;
            % Params.targetDifference = zeros(size(P*ones(2*24,1))) ;
        end
        
        b = [b ; Params.targetDifference ];   
        
        A = [ AA; P ] ;
    else

        A = AA;
    
    end

end

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currents = cgls( A'*A, ... % least squares operator
                            A'*b, ... % effective solution vector
                            zeros( [solutionVectorLength 1] ), ... % initial model (currents) guess
                           CgParams ) ;

currents = Shim.Model.currents ;

if Params.isSolvingAugmentedSystem

    [currents, currentsExpired] = splitcurrentvector( currents ) ;

end


% ------- 
% TODO: CHECK SOL. VECTOR FOR CONDITIONS 
% + allow more optional params for optimization 
isCurrentSolutionOk = false ;

if ~Params.isReturningPseudoInverse && ~isCurrentSolutionOk
    
    Options = optimset(...
        'DerivativeCheck','off',...
        'GradObj','on',...
        'Display', 'off',... %'iter-detailed',...
        'MaxFunEvals',36000,...
        'TolX',1e-11,...
        'TolCon',1E-8);

    tic
    if ~Params.isSolvingAugmentedSystem
   
        [currents] = fmincon( ...
            @shimcost,...
            zeros( solutionVectorLength, 1),...
            [],...
            [],...
            [],...
            [],...
            -Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            [],...
            Options);

        Shim.Model.currents = currents ;

    else

        [currents] = fmincon( ...
            @shimcost,...
            zeros( solutionVectorLength, 1),...
            [],...
            [],...
            [],...
            [],...
            -Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            [],...
            Options);
        
        [currents, currentsExpired] = splitcurrentvector( currents ) ;

    end

    toc
    
end


function [f, df] = shimcost( currents )
    
     y = A*currents - b;
     f = y'*y;
     df = 2*A'*y;
    
end


function [currentsInspired, currentsExpired] = splitcurrentvector( currents ) 
%SPLITCURRENTVECTOR
%
% De-concatenates vertically stacked/agumented current vector
    
    assert( length(currents) == 2*Specs.Amp.nActiveChannels ) ;
    currentsInspired = currents(1:Specs.Amp.nActiveChannels) ;
    currentsExpired  = currents((Specs.Amp.nActiveChannels+1):end) ;

end


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
methods(Static=true, Hidden=true)
% =========================================================================
function  [ Params ] = assigndefaultparameters( Params )
%ASSIGNDEFAULTPARAMETERS  
% 
% Params = ASSIGNDEFAULTPARAMETERS( Params )
% 
% Add default parameters fields to Params without replacing values (unless empty)
%
% DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
% DEFAULT_PROBESPECS = [] ;
%
% DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;

DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ancha_admin/Documents/Acdc/Calibration/data/AcdcReferenceMaps20171107.mat';
DEFAULT_PROBESPECS              = [] ;

DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_PROBESPECS ;
end

if ~myisfield( Params, 'isInterpolatingReferenceMaps' ) || isempty(Params.isInterpolatingReferenceMaps)
   Params.isInterpolatingReferenceMaps = DEFAULT_ISINTERPOLATINGREFERENCEMAPS ;
end

end
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