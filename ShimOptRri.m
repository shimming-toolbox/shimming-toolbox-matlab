classdef ShimOptRri < ShimOpt
%SHIMOPTRRI - Shim Optimization for RRI 24 channel array (aka spine shim)
%
% .......
% 
% Usage
%
% Shim = ShimOptRri( Params )
% 
% Defaults 
% 
% Params.pathToShimReferenceMaps =
%   '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007.mat'
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
% Updated::20170214::ryan.topfer@polymtl.ca
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
function Shim = ShimOptRri( Params )
%SHIMOPTRRI - Shim Optimization

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOptRri.assigndefaultparameters( Params )


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
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]

Specs = ShimSpecsRri();

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

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end
%
% if ~myisfield(Params, 'controlOffsetAndLinear')
%     Params.controlOffsetAndLinear = false;
% end

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
CgParams.tolerance     = 1E-6 ;
CgParams.maxIterations = 10000 ;    



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
        
        % augmented again :     
        % the added vector is the target inspired-expired difference-field
        % ------
        % changed 20171005 by RT :
        if ~myisfield(Params, 'targetFieldDifference') || isempty( Params.targetFieldDifference ) 
            % Previously was targetFieldDifference was a vector of zeros :
            Params.targetFieldDifference = zeros(size(M*FieldExpired.img(:))) ; 
            % riro = Shim.Field.img - FieldExpired.img ;
            % Params.targetFieldDifference = M*riro(:) ;
        end
        % -----
        % dbstop in ShimOptRri at 262
        b = [b ; Params.targetFieldDifference ];   
        
        A = [ AA; P ] ;
    else

        A = AA;
    
    end

    % M = repmat( M, [2 1] ) ;

end



% if Params.controlOffsetAndLinear
%     
%     % compute the shimming offset (Hz)
%     
%     offset = sum(b)/sum(M*ones(size(b)));
%     
%     %compute the desired linear shimming (Hz/mm)
%     
%     [X,Y,Z] = Shim.Field.getvoxelpositions();
%     
%     xx = M*X(:);
%     yy = M*Y(:);
%     zz = M*Z(:);
%     
%     dBdx = xx'*b./(xx'*xx);
%     dBdy = yy'*b./(yy'*yy);
%     dBdz = zz'*b./(zz'*zz);
%     
%     % Create and write on test file
%     
%     fid = fopen([Params.dataLoadDir datestr(now, 30) '-values_for_offline_shimming.txt'], 'w+');
%     fprintf(fid, '%f\n', offset);
%     fprintf(fid, '%f\n', dBdx);
%     fprintf(fid, '%f\n', dBdy);
%     fprintf(fid, '%f\n', dBdz);
%     fclose(fid);
% end

% dbstop in ShimOptRri at 296

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
    
    [X0,X1,X2,X3] = ShimComRri.getchanneltobankmatrices( ) ;

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
            @first_order_norm,...
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
            @first_order_norm_augmented,...
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

function [C, Ceq] = first_order_norm( currents )
% C(x) <= 0
% (e.g. x = currents)

    C   = 0; 
    Ceq = [];
    waterLevel = 1E-8; % small number (relative to |currents|) for stability

    % split currents up into banks
    % -------
    % bank 0
    i0 = X0*currents;

    % -------
    % bank 1
    i1 = X1*currents;

    % -------
    % bank 2
    i2 = X2*currents;

    % -------
    % bank 3
    i3 = X3*currents;

    % Overall abs current cannot exceed Specs.Amp.maxCurrentPerBank (e.g. 20 A)
    % This condition is redundant given the following 2 on pos/neg currents
    C(1) = sum( abs(i0) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(2) = sum( abs(i1) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(3) = sum( abs(i2) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(4) = sum( abs(i3) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 

    % pos. current cannot exceed Specs.Amp.maxCurrentPerRail (e.g. + 10 A) 
    C(5) = abs(sum( ((i0>0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(6) = abs(sum( ((i1>0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(7) = abs(sum( ((i2>0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(8) = abs(sum( ((i3>0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    
    % neg. current cannot be below Specs.Amp.maxCurrentPerRail (e.g. - 10 A) 
    C(9)  = abs(sum( ((i0<0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(10) = abs(sum( ((i1<0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(11) = abs(sum( ((i2<0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(12) = abs(sum( ((i3<0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 

end
    
function [C, Ceq] = first_order_norm_augmented( currents )
% C(x) <= 0
% (e.g. x = currents)

Ceq = [];

[currentsInspired, currentsExpired] = splitcurrentvector( currents ) ;

[Cinspired, ~] = first_order_norm( currentsInspired ) ;
[Cexpired,  ~] = first_order_norm( currentsExpired ) ;

C = [Cinspired; Cexpired];

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

% DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007.mat';
% DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170410.mat';
% DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170418.mat';
DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170706.mat';
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

