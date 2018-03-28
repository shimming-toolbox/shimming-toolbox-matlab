classdef ShimOptRri < ShimOpt
%SHIMOPTRRI - Shim Optimization for RRI 24 channel array (aka spine shim)
%
% ShimOptRri is a ShimOpt subclass 
%     
% =========================================================================
% Updated::20180325::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
% .....
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
function Shim = ShimOptRri( Params, Field )
%SHIMOPTRRI - Shim Optimization

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOptRri.assigndefaultparameters( Params )


ShimUse.customdisplay(['\n Preparing for shim ...  \n\n'...
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

Shim.Aux = ShimOpt_IUGM_Prisma_fit( ) ; % Params input??

end
% =========================================================================
function [currents] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESTATICSHIMCURRENTS 
%
% currents  = OPTIMIZESHIMCURRENTS( Shim, Params )
% currents  = OPTIMIZESHIMCURRENTS( Shim, Params, FieldExpired )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]

Specs = ShimSpecsRri();

DEFAULT_REGULARIZATIONPARAMETER     = 0;
DEFAULT_ISRETURNINGPSEUDOINVERSE    = 0;

if nargin < 2 
    Params.dummy = [];
end

% TODO (if needed): define RRI system-specific Params

if nargin < 3
    currents = @optimizeshimcurrents.ShimOpt( Shim, Specs, Params, @checknonlinearconstraints ) ;
end

function [C, Ceq] = checknonlinearconstraints( currents )
%CHECKNONLINEARCONSTRAINTS 
%
% Check current solution satisfies nonlinear system constraints
% 
% i.e. this is the C(x) function in FMINCON (see DOC)
%
% C(x) <= 0
%
% (e.g. x = currents)
    
    [X0, X1, X2, X3] = ShimComRri.getchanneltobankmatrices( ) ;

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
    
    % check on abs current per channel
    C((end+1):(end+length(currents))) = abs(currents) - Params.maxCurrentPerChannel ;
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

DEFAULT_INSTITUTIONNAME = 'IUGM' ;
DEFAULT_STATIONNAME     = 'MRC35049' ;

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_PROBESPECS ;
end

if ~myisfield( Params, 'isInterpolatingReferenceMaps' ) || isempty(Params.isInterpolatingReferenceMaps)
   Params.isInterpolatingReferenceMaps = DEFAULT_ISINTERPOLATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'InstitutionName' ) || isempty(Params.InstitutionName)
   Params.InstitutionName = DEFAULT_INSTITUTIONNAME ;
end

if ~myisfield( Params, 'StationName' ) || isempty(Params.StationName)
   Params.StationName = DEFAULT_STATIONNAME ;
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

