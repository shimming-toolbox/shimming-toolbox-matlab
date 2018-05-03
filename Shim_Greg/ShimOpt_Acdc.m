classdef ShimOpt_Acdc < ShimOpt
%SHIMOPT_ACDC - Shim Optimization for Ac/Dc 8 channel array (cervical spine shim)
%
% ShimOpt_Acdc is a ShimOpt subclass 
%     
% =========================================================================
% Updated::20170328::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
%
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
function Shim = ShimOpt_Acdc( Params, Field )
%SHIMOPT_ACDC - Shim Optimization

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_Acdc.assigndefaultparameters( Params ) ;

% .......
% Load shim basis if provided 
if ~isempty(Params.pathToShimReferenceMaps)
    
ShimUse.customdisplay(['\n Preparing for shim ...  \n\n'...
        'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

    % Loads .mat: struct with fields
    %
    % .img 
    %       linear dB/dI 'current-to-field' operator
    % .Hdr
    %       dicom Hdr from the reference map acquisition 
    RefMaps = load( Params.pathToShimReferenceMaps ) ; 

    %%-----
    % dB/dI linear 'Current-to-Field' operator
    Shim.img              = RefMaps.Shim.img ;
    Shim.Hdr              = RefMaps.Shim.Hdr ;

elseif ~isempty( Params.dataLoadDirectories )
   
   [ Shim.img, Shim.Hdr ] = ShimOpt.mapdbdi( Params ) ;

end


Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 

if (nargin == 2) && (~isempty(Field))
    
    Shim.setoriginalfield( Field ) ;

end

%-------
% associate host MRI 
if strcmp( Params.InstitutionName, 'IUGM' ) && strcmp( Params.StationName, 'MRC35049' ) 

    Shim.Aux = ShimOpt_IUGM_Prisma_fit( ) ; % Params input??

end

end
% =========================================================================
function [currents] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% currents = OPTIMIZESHIMCURRENTS( Shim, Params )
% [currentsInspired, currentsExpired] = OPTIMIZESHIMCURRENTS( Shim, Params, FieldExpired )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: determined by class ShimSpecsAcdc.Amp.maxCurrentPerChannel]
 

Specs = ShimSpecsAcdc();

DEFAULT_REGULARIZATIONPARAMETER     = 0;
DEFAULT_ISRETURNINGPSEUDOINVERSE    = 0;

if nargin < 2 
    Params.dummy = [];
end

% TODO (if needed): define AcDc system-specific Params

currents = optimizeshimcurrents@ShimOpt( Shim, Specs, Params, @checknonlinearconstraints ) ;

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
    
    Ceq = [];
    % check on abs current per channel
    C = abs(currents) - Params.maxCurrentPerChannel ;
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

DEFAULT_PATHTOSHIMREFERENCEMAPS = '~/Projects/Shimming/Acdc/Calibration/data/ShimReferenceMaps_Acdc_20180326.mat';
DEFAULT_PROBESPECS              = [] ;

DEFAULT_INSTITUTIONNAME = 'IUGM' ;
DEFAULT_STATIONNAME     = 'MRC35049' ;

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_PROBESPECS ;
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
