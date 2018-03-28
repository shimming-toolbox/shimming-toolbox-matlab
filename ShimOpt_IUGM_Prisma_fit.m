classdef ShimOpt_IUGM_Prisma_fit < ShimOpt
%SHIMOPTUNFPRISMA - Shim Optimization for Prisma-fit @ UNF 
%
% .......
% 
% Usage
%
% Shim = ShimOpt_IUGM_Prisma_fit( Params )
% 
% Defaults 
% 
% Params.pathToShimReferenceMaps = ...
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
% Updated::20170301::ryan.topfer@polymtl.ca
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
function Shim = ShimOpt_IUGM_Prisma_fit( Params, Field )
%SHIMOPTACDC - Shim Optimization

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_IUGM_Prisma_fit.assigndefaultparameters( Params ) ;

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

Shim.Aux = [] ;

% DICOM .Hdr.Private_0019_1013 
%
% Describes the absolute table position (of the reference maps). 
% Used in ShimOpt.interpolatetoimggrid() to account for variable table positioning,
% but for the scanner shims the field shift doesn't depend on table position, so:
Shim.Hdr.Private_0019_1013 = [] ;

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
%       [default: determined by class ShimSpecsPrisma.Amp.maxCurrentPerChannel]
 

Specs = ShimSpecsPrisma();

DEFAULT_REGULARIZATIONPARAMETER     = 0;
DEFAULT_ISRETURNINGPSEUDOINVERSE    = true; % THIS UNTIL MAX SHIM CURRENTS ARE KNOWN

if nargin < 2 
    Params.dummy = [];
end

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

DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMaps_IUGM_Prisma_fit_20180325';
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

