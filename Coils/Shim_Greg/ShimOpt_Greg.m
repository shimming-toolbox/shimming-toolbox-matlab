classdef ShimOpt_Greg < ShimOpt
%SHIMOPT_GREG - Shim Optimization for Ac/Dc 8-channel array (c-spine shim)
%
% =========================================================================


% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt_Greg( varargin )

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;

[ Field, Params, Specs] = ShimOpt.parseinput( varargin ) ;

if isempty(Specs)
    Shim.System.Specs = ShimSpecs_Greg() ; 
else
    Shim.System.Specs = Specs ;
end

Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

Params = ShimOpt_Greg.assigndefaultparameters( Params ) ;

if Params.isCalibratingReferenceMaps
    
    Params = ShimOpt_Greg.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_Greg.calibratereferencemaps( Params ) ;

elseif ~isempty( Params.pathToShimReferenceMaps )

   [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt.loadshimreferencemaps( Params.pathToShimReferenceMaps ) ; 
    
    Shim.Ref.img = Shim.img ;
    Shim.Ref.Hdr = Shim.Hdr ;

end

%-------
% associate host MRI 
if strcmp( Params.InstitutionName, 'IUGM' ) && strcmp( Params.StationName, 'MRC35049' ) 
    Shim.Aux = ShimOpt_IUGM_Prisma_fit(  ) ; % Params input??
end

if ~isempty( Field )
    Shim.setoriginalfield( Field ) ;
end

end
% =========================================================================
function [Corrections] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Corrections = OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%
%   .maxCorrectionPerChannel
%       [default: determined by ShimSpecs property: .Amp.maxCurrentPerChannel]
%
%   .minCorrectionPerChannel
%       [default: -.maxCorrectionPerChannel]

if nargin < 2 
    Params.dummy = [];
end

% if ~myisfield( Params, 'maxCorrectionPerChannel') || isempty( Params.maxCorrectionPerChannel ) 
%     Params.maxCorrectionPerChannel = Shim.System.Specs.Amp.maxCurrentPerChannel ; 
% end
%
% if ~myisfield( Params, 'minCorrectionPerChannel') || isempty( Params.minCorrectionPerChannel ) 
%     Params.minCorrectionPerChannel = -Params.maxCorrectionPerChannel ; 
% end

Corrections = optimizeshimcurrents@ShimOpt( Shim, Params ) ;

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

DEFAULT_ISCALIBRATINGREFERENCEMAPS = false ;
DEFAULT_PATHTOSHIMREFERENCEMAPS    = [ shimbindir() 'ShimReferenceMaps_Greg' ] ;

DEFAULT_INSTITUTIONNAME = 'IUGM' ;
DEFAULT_STATIONNAME     = 'MRC35049' ;

if ~myisfield( Params, 'isCalibratingReferenceMaps' ) || isempty(Params.isCalibratingReferenceMaps)
   Params.isCalibratingReferenceMaps = DEFAULT_ISCALIBRATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   
    if Params.isCalibratingReferenceMaps
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ shimbindir() 'ShimReferenceMaps_Greg'  ] ;
    else
        Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
    end
end

if ~myisfield( Params, 'InstitutionName' ) || isempty(Params.InstitutionName)
   Params.InstitutionName = DEFAULT_INSTITUTIONNAME ;
end

if ~myisfield( Params, 'StationName' ) || isempty(Params.StationName)
   Params.StationName = DEFAULT_STATIONNAME ;
end

end
% =========================================================================
function Params = declarecalibrationparameters( Params )
%DECLARECALIBRATIONPARAMETERS
% 
% Initializes parameters for shim reference map construction (aka shim calibration)

disp( ['Preparing shim calibration...' ] )        

%% -----

DEFAULTS.unwrapper = 'AbdulRahman_2007' ;
DEFAULTS.threshold = 0 ;

Params.dummy = [] ;
Params = assignifempty( Params, DEFAULTS ) ;
%

Specs             = ShimSpecs_Greg() ;
Params.nChannels  = Specs.Amp.nActiveChannels ;

Params.nCurrents  = 2 ;

% for multi-coil systems, entries of Params.currents need to be manually
Params.currents = zeros( Params.nChannels, Params.nCurrents ) ; 
Params.currents = repmat( [0.5 -0.5], [Params.nChannels 1] ) ; , % [units: A]

%% -----
% define image data directories:

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories    = cell( Params.nCurrents, 2, Params.nChannels ) ;


%% ------
% channel 1 : 
Params.dataLoadDirectories{1,1,1} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/33-fm2d_CH1_+0.5A/' ;
Params.dataLoadDirectories{1,2,1} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/34-fm2d_CH1_+0.5A/' ;
Params.dataLoadDirectories{2,1,1} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/35-fm2d_CH1_-0.5A/' ;
Params.dataLoadDirectories{2,2,1} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/36-fm2d_CH1_-0.5A/' ;

%% ------
% channel 2 : 
Params.dataLoadDirectories{1,1,2} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/37-fm2d_CH2_+0.5A/' ;
Params.dataLoadDirectories{1,2,2} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/38-fm2d_CH2_+0.5A/' ;
Params.dataLoadDirectories{2,1,2} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/39-fm2d_CH2_-0.5A/' ;
Params.dataLoadDirectories{2,2,2} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/40-fm2d_CH2_-0.5A/' ;

%% ------
% channel 3 : 
Params.dataLoadDirectories{1,1,3} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/41-fm2d_CH3_+0.5A/' ;
Params.dataLoadDirectories{1,2,3} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/42-fm2d_CH3_+0.5A/' ;
Params.dataLoadDirectories{2,1,3} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/43-fm2d_CH3_-0.5A/' ;
Params.dataLoadDirectories{2,2,3} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/44-fm2d_CH3_-0.5A/' ;

%% ------
% channel 4 : 
Params.dataLoadDirectories{1,1,4} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/45-fm2d_CH4_+0.5A/' ;
Params.dataLoadDirectories{1,2,4} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/46-fm2d_CH4_+0.5A/' ;
Params.dataLoadDirectories{2,1,4} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/47-fm2d_CH4_-0.5A/' ;
Params.dataLoadDirectories{2,2,4} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/48-fm2d_CH4_-0.5A/' ;

%% ------
% channel 5 : 
Params.dataLoadDirectories{1,1,5} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/49-fm2d_CH5_+0.5A/' ;
Params.dataLoadDirectories{1,2,5} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/50-fm2d_CH5_+0.5A/' ;
Params.dataLoadDirectories{2,1,5} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/55-fm2d_CH5_-0.5A/' ;
Params.dataLoadDirectories{2,2,5} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/56-fm2d_CH5_-0.5A/' ;

%% ------
% channel 6 : 
Params.dataLoadDirectories{1,1,6} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/53-fm2d_CH6_+0.5A/' ;
Params.dataLoadDirectories{1,2,6} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/54-fm2d_CH6_+0.5A/' ;
Params.dataLoadDirectories{2,1,6} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/57-fm2d_CH6_-0.5A/' ;
Params.dataLoadDirectories{2,2,6} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/58-fm2d_CH6_-0.5A/' ;

%% ------
% channel 7 : 
Params.dataLoadDirectories{1,1,7} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/59-fm2d_CH7_+0.5A/' ;
Params.dataLoadDirectories{1,2,7} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/60-fm2d_CH7_+0.5A/' ;
Params.dataLoadDirectories{2,1,7} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/61-fm2d_CH7_-0.5A/' ;
Params.dataLoadDirectories{2,2,7} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/62-fm2d_CH7_-0.5A/' ;

%% ------
% channel 8 : 
Params.dataLoadDirectories{1,1,8} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/63-fm2d_CH8_+0.5A/' ;
Params.dataLoadDirectories{1,2,8} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/64-fm2d_CH8_+0.5A/' ;
Params.dataLoadDirectories{2,1,8} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/65-fm2d_CH8_-0.5A/' ;
Params.dataLoadDirectories{2,2,8} = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/66-fm2d_CH8_-0.5A/' ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
