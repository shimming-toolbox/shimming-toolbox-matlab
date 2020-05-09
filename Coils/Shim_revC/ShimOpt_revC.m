classdef ShimOpt_revC < ShimOpt
%ShimOpt_revC - Shim Optimization for Ac/Dc 8 channel array (cervical spine shim)
%
% ShimOpt_revC is a ShimOpt subclass 
%     
% =========================================================================
% Updated::20180726::ryan.topfer@polymtl.ca
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
function Shim = ShimOpt_revC( Params, Field )
%ShimOpt_revC - Shim Optimization

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;
Shim.System.Specs    = ShimSpecs_revC() ; 
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_revC.assigndefaultparameters( Params ) ;

if Params.isCalibratingReferenceMaps
    
    Params = ShimOpt_revC.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_revC.calibratereferencemaps( Params ) ;

elseif ~isempty( Params.pathToShimReferenceMaps )

    [ Shim.img, Shim.Hdr ] = ShimOpt.loadshimreferencemaps( Params.pathToShimReferenceMaps )   
    
    Shim.Ref.img = Shim.img ;
    Shim.Ref.Hdr = Shim.Hdr ;

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
function [Corrections] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Corrections = OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%
%   .maxCorrectionPerChannel
%       [default: determined by ShimSpecs_revC property: .Amp.maxCurrentPerChannel]
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
        Params.pathToShimReferenceMaps = [ '~/Projects/Shimming/Static/Calibration/Data/' ...
                        'ShimReferenceMaps_' 'Greg_' today ] ;
    else
        Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
    end
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
function Params = declarecalibrationparameters( Params )
%DECLARECALIBRATIONPARAMETERS
% 
% Initializes parameters for shim reference map construction (aka shim calibration)

Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;
Params.nEchoes    = 5 ; 

Params.currents = repmat( [0.5 -0.5], [Params.nChannels 1] ) ; , % [units: A]

% magnitude threshold to define phase unwrapping region
Params.threshold       = 0.001 ;

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

tmp = { ...
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/33-fm2d_CH1_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/34-fm2d_CH1_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/35-fm2d_CH1_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/36-fm2d_CH1_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/37-fm2d_CH2_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/38-fm2d_CH2_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/39-fm2d_CH2_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/40-fm2d_CH2_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/41-fm2d_CH3_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/42-fm2d_CH3_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/43-fm2d_CH3_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/44-fm2d_CH3_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/45-fm2d_CH4_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/46-fm2d_CH4_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/47-fm2d_CH4_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/48-fm2d_CH4_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/49-fm2d_CH5_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/50-fm2d_CH5_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/55-fm2d_CH5_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/56-fm2d_CH5_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/53-fm2d_CH6_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/54-fm2d_CH6_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/57-fm2d_CH6_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/58-fm2d_CH6_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/59-fm2d_CH7_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/60-fm2d_CH7_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/61-fm2d_CH7_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/62-fm2d_CH7_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/63-fm2d_CH8_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/64-fm2d_CH8_+0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/65-fm2d_CH8_-0.5A/' ;
    '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_60p/66-fm2d_CH8_-0.5A/' ; } ;

Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

for iChannel = 1 : Params.nChannels
    for iCurrent = 1 : Params.nCurrents
        for imgType = 1 : 2 % 1=mag, 2=phase

        iDir = (iChannel-1)*(Params.nCurrents*2) + 1 ;

        dicomSubDirs  = dir( [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo*/'] ) ;
        nDicomSubDirs = length( dicomSubDirs ) ;
        assert( nDicomSubDirs == 5 )
        
        Params.dataLoadDirectories{ 1, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_2.53/' ] ;
        Params.dataLoadDirectories{ 2, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_6.8/' ] ;
        Params.dataLoadDirectories{ 3, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_12.04/' ] ;
        Params.dataLoadDirectories{ 4, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_17.28/' ] ;
        Params.dataLoadDirectories{ 5, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_22.52/' ] ;

        end
    end 
end

Params.unwrapper = 'AbdulRahman_2007' ;        

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
