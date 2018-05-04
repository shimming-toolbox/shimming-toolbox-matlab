classdef ShimOpt_Acdc < ShimOpt
%SHIMOPT_ACDC - Shim Optimization for Ac/Dc 8 channel array (cervical spine shim)
%
% ShimOpt_Acdc is a ShimOpt subclass 
%     
% =========================================================================
% Updated::20180503::ryan.topfer@polymtl.ca
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

Params = ShimOpt_IUGM_Prisma_fit.assigndefaultparameters( Params ) ;

if Params.isCalibratingReferenceMaps
    
    Params = ShimOpt_Acdc.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_Acdc.calibratereferencemaps( Params ) ;

elseif ~isempty( Params.pathToShimReferenceMaps )

    [ Shim.img, Shim.Hdr ] = ShimOpt.loadshimreferencemaps( Params.pathToShimReferenceMaps )   

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
   
    if Params.isCalibratingReferenceMaps
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ '~/Projects/Shimming/Static/Calibration/Data/' ...
                        'ShimReferenceMaps_' 'Acdc_' today ] ;
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

% error('if Params is empty, then replace w/following Params :' )
% ...

Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;

Params.currents = [-0.4 0.4; %ch1, [units: A]
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; %...
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; % 
                   -0.4 0.4;] ; %ch8

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories = cell( Params.nCurrents, 2, Params.nChannels ) ;

% tmp = { ...
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/149-eld_mapping_shim0_axial_fovPhase100perc_phaseOver0perc_S76_DIS3D/echo_7.38/'  ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/151-eld_mapping_shim0_axial_fovPhase100perc_phaseOver0perc_S78_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/153-gre_field_mapping_A11_minus30_S80_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/155-gre_field_mapping_A11_minus30_S82_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/157-gre_field_mapping_A11_plus30_S84_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/159-gre_field_mapping_A11_plus30_S86_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/161-gre_field_mapping_B11_minus30_S88_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/163-gre_field_mapping_B11_minus30_S90_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/165-gre_field_mapping_B11_plus30_S92_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/167-gre_field_mapping_B11_plus30_S94_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/169-gre_field_mapping_A10_minus30_S96_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/171-gre_field_mapping_A10_minus30_S98_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/173-gre_field_mapping_A10_plus30_S101_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/175-gre_field_mapping_A10_plus30_S103_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/177-gre_field_mapping_A20_minus600_S105_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/179-gre_field_mapping_A20_minus600_S107_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/181-gre_field_mapping_A20_plus600_S109_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/183-gre_field_mapping_A20_plus600_S111_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/185-gre_field_mapping_A21_minus600_S113_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/187-gre_field_mapping_A21_minus600_S115_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/189-gre_field_mapping_A21_plus600_S117_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/191-gre_field_mapping_A21_plus600_S119_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/193-gre_field_mapping_B21_minus600_S121_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/195-gre_field_mapping_B21_minus600_S123_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/197-gre_field_mapping_B21_plus600_S125_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/199-gre_field_mapping_B21_plus600_S127_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/201-gre_field_mapping_A22_minus800_S129_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/203-gre_field_mapping_A22_minus800_S131_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/209-gre_field_mapping_A22_plus600_S137_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/211-gre_field_mapping_A22_plus600_S139_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/213-gre_field_mapping_B22_minus600_S141_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/215-gre_field_mapping_B22_minus600_S143_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/217-gre_field_mapping_B22_plus600_S145_DIS3D/echo_7.38/' ;
%     '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/219-gre_field_mapping_B22_plus600_S147_DIS3D/echo_7.38/' ; } ;

% 1st 2 directories correspond to the baseline shim 
Params.dataLoadDirectories{1,1,1} = tmp{1} ;
Params.dataLoadDirectories{1,2,1} = tmp{2} ;

nImgPerCurrent = 2 ; % = 1 mag image + 1 phase

disp( ['Preparing shim calibration...' ] )        

for iChannel = 1 : Params.nChannels
    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
    
    for iCurrent = 1 : Params.nCurrents 
        Params.dataLoadDirectories{ iCurrent, 1, iChannel + 1} = tmp{ nImgPerCurrent*(Params.nCurrents*iChannel + iCurrent) -3 } ;
        Params.dataLoadDirectories{ iCurrent, 2, iChannel + 1} = tmp{ nImgPerCurrent*(Params.nCurrents*iChannel + iCurrent) -2 } ;
    end
end


Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img/max(Mag.img(:)) > 0.01 ; % region of reliable SNR for unwrapping

Params.reliabilityMask(:,1:5,:)    = 0 ; % based on visual inspection of magnitude (there is aliasing in phase-encore dir)
Params.reliabilityMask(:,71:end,:) = 0 ; 

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize   = voxelSize ;
Params.Extension.expansionOrder = 2 ;
Params.Extension.radius     = 8 ;

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
