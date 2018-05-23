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
Shim.System.Specs    = ShimSpecsAcdc() ; 
Shim.System.currents = zeros( 1, Shim.System.Specs.Amp.nActiveChannels ) ; 

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_Acdc.assigndefaultparameters( Params ) ;

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

DEFAULT_ISOPTIMIZINGAUX          = true ;

if nargin < 2 
    Params.dummy = [];
end

if ~myisfield(Params, 'isOptimizingAux') || isempty( Params.isOptimizingAux )
    Params.isOptimizingAux = DEFAULT_ISOPTIMIZINGAUX ;
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Shim.System.Specs.Amp.maxCurrentPerChannel ; 
    if Params.isOptimizingAux
        Params.maxCurrentPerChannel = [ Params.maxCurrentPerChannel Shim.Aux.System.Specs.Amp.maxCurrentPerChannel ]; 
    end
end

if ~myisfield(Params, 'minCurrentPerChannel') || isempty( Params.minCurrentPerChannel ) 
    Params.minCurrentPerChannel = -Params.maxCurrentPerChannel ; 
    if Params.isOptimizingAux
        Params.maxCurrentPerChannel = [ Params.minCurrentPerChannel -Shim.Aux.System.Specs.Amp.maxCurrentPerChannel ]; 
    end
end

currents = optimizeshimcurrents@ShimOpt( Shim, Params, @checknonlinearconstraints ) ;

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

DEFAULT_ISCALIBRATINGREFERENCEMAPS = false ;
DEFAULT_PATHTOSHIMREFERENCEMAPS = '~/Projects/Shimming/Acdc/Calibration/data/ShimReferenceMaps_Acdc_20180517.mat';
DEFAULT_PROBESPECS              = [] ;

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

tmp = { ...
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/218-gre_field_mapping_ch0_0mA_S84_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/216-gre_field_mapping_ch0_0mA_S86_DIS3D /echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/214-gre_field_mapping_ch1_-400mA_S88_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/212-gre_field_mapping_ch1_-400mA_S90_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/210-gre_field_mapping_ch1_+400mA_S92_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/208-gre_field_mapping_ch1_+400mA_S94_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/206-gre_field_mapping_ch2_-400mA_S96_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/204-gre_field_mapping_ch2_-400mA_S98_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/203-gre_field_mapping_ch2_+400mA_S100_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/201-gre_field_mapping_ch2_+400mA_S102_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/199-gre_field_mapping_ch3_-400mA_S104_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/197-gre_field_mapping_ch3_-400mA_S106_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/195-gre_field_mapping_ch3_+400mA_S108_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/193-gre_field_mapping_ch3_+400mA_S110_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/191-gre_field_mapping_ch4_-400mA_S112_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/189-gre_field_mapping_ch4_-400mA_S114_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/187-gre_field_mapping_ch4_+400mA_S116_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/185-gre_field_mapping_ch4_+400mA_S118_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/183-gre_field_mapping_ch5_-400mA_S120_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/181-gre_field_mapping_ch5_-400mA_S122_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/179-gre_field_mapping_ch5_+400mA_S124_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/177-gre_field_mapping_ch5_+400mA_S126_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/175-gre_field_mapping_ch6_-400mA_S128_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/173-gre_field_mapping_ch6_-400mA_S130_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/171-gre_field_mapping_ch6_+400mA_S132_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/169-gre_field_mapping_ch6_+400mA_S134_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/167-gre_field_mapping_ch7_-400mA_S136_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/165-gre_field_mapping_ch7_-400mA_S138_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/163-gre_field_mapping_ch7_+400mA_S140_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/161-gre_field_mapping_ch7_+400mA_S142_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/159-gre_field_mapping_ch8_-400mA_S144_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/157-gre_field_mapping_ch8_-400mA_S146_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/155-gre_field_mapping_ch8_+400mA_S148_DIS3D/echo_7.38/' ;
        '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/153-gre_field_mapping_ch8_+400mA_S150_DIS3D/echo_7.38/' ; } ;
%
%
% tmp = { ...
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/81-gre_field_mapping_shim0_FA70_S21_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/82-gre_field_mapping_shim0_FA70_S22_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/83-gre_field_mapping_shim_ch1_-400mA_S23_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/84-gre_field_mapping_shim_ch1_-400mA_S24_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/85-gre_field_mapping_shim_ch1_+400mA_S25_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/86-gre_field_mapping_shim_ch1_+400mA_S26_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/87-gre_field_mapping_shim_ch2_-400mA_S27_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/88-gre_field_mapping_shim_ch2_-400mA_S28_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/89-gre_field_mapping_shim_ch2_+400mA_S29_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/90-gre_field_mapping_shim_ch2_+400mA_S30_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/91-gre_field_mapping_shim_ch3_-400mA_S31_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/92-gre_field_mapping_shim_ch3_-400mA_S32_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/93-gre_field_mapping_shim_ch3_+400mA_S33_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/94-gre_field_mapping_shim_ch3_+400mA_S34_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/95-gre_field_mapping_shim_ch4_-400mA_S35_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/96-gre_field_mapping_shim_ch4_-400mA_S36_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/97-gre_field_mapping_shim_ch4_+400mA_S37_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/98-gre_field_mapping_shim_ch4_+400mA_S38_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/100-gre_field_mapping_shim_ch5_-400mA_S39_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/101-gre_field_mapping_shim_ch5_-400mA_S40_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/102-gre_field_mapping_shim_ch5_+400mA_S41_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/103-gre_field_mapping_shim_ch5_+400mA_S42_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/104-gre_field_mapping_shim_ch6_-400mA_S43_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/105-gre_field_mapping_shim_ch6_-400mA_S44_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/106-gre_field_mapping_shim_ch6_+400mA_S45_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/107-gre_field_mapping_shim_ch6_+400mA_S46_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/108-gre_field_mapping_shim_ch7_-400mA_S47_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/109-gre_field_mapping_shim_ch7_-400mA_S48_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/110-gre_field_mapping_shim_ch7_+400mA_S49_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/111-gre_field_mapping_shim_ch7_+400mA_S50_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/112-gre_field_mapping_shim_ch8_-400mA_S51_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/113-gre_field_mapping_shim_ch8_-400mA_S52_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/114-gre_field_mapping_shim_ch8_+400mA_S53_DIS3D/echo_7.38';
%     '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/115-gre_field_mapping_shim_ch8_+400mA_S54_DIS3D/echo_7.38'; } ;

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
