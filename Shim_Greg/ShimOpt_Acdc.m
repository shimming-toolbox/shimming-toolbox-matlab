classdef ShimOpt_Acdc < ShimOpt
%SHIMOPT_ACDC - Shim Optimization for Ac/Dc 8 channel array (cervical spine shim)
%
% ShimOpt_Acdc is a ShimOpt subclass 
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
function Shim = ShimOpt_Acdc( Params, Field )
%SHIMOPT_ACDC - Shim Optimization

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;
Shim.System.Specs    = ShimSpecsAcdc() ; 
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_Acdc.assigndefaultparameters( Params ) ;

if Params.isCalibratingReferenceMaps
    
    Params = ShimOpt_Acdc.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_Acdc.calibratereferencemaps( Params ) ;

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
%       [default: determined by ShimSpecsAcdc property: .Amp.maxCurrentPerChannel]
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
DEFAULT_PATHTOSHIMREFERENCEMAPS = '~/Projects/Shimming/Acdc/Calibration/data/ShimReferenceMaps_Acdc_20180606.mat';
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
Params.nEchoes    = 2 ;

Params.currents = [0 0.2; %ch1, [units: A]
                   0 0.2; %
                   0 0.2; %
                   0 0.2; %...
                   0 0.2; %
                   0 0.2; %
                   0 0.2; % 
                   0 0.2;] ; %ch8

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

tmp = { ...
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/28-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/28-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/29-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/29-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/26-fm3d_shimON_ch1_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/26-fm3d_shimON_ch1_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/27-fm3d_shimON_ch1_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/27-fm3d_shimON_ch1_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/28-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/28-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/29-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/29-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/30-fm3d_shimON_ch2_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/30-fm3d_shimON_ch2_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/31-fm3d_shimON_ch2_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/31-fm3d_shimON_ch2_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/34-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/34-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/35-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/35-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/32-fm3d_shimON_ch3_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/32-fm3d_shimON_ch3_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/33-fm3d_shimON_ch3_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/33-fm3d_shimON_ch3_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/34-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/34-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/35-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/35-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/36-fm3d_shimON_ch4_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/36-fm3d_shimON_ch4_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/37-fm3d_shimON_ch4_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/37-fm3d_shimON_ch4_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/40-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/40-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/41-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/41-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/38-fm3d_shimON_ch5_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/38-fm3d_shimON_ch5_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/39-fm3d_shimON_ch5_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/39-fm3d_shimON_ch5_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/40-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/40-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/41-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/41-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/42-fm3d_shimON_ch6_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/42-fm3d_shimON_ch6_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/43-fm3d_shimON_ch6_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/43-fm3d_shimON_ch6_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/46-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/46-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/47-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/47-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/44-fm3d_shimON_ch7_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/44-fm3d_shimON_ch7_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/45-fm3d_shimON_ch7_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/45-fm3d_shimON_ch7_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/46-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/46-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/47-fm3d_shimOFF/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/47-fm3d_shimOFF/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/48-fm3d_shimON_ch8_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/48-fm3d_shimON_ch8_200mA/echo_10.21';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/49-fm3d_shimON_ch8_200mA/echo_4.26';
    '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_32p/49-fm3d_shimON_ch8_200mA/echo_10.21'; } ;

% 1st 4 directories correspond to the baseline shim 
Params.dataLoadDirectories{1,1,1} = tmp{1} ;
Params.dataLoadDirectories{2,1,1} = tmp{2} ;
Params.dataLoadDirectories{1,2,1} = tmp{3} ;
Params.dataLoadDirectories{2,2,1} = tmp{4} ;

nImgPerCurrent = 4 ; % = 2 mag image + 2 phase

disp( ['Preparing shim calibration...' ] )        

for iChannel = 1 : Params.nChannels

    for iCurrent = 1 : Params.nCurrents 
            
            iImg = nImgPerCurrent*(Params.nCurrents*iChannel + iCurrent - 3 ) + 1 ;
             
            % mag
            Params.dataLoadDirectories{ 1, 1, iCurrent, iChannel + 1} = tmp{ iImg } ;
            Params.dataLoadDirectories{ 2, 1, iCurrent, iChannel + 1} = tmp{ iImg + 1 } ;
            % phase
            Params.dataLoadDirectories{ 1, 2, iCurrent, iChannel + 1} = tmp{ iImg + 2 } ;
            Params.dataLoadDirectories{ 2, 2, iCurrent, iChannel + 1} = tmp{ iImg + 3 } ;
    end
end

    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/104-gre_fm_shimOFF_S67_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/105-gre_fm_shimOFF_S68_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/106-gre_fm_ch1_-400mA_S69_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/107-gre_fm_ch1_-400mA_S70_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/108-gre_fm_ch1_+400mA_S71_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/109-gre_fm_ch1_+400mA_S72_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/110-gre_fm_ch2_-400mA_S73_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/111-gre_fm_ch2_-400mA_S74_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/112-gre_fm_ch2_+400mA_S75_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/113-gre_fm_ch2_+400mA_S76_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/114-gre_fm_ch3_-400mA_S77_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/115-gre_fm_ch3_-400mA_S78_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/116-gre_fm_ch3_+400mA_S79_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/117-gre_fm_ch3_+400mA_S80_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/118-gre_fm_ch4_-400mA_S81_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/119-gre_fm_ch4_-400mA_S82_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/120-gre_fm_ch4_+400mA_S83_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/121-gre_fm_ch4_+400mA_S84_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/122-gre_fm_ch5_-400mA_S85_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/123-gre_fm_ch5_-400mA_S86_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/124-gre_fm_ch5_+400mA_S87_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/125-gre_fm_ch5_+400mA_S88_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/126-gre_fm_ch6_-400mA_S89_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/127-gre_fm_ch6_-400mA_S90_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/128-gre_fm_ch6_+400mA_S91_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/129-gre_fm_ch6_+400mA_S92_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/132-gre_fm_ch7_-400mA_S95_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/133-gre_fm_ch7_-400mA_S96_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/134-gre_fm_ch7_+400mA_S97_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/135-gre_fm_ch7_+400mA_S98_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/136-gre_fm_ch8_-400mA_S100_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/137-gre_fm_ch8_-400mA_S101_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/138-gre_fm_ch8_+400mA_S102_DIS3D/echo_7.38/' ;
    % '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/139-gre_fm_ch8_+400mA_S103_DIS3D/echo_7.38/' ; } ;

% Params.currents = [-0.4 0.4; %ch1, [units: A]
%                    -0.4 0.4; %
%                    -0.4 0.4; %
%                    -0.4 0.4; %...
%                    -0.4 0.4; %
%                    -0.4 0.4; %
%                    -0.4 0.4; % 
%                    -0.4 0.4;] ; %ch8
%
% % 2 columns: [ MAG | PHASE ] ;
% Params.dataLoadDirectories = cell( Params.nCurrents, 2, Params.nChannels ) ;
%
% tmp = { ...
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/104-gre_fm_shimOFF_S67_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/105-gre_fm_shimOFF_S68_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/106-gre_fm_ch1_-400mA_S69_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/107-gre_fm_ch1_-400mA_S70_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/108-gre_fm_ch1_+400mA_S71_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/109-gre_fm_ch1_+400mA_S72_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/110-gre_fm_ch2_-400mA_S73_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/111-gre_fm_ch2_-400mA_S74_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/112-gre_fm_ch2_+400mA_S75_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/113-gre_fm_ch2_+400mA_S76_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/114-gre_fm_ch3_-400mA_S77_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/115-gre_fm_ch3_-400mA_S78_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/116-gre_fm_ch3_+400mA_S79_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/117-gre_fm_ch3_+400mA_S80_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/118-gre_fm_ch4_-400mA_S81_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/119-gre_fm_ch4_-400mA_S82_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/120-gre_fm_ch4_+400mA_S83_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/121-gre_fm_ch4_+400mA_S84_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/122-gre_fm_ch5_-400mA_S85_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/123-gre_fm_ch5_-400mA_S86_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/124-gre_fm_ch5_+400mA_S87_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/125-gre_fm_ch5_+400mA_S88_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/126-gre_fm_ch6_-400mA_S89_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/127-gre_fm_ch6_-400mA_S90_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/128-gre_fm_ch6_+400mA_S91_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/129-gre_fm_ch6_+400mA_S92_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/132-gre_fm_ch7_-400mA_S95_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/133-gre_fm_ch7_-400mA_S96_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/134-gre_fm_ch7_+400mA_S97_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/135-gre_fm_ch7_+400mA_S98_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/136-gre_fm_ch8_-400mA_S100_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/137-gre_fm_ch8_-400mA_S101_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/138-gre_fm_ch8_+400mA_S102_DIS3D/echo_7.38/' ;
%     '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_30p/139-gre_fm_ch8_+400mA_S103_DIS3D/echo_7.38/' ; } ;
%
% tmp = { ...
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/218-gre_field_mapping_ch0_0mA_S84_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/216-gre_field_mapping_ch0_0mA_S86_DIS3D /echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/214-gre_field_mapping_ch1_-400mA_S88_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/212-gre_field_mapping_ch1_-400mA_S90_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/210-gre_field_mapping_ch1_+400mA_S92_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/208-gre_field_mapping_ch1_+400mA_S94_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/206-gre_field_mapping_ch2_-400mA_S96_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/204-gre_field_mapping_ch2_-400mA_S98_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/202-gre_field_mapping_ch2_+400mA_S101_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/200-gre_field_mapping_ch2_+400mA_S103_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/198-gre_field_mapping_ch3_-400mA_S105_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/196-gre_field_mapping_ch3_-400mA_S107_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/194-gre_field_mapping_ch3_+400mA_S109_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/192-gre_field_mapping_ch3_+400mA_S111_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/190-gre_field_mapping_ch4_-400mA_S113_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/188-gre_field_mapping_ch4_-400mA_S115_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/186-gre_field_mapping_ch4_+400mA_S117_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/184-gre_field_mapping_ch4_+400mA_S119_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/182-gre_field_mapping_ch5_-400mA_S121_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/180-gre_field_mapping_ch5_-400mA_S123_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/178-gre_field_mapping_ch5_+400mA_S125_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/176-gre_field_mapping_ch5_+400mA_S127_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/174-gre_field_mapping_ch6_-400mA_S129_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/172-gre_field_mapping_ch6_-400mA_S131_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/170-gre_field_mapping_ch6_+400mA_S133_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/168-gre_field_mapping_ch6_+400mA_S135_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/166-gre_field_mapping_ch7_-400mA_S137_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/164-gre_field_mapping_ch7_-400mA_S139_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/162-gre_field_mapping_ch7_+400mA_S141_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/160-gre_field_mapping_ch7_+400mA_S143_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/158-gre_field_mapping_ch8_-400mA_S145_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/156-gre_field_mapping_ch8_-400mA_S147_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/154-gre_field_mapping_ch8_+400mA_S149_DIS3D/echo_7.38/' ;
%         '/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_28p/152-gre_field_mapping_ch8_+400mA_S151_DIS3D/echo_7.38/' ; } ;
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
%
% % 1st 2 directories correspond to the baseline shim 
% Params.dataLoadDirectories{1,1,1} = tmp{1} ;
% Params.dataLoadDirectories{1,2,1} = tmp{2} ;
%
% nImgPerCurrent = 2 ; % = 1 mag image + 1 phase
%
% % dbstop in ShimOpt_Acdc at 326
% disp( ['Preparing shim calibration...' ] )        
%
% for iChannel = 1 : Params.nChannels
%     disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
%     
%     for iCurrent = 1 : Params.nCurrents 
%         Params.dataLoadDirectories{ iCurrent, 1, iChannel + 1} = tmp{ nImgPerCurrent*(Params.nCurrents*iChannel + iCurrent) -3 } ;
%         Params.dataLoadDirectories{ iCurrent, 2, iChannel + 1} = tmp{ nImgPerCurrent*(Params.nCurrents*iChannel + iCurrent) -2 } ;
%     end
% end

Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;

Params.Filtering.isFiltering  = false ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

mag = Mag.img ;

disp(['Loading magnitude images to determine region of sufficient SNR for shim field assessment' ] )        

for iChannel = 1 : Params.nChannels
    for iCurrent = 1 : Params.nCurrents
        disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
        iImg = Params.nCurrents*iChannel + iCurrent - 1 ;
        Tmp  = MaRdI( Params.dataLoadDirectories{ 1, 1, iCurrent, iChannel + 1 } ) ;
        mag(:,:,:, iImg) = Tmp.img ;
    end
end

mag = mean( mag, 4 ) ;
Params.reliabilityMask = mag/max(mag(:)) > 0.01 ; % region of reliable SNR for unwrapping
Params.reliabilityMask = shaver( Params.reliabilityMask, [3 3 1] ) ;

Params.Extension.isExtending = false ; % harmonic field extrapolation 
Params.Extension.voxelSize   = voxelSize ;
Params.Extension.expansionOrder = 1 ;
Params.Extension.radius     = 6 ;

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
