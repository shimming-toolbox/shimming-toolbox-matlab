classdef ShimOpt_Rriyan < ShimOpt
%SHIMOPT_RRIYAN - Shim Optimization for RRI 24 channel array (aka spine shim)
%
% ShimOpt_RRIYAN is a ShimOpt subclass 
%     
% =========================================================================
% Updated::20190528::ryan.topfer@polymtl.ca
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
function Shim = ShimOpt_Rriyan( Params, Field )
%SHIMOPT_RRYAN - Shim Optimization

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;
Shim.System.Specs    = ShimSpecs_Rriyan();
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_Rriyan.assigndefaultparameters( Params )

if Params.isCalibratingReferenceMaps
    
    Params = ShimOpt_Rriyan.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_Rriyan.calibratereferencemaps( Params ) ;

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
%OPTIMIZESTATICSHIMCURRENTS 


DEFAULT_ISOPTIMIZINGAUX             = true ;

if nargin < 2 
    Params.dummy = [];
end

if ~myisfield(Params, 'isOptimizingAux') || isempty( Params.isOptimizingAux )
    Params.isOptimizingAux = DEFAULT_ISOPTIMIZINGAUX ;
end

if nargin < 3
    Corrections = optimizeshimcurrents@ShimOpt( Shim, Params, @checknonlinearconstraints ) ;
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

DEFAULT_ISCALIBRATINGREFERENCEMAPS = false ;
DEFAULT_PATHTOSHIMREFERENCEMAPS    = [ shimbindir() 'ShimReferenceMaps_Rriyan' ] ;
DEFAULT_PROBESPECS                 = [] ;

DEFAULT_INSTITUTIONNAME = 'IUGM' ;
DEFAULT_STATIONNAME     = 'MRC35049' ;

if ~myisfield( Params, 'isCalibratingReferenceMaps' ) || isempty(Params.isCalibratingReferenceMaps)
   Params.isCalibratingReferenceMaps = DEFAULT_ISCALIBRATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   
    if Params.isCalibratingReferenceMaps
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ shimbindir() 'ShimReferenceMaps_Rriyan'  ] ;
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

Params.nChannels  = 24 ;
Params.nCurrents  = 2 ;
Params.nEchoes    = 1 ; % nEchoes = 1 for phase *difference* images

Params.currents = repmat( [0.5 -0.5], [Params.nChannels 1] ); 

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

tmp = { ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/111-gre_field_mapping_ch24_+0.5A_S110_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/112-gre_field_mapping_ch24_+0.5A_S109_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/113-gre_field_mapping_ch24_-0.5A_S108_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/114-gre_field_mapping_ch24_-0.5A_S107_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/115-gre_field_mapping_ch23_+0.5A_S106_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/116-gre_field_mapping_ch23_+0.5A_S105_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/117-gre_field_mapping_ch23_-0.5A_S104_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/118-gre_field_mapping_ch23_-0.5A_S103_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/119-gre_field_mapping_ch22_+0.5A_S102_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/120-gre_field_mapping_ch22_+0.5A_S101_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/121-gre_field_mapping_ch22_-0.5A_S100_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/122-gre_field_mapping_ch22_-0.5A_S98_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/123-gre_field_mapping_ch21_+0.5A_S97_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/124-gre_field_mapping_ch21_+0.5A_S96_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/125-gre_field_mapping_ch21_-0.5A_S95_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/126-gre_field_mapping_ch21_-0.5A_S94_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/127-gre_field_mapping_ch20_+0.5A_S93_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/128-gre_field_mapping_ch20_+0.5A_S92_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/129-gre_field_mapping_ch20_-0.5A_S91_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/130-gre_field_mapping_ch20_-0.5A_S90_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/131-gre_field_mapping_ch19_+0.5A_S89_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/132-gre_field_mapping_ch19_+0.5A_S88_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/133-gre_field_mapping_ch19_-0.5A_S87_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/134-gre_field_mapping_ch19_-0.5A_S86_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/135-gre_field_mapping_ch18_+0.5A_S85_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/136-gre_field_mapping_ch18_+0.5A_S84_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/137-gre_field_mapping_ch18_-0.5A_S83_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/138-gre_field_mapping_ch18_-0.5A_S82_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/139-gre_field_mapping_ch17_+0.5A_S81_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/140-gre_field_mapping_ch17_+0.5A_S80_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/141-gre_field_mapping_ch17_-0.5A_S79_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/142-gre_field_mapping_ch17_-0.5A_S78_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/143-gre_field_mapping_ch16_+0.5A_S77_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/144-gre_field_mapping_ch16_+0.5A_S76_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/145-gre_field_mapping_ch16_-0.5A_S75_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/146-gre_field_mapping_ch16_-0.5A_S74_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/147-gre_field_mapping_ch15_+0.5A_S73_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/148-gre_field_mapping_ch15_+0.5A_S72_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/149-gre_field_mapping_ch15_-0.5A_S71_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/150-gre_field_mapping_ch15_-0.5A_S70_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/151-gre_field_mapping_ch14_+0.5A_S69_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/152-gre_field_mapping_ch14_+0.5A_S68_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/153-gre_field_mapping_ch14_-0.5A_S67_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/154-gre_field_mapping_ch14_-0.5A_S66_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/155-gre_field_mapping_ch13_+0.5A_S65_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/156-gre_field_mapping_ch13_+0.5A_S64_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/157-gre_field_mapping_ch13_-0.5A_S63_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/158-gre_field_mapping_ch13_-0.5A_S62_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/159-gre_field_mapping_ch12_+0.5A_S61_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/160-gre_field_mapping_ch12_+0.5A_S60_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/161-gre_field_mapping_ch12_-0.5A_S59_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/162-gre_field_mapping_ch12_-0.5A_S58_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/163-gre_field_mapping_ch11_+0.5A_S57_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/164-gre_field_mapping_ch11_+0.5A_S56_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/165-gre_field_mapping_ch11_-0.5A_S55_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/166-gre_field_mapping_ch11_-0.5A_S54_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/167-gre_field_mapping_ch10_+0.5A_S53_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/168-gre_field_mapping_ch10_+0.5A_S52_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/169-gre_field_mapping_ch10_-0.5A_S51_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/170-gre_field_mapping_ch10_-0.5A_S50_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/171-gre_field_mapping_ch09_+0.5A_S49_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/172-gre_field_mapping_ch09_+0.5A_S48_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/173-gre_field_mapping_ch09_-0.5A_S47_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/174-gre_field_mapping_ch09_-0.5A_S46_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/175-gre_field_mapping_ch08_+0.5A_S45_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/176-gre_field_mapping_ch08_+0.5A_S44_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/177-gre_field_mapping_ch08_-0.5A_S43_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/178-gre_field_mapping_ch08_-0.5A_S42_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/179-gre_field_mapping_ch07_+0.5A_S41_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/180-gre_field_mapping_ch07_+0.5A_S40_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/181-gre_field_mapping_ch07_-0.5A_S39_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/182-gre_field_mapping_ch07_-0.5A_S38_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/185-gre_field_mapping_ch06_+0.5A_S35_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/186-gre_field_mapping_ch06_+0.5A_S34_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/187-gre_field_mapping_ch06_-0.5A_S33_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/188-gre_field_mapping_ch06_-0.5A_S32_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/189-gre_field_mapping_ch05_+0.5A_S31_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/190-gre_field_mapping_ch05_+0.5A_S30_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/191-gre_field_mapping_ch05_-0.5A_S29_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/192-gre_field_mapping_ch05_-0.5A_S28_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/193-gre_field_mapping_ch04_+0.5A_S27_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/194-gre_field_mapping_ch04_+0.5A_S26_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/195-gre_field_mapping_ch04_-0.5A_S25_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/196-gre_field_mapping_ch04_-0.5A_S24_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/197-gre_field_mapping_ch03_+0.5A_S23_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/198-gre_field_mapping_ch03_+0.5A_S22_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/199-gre_field_mapping_ch03_-0.5A_S21_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/200-gre_field_mapping_ch03_-0.5A_S20_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/201-gre_field_mapping_ch02_+0.5A_S19_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/202-gre_field_mapping_ch02_+0.5A_S18_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/203-gre_field_mapping_ch02_-0.5A_S17_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/204-gre_field_mapping_ch02_-0.5A_S16_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/205-gre_field_mapping_ch01_+0.5A_S15_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/206-gre_field_mapping_ch01_+0.5A_S14_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/207-gre_field_mapping_ch01_-0.5A_S13_DIS3D/echo_7.38/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/208-gre_field_mapping_ch01_-0.5A_S12_DIS3D/echo_7.38/' ; } ;

tmp = flipud( tmp ) ;

Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

for iChannel = 1 : Params.nChannels
    for iCurrent = 1 : Params.nCurrents
        for imgType = 1 : 2 % 1=mag, 2=phase
        
            iDir = (iChannel-1)*(Params.nCurrents*2) + 1 ;

            Params.dataLoadDirectories{ 1, imgType, iCurrent, iChannel } = [ tmp{ iDir + (imgType -1) + 2*(iCurrent-1) } ] ;

        end
    end 
end

Params.unwrapper = 'AbdulRahman_2007' ;        

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

