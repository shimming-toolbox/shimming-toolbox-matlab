classdef ShimOpt_IUGM_Prisma_fit < ShimOpt
%SHIMOPTUNFPRISMA - Shim Optimization for Prisma-fit @ UNF 
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
function Shim = ShimOpt_IUGM_Prisma_fit( Params, Field )
%SHIMOPTACDC - Shim Optimization

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

    Params = ShimOpt_IUGM_Prisma_fit.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_IUGM_Prisma_fit.calibratereferencemaps( Params ) ;

elseif ~isempty(Params.pathToShimReferenceMaps)
   
   [ Shim.img, Shim.Hdr ] = ShimOpt.loadshimreferencemaps( Params.pathToShimReferenceMaps ) ; 

    % DICOM Hdr.Private_0019_1013 describes the absolute table position (of the
    % reference maps).  It's used in ShimOpt.interpolatetoimggrid() to account for
    % variable table positioning, but for the scanner shims the field shift
    % doesn't depend on table position, so, (in case it's not already empty) :
    Shim.Hdr.Private_0019_1013 = [] ;

end

Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 


if (nargin == 2) && (~isempty(Field))
    
    Shim.setoriginalfield( Field ) ;

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
% DEFAULT_ISCALIBRATINGREFERENCEMAPS = false ;
%
% DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
%
% DEFAULT_PROBESPECS = [] ;


DEFAULT_ISCALIBRATINGREFERENCEMAPS = false ;
DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMaps_IUGM_Prisma_fit_20180502';
DEFAULT_PROBESPECS              = [] ;

if ~myisfield( Params, 'isCalibratingReferenceMaps' ) || isempty(Params.isCalibratingReferenceMaps)
   Params.isCalibratingReferenceMaps = DEFAULT_ISCALIBRATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   
    if Params.isCalibratingReferenceMaps
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ '~/Projects/Shimming/Static/Calibration/Data/' ...
                        'ShimReferenceMaps_IUGM_Prisma_fit_' today ] ;
    else
        Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
    end

end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_PROBESPECS ;
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

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories = cell( Params.nCurrents, 2, Params.nChannels ) ;

Params.currents = zeros( Params.nChannels, Params.nCurrents ) ; 
% will read shim current offsets (relative to baseline 'tune-up' values)
% directly from Siemens DICOM Hdr below, but they should be:
% Params.currents = [ -30 30; % A11
%                     -30 30; % B11
%                     -30 30; % A10
%                     -600 600; % A20
%                     -600 600; % A21
%                     -600 600; % B21
%                     -600 600; % A22
%                     -600 600;] ; % B22
tmp = { ...
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/149-eld_mapping_shim0_axial_fovPhase100perc_phaseOver0perc_S76_DIS3D/echo_7.38/'  ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/151-eld_mapping_shim0_axial_fovPhase100perc_phaseOver0perc_S78_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/153-gre_field_mapping_A11_minus30_S80_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/155-gre_field_mapping_A11_minus30_S82_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/157-gre_field_mapping_A11_plus30_S84_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/159-gre_field_mapping_A11_plus30_S86_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/161-gre_field_mapping_B11_minus30_S88_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/163-gre_field_mapping_B11_minus30_S90_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/165-gre_field_mapping_B11_plus30_S92_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/167-gre_field_mapping_B11_plus30_S94_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/169-gre_field_mapping_A10_minus30_S96_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/171-gre_field_mapping_A10_minus30_S98_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/173-gre_field_mapping_A10_plus30_S101_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/175-gre_field_mapping_A10_plus30_S103_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/177-gre_field_mapping_A20_minus600_S105_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/179-gre_field_mapping_A20_minus600_S107_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/181-gre_field_mapping_A20_plus600_S109_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/183-gre_field_mapping_A20_plus600_S111_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/185-gre_field_mapping_A21_minus600_S113_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/187-gre_field_mapping_A21_minus600_S115_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/189-gre_field_mapping_A21_plus600_S117_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/191-gre_field_mapping_A21_plus600_S119_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/193-gre_field_mapping_B21_minus600_S121_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/195-gre_field_mapping_B21_minus600_S123_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/197-gre_field_mapping_B21_plus600_S125_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/199-gre_field_mapping_B21_plus600_S127_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/201-gre_field_mapping_A22_minus800_S129_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/203-gre_field_mapping_A22_minus800_S131_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/209-gre_field_mapping_A22_plus600_S137_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/211-gre_field_mapping_A22_plus600_S139_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/213-gre_field_mapping_B22_minus600_S141_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/215-gre_field_mapping_B22_minus600_S143_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/217-gre_field_mapping_B22_plus600_S145_DIS3D/echo_7.38/' ;
    '~/Projects/Shimming/Static/Calibration/Data/acdc_21p/219-gre_field_mapping_B22_plus600_S147_DIS3D/echo_7.38/' ; } ;

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
       
        % for calibration of Siemens (e.g. Prisma) scanner shims only : 
        % load one of the images for each 'current' to get the shim values directly from the Siemens Hdr
        Img = MaRdI( Params.dataLoadDirectories{ iCurrent, 1, iChannel +1 }  ) ; % mag
        [f0,g0,s0] = Img.adjvalidateshim( ) ;
        shimValues = ShimOpt_IUGM_Prisma_fit.converttomultipole( [g0 ; s0] ) ; % convert to the 'multipole units' of the 3D shim card (Siemens console GUI)
        Params.currents( iChannel, iCurrent ) = shimValues( iChannel ) ; % TODO : consistent approach to units, since these aren't in amps...
    end
end

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = (Mag.img/max(Mag.img(:))) > 0.1 ; % region of reliable SNR for unwrapping


Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize = voxelSize ;
Params.Extension.radius     = 8 ;
Params.Extension.expansionOrder = 2 ;

Params.unwrapper = 'AbdulRahman_2007' ;        

end
% =========================================================================
function [ img, Hdr ] = calibratereferencemaps( Params )
%CALIBRATEREFERENCEMAPS
% 
% Wraps to ShimOpt.mapdbdi( )
% 
% [ img, Hdr ] = CALIBRATEREFERENCEMAPS( Params )
%
% Differences for ShimOpt_IUGM_Prisma_fit :
%   
%   1. After the actual mapping of the shim fields, the '0th' order frequency
%   offset is added to the array stack as img(:,:,:,1) (i.e. gradient terms
%   become img(:,:,:,2) for X, img(:,:,:,3) for Y, etc.)
%
%   2. Hdr.Private_0019_1013 is made empty 
%       This field describes the absolute table position of the reference maps.
%       It is used (i.e. for multi-coil shim arrays, for which the shim
%       position, relative to ISO, varies with the table). Having this empty
%       will indicate that values of these shim fields are fixed relative to ISO.

[ img, Hdr ] = ShimOpt.mapdbdi( Params ) ;

% insert channel (not measured) for "0th order" (Larmor transmit frequency):
tmp = zeros( size(img) + [0 0 0 1] ) ;
tmp(:,:,:,1)     = double( Hdr.MaskingImage ); % assigning the same support as the other channels
tmp(:,:,:,2:end) = img ;

img = tmp;

% DICOM .Hdr.Private_0019_1013 
%
% Describes the absolute table position of the reference maps.  This is
% used in ShimOpt.interpolatetoimggrid( ) to account for variable table
% positioning relative to ISO, but, since the scanner shims are fixed in place:
Hdr.Private_0019_1013 = [] ;

disp(['Saving shim reference maps for future use: '])
disp( Params.pathToShimReferenceMaps ) ;

save( Params.pathToShimReferenceMaps, 'img', 'Hdr' ) ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ shimValues  ] = converttomultipole( shimValues )
%CONVERTTOMULTIPOLE
% 
% shimValues = CONVERTTOMULTIPOLE( shimValues )
%
% Shim values stored in MrProt (private Siemens DICOM.Hdr) are in units of 
% DAC counts for the gradient offsets and in units of mA for the 2nd order shims.
% CONVERTTOMULTIPOLE uses the information given by the Siemens commandline tool
%   AdjValidate -shim -info
% to convert a vector of shim settings in those units into the "multipole" values
% which are used in the Siemens GUI display (i.e. Shim3d)
%
%TODO
%   Refactor and move the method to ShimCom_IUGM_Prisma_fit() 

nChannels = numel( shimValues ) ;

if nChannels == 3 
    % input shimValues are gradient offsets [units : DAC counts]
    % output shimValues units : micro-T/m]
    
    shimValues(1) = 2300*shimValues(1)/14436 ;
    shimValues(2) = 2300*shimValues(2)/14265 ;
    shimValues(3) = 2300*shimValues(3)/14045 ;

elseif nChannels == 5
    % input shimValues are for the 2nd order shims [units : mA]
    % output shimValues units : micro-T/m^2]

    shimValues(1) = 4959.01*shimValues(1)/9998 ;
    shimValues(2) = 3551.29*shimValues(2)/9998 ;
    shimValues(3) = 3503.299*shimValues(3)/9998 ;
    shimValues(4) = 3551.29*shimValues(4)/9998 ;
    shimValues(5) = 3487.302*shimValues(5)/9998 ;

elseif nChannels == 8

    shimValues(1) = 2300*shimValues(1)/14436 ;
    shimValues(2) = 2300*shimValues(2)/14265 ;
    shimValues(3) = 2300*shimValues(3)/14045 ;

    shimValues(4) = 4959.01*shimValues(4)/9998 ;
    shimValues(5) = 3551.29*shimValues(5)/9998 ;
    shimValues(6) = 3503.299*shimValues(6)/9998 ;
    shimValues(7) = 3551.29*shimValues(7)/9998 ;
    shimValues(8) = 3487.302*shimValues(8)/9998 ;

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
