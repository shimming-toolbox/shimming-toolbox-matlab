classdef ShimOpt_IUGM_Prisma_fit < ShimOpt
%SHIMOPT_IUGM_PRISMA_FIT - Shim Optimization for Prisma-fit @ UNF 
%     
% ShimOpt_IUGM_Prisma_fit is a ShimOpt subclass. See ShimOpt documentation for
% usage.
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

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
Shim.System.Specs    = ShimSpecs_IUGM_Prisma_fit();
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt_IUGM_Prisma_fit.assigndefaultparameters( Params ) ;

if Params.isCalibratingReferenceMaps

    Params = ShimOpt_IUGM_Prisma_fit.declarecalibrationparameters( Params ) ;
    [ Shim.img, Shim.Hdr ] = ShimOpt_IUGM_Prisma_fit.calibratereferencemaps( Params ) ;

elseif ~isempty(Params.pathToShimReferenceMaps)
   
   [ Shim.img, Shim.Hdr ] = ShimOpt.loadshimreferencemaps( Params.pathToShimReferenceMaps ) ; 
    
    Shim.Ref.img = Shim.img ;
    Shim.Ref.Hdr = Shim.Hdr ;

end

Params.TrackerSpecs.state = 'inert' ;
Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 

if (nargin == 2) && (~isempty(Field))
    
    Shim.setoriginalfield( Field ) ;

end

end
% =========================================================================
function [] = interpolatetoimggrid( Shim, Field )
%INTERPOLATETOIMGGRID 
%
% [] = INTERPOLATETOIMGGRID( Shim, Field )
%
% Interpolates Shim.img (reference maps) to the grid (voxel positions) of
% MaRdI-type Img
% 
% i.e.
%
%   [X,Y,Z] = Field.getvoxelpositions ;
%   Shim.resliceimg( X, Y, Z ) ;
%
% NOTE
%
%   The patient coordinate system is defined by the initial (laser) placement
%   of the subject. After the 1st localizer (for which the Z=0 position will
%   correspond to isocenter), it is likely that the operator will choose a
%   particular FOV for the following scans, thereby repositioning the table by
%   a certain amount ( Field.Hdr.Img.ImaRelTablePosition ).  i.e. Isocenter has
%   been moved from Z=0 to Z = Field.Hdr.Img.ImaRelTablePosition.
% 
%   For our multi-coil shim arrays, the shim moves along with the table (as
%   does the patient coordinate system), so a shim field shift at initial
%   location r' = (x',y',z') will continue to be exactly that.
%
%   The scanner shims, on the other hand, are fixed relative to isocenter. So a
%   shim field shift induced at initial table position r', will now instead be
%   induced at r' + Field.Hdr.Img.ImaRelTablePosition.

[X, Y, Z]    = Field.getvoxelpositions ;
[X0, Y0, Z0] = Shim.getvoxelpositions ;

dR = Field.isocenter() ; 
assert( dR(1) == 0, 'Table shifted in L/R direction?' ) ;
assert( dR(2) == 0, 'Table shifted in A/P direction?' ) ;

if ( dR(3) ~= 0 ) % field positions originally at Z0 have been shifted
    % NOTE
    %   tablePosition is increasingly negative the more it is into the scanner.
    %   the opposite is true for the z-coordinate of a voxel in the dicom
    %   reference system.
    warning('Correcting for table shift with respect to shim reference images')
    Z0 = Z0 + dR(3) ;
    Shim.Hdr.ImagePositionPatient(3) = Shim.Hdr.ImagePositionPatient(3) + dR(3) ;   
end

% -------
% check if voxel positions already happen to coincide. if they do, don't interpolate (time consuming).
if ~MaRdI.compareimggrids( X, Y, Z, X0, Y0, Z0 )
    Shim.resliceimg( X, Y, Z ) ;
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
%   .maxCurrentPerChannel
%       [default: determined by class ShimSpecs.Amp.maxCurrentPerChannel]
 
if nargin < 2 
    Params.dummy = [];
end

Corrections = optimizeshimcurrents@ShimOpt( Shim, Params  ) ;

function [C, Ceq] = checknonlinearconstraints( corrections )
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
    C = abs( corrections ) - Params.maxCurrentPerChannel ;

end

end
% =========================================================================
function [] = setoriginalfield( Shim, Field )
%SETORIGINALFIELD 
%
% [] = SETORIGINALFIELD( Shim, Field )
%
% Sets Shim.Field
%
% Field is a FieldEval type object with .img in Hz

Shim.Field = Field.copy() ;

Shim.interpolatetoimggrid( Shim.Field ) ;
Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

% get the original shim offsets
[f0, g0, s0]  = Shim.Field.adjvalidateshim() ;
Shim.System.currents            =  [ ShimOpt_IUGM_Prisma_fit.converttomultipole( [g0 ; s0] ) ] ; 
Shim.System.Tx.imagingFrequency = f0 ;

% if ~isempty( Shim.Aux ) && ~isempty( Shim.Aux.Shim ) 
%     Shim.Aux.Shim.Field = Shim.Field ;
%     Shim.Aux.Shim.interpolatetoimggrid( Shim.Field ) ;
%     Shim.Aux.Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;
% end

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
DEFAULT_PATHTOSHIMREFERENCEMAPS    = [ shimbindir() 'ShimReferenceMaps_IUGM_Prisma_fit' ] ;
DEFAULT_PROBESPECS                 = [] ;

if ~myisfield( Params, 'isCalibratingReferenceMaps' ) || isempty(Params.isCalibratingReferenceMaps)
   Params.isCalibratingReferenceMaps = DEFAULT_ISCALIBRATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   
    if Params.isCalibratingReferenceMaps
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ shimbindir() 'ShimReferenceMaps_IUGM_Prisma_fit' ] ;
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
Params.nEchoes    = 5 ; % nEchoes = 1 for phase *difference* images

% 2 columns: [ MAG | PHASE ] ;
Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

tmp = { ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/03-fm3d_shimX_chXX_XXmA_test/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/04-fm3d_shimX_chXX_XXmA_test/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/05-fm3d_A11_plus30/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/06-fm3d_A11_plus30/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/07-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/08-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/09-fm3d_B11_plus30/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/10-fm3d_B11_plus30/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/11-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/12-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/13-fm3d_A01_plus30/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/14-fm3d_A01_plus30/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/15-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/16-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/17-fm3d_A20_plus600/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/18-fm3d_A20_plus600/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/19-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/20-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/21-fm3d_A21_plus600/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/22-fm3d_A21_plus600/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/23-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/24-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/25-fm3d_B21_plus600/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/26-fm3d_B21_plus600/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/27-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/28-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/29-fm3d_A22_plus1000/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/30-fm3d_A22_plus1000/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/31-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/32-fm3d_baseline/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/33-fm3d_B22_plus1000/' ;
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/34-fm3d_B22_plus1000/' ; } ;

Params.dataLoadDirectories = cell( Params.nEchoes, 2, Params.nCurrents, Params.nChannels ) ;

for iChannel = 1 : Params.nChannels
    for iCurrent = 1 : Params.nCurrents
        for imgType = 1 : 2 % 1=mag, 2=phase

        iDir = (iChannel-1)*(Params.nCurrents*2) + 1 ;

        dicomSubDirs  = dir( [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo*/'] ) ;
        nDicomSubDirs = length( dicomSubDirs ) ;
        assert( nDicomSubDirs == 5 )

        Params.dataLoadDirectories{ 1, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_1.68/' ] ;
        Params.dataLoadDirectories{ 2, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_5.68/' ] ;
        Params.dataLoadDirectories{ 3, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_10.68/' ] ;
        Params.dataLoadDirectories{ 4, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_15.68/' ] ;
        Params.dataLoadDirectories{ 5, imgType, iCurrent, iChannel } = [ tmp{iDir + (imgType -1) + 2*(iCurrent-1) } 'echo_20.68/' ] ;

        end
    end 
end

disp( ['Preparing shim calibration...' ] )        

% Read shim current offsets (relative to baseline 'tune-up' values) directly from Siemens DICOM Hdr
Params.currents = zeros( Params.nChannels, Params.nCurrents ) ; 

for iChannel = 1 : Params.nChannels
    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        

    for iCurrent = 1 : Params.nCurrents 
        % for calibration of Siemens (e.g. Prisma) scanner shims only : 
        % load one of the images for each 'current' to get the shim values directly from the Siemens Hdr
        Img = MaRdI( Params.dataLoadDirectories{ 1, 1, iCurrent, iChannel }  ) ; % mag
        [f0,g0,s0] = Img.adjvalidateshim( ) ;
        shimValues = ShimOpt_IUGM_Prisma_fit.converttomultipole( [g0 ; s0] ) ; % convert to the 'multipole units' of the 3D shim card (Siemens console GUI)
        Params.currents( iChannel, iCurrent ) = shimValues( iChannel ) ; % TODO : consistent approach to units, since these aren't in amps...
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
