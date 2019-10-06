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
function Shim = ShimOpt_IUGM_Prisma_fit( varargin )
%SHIMOPT - Shim Optimization

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;
Shim.System.Specs    = ShimSpecs_IUGM_Prisma_fit();
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Params = ShimOpt_IUGM_Prisma_fit.assigndefaultparameters( Params ) ;

if Params.isCalibratingReferenceMaps

    [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt_IUGM_Prisma_fit.calibratereferencemaps( Params ) ;

elseif ~isempty(Params.pathToShimReferenceMaps)
   
   [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt.loadshimreferencemaps( Params.pathToShimReferenceMaps ) ; 
    
    Shim.Ref.img = Shim.img ;
    Shim.Ref.Hdr = Shim.Hdr ;
end

if ~isempty( Field ) 
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

DEFAULTS.isCalibratingReferenceMaps = false ;
DEFAULTS.pathToShimReferenceMaps    = [ shimbindir() 'ShimReferenceMaps_IUGM_Prisma_fit' ] ;

if ~myisfieldfilled( Params, 'isCalibratingReferenceMaps' ) 
   Params.isCalibratingReferenceMaps = DEFAULTS.isCalibratingReferenceMaps ;
end

if ~myisfieldfilled( Params, 'pathToShimReferenceMaps' ) 
    if Params.isCalibratingReferenceMaps
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ shimbindir() 'ShimReferenceMaps_IUGM_Prisma_fit' ] ;
    else
        Params.pathToShimReferenceMaps = DEFAULTS.pathToShimReferenceMaps ;
    end
end

end
% =========================================================================
function [ img, Hdr, Interpolant ] = calibratereferencemaps( Params )
%CALIBRATEREFERENCEMAPS

% Wraps to .mapdbdi( ) and writes output to disk
% 
% [ img, Hdr ] = CALIBRATEREFERENCEMAPS( Params )

Params = ShimOpt_IUGM_Prisma_fit.declarecalibrationparameters( Params ) ;

[ img, Hdr, Interpolant ] = ShimOpt_IUGM_Prisma_fit.mapdbdi( Params ) ;

disp(['Saving shim reference maps for future use: '])
disp( Params.pathToShimReferenceMaps ) ;

save( Params.pathToShimReferenceMaps, 'img', 'Hdr', 'Interpolant' ) ;

end
% =========================================================================
function [ dataLoadDirectories ] = getdataloaddirectories( seriesName )
%GETDATALOADDIRECTORIES  
%
% [ dataLoadDirectories ] = GETDATALOADDIRECTORIES( seriesName )
%
% seriesName can be 'acdc_39p' or 'acdc_82p'

switch seriesName

    case 'acdc_39p'
        % 2 columns: [ MAG | PHASE ] ;
        dataLoadDirectories    = cell( 2, 2, 8 ) ;

        %% ------
        % channel 1 : "A11" (x-gradient)

        % current 1:
        dataLoadDirectories{1,1,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/03-fm3d_shimX_chXX_XXmA_test/' ;
        dataLoadDirectories{1,2,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/04-fm3d_shimX_chXX_XXmA_test/' ;
        dataLoadDirectories{2,1,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/05-fm3d_A11_plus30/' ;
        dataLoadDirectories{2,2,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/06-fm3d_A11_plus30/' ;

        %% ------
        % channel 2 : "B11" (y-gradient)
        dataLoadDirectories{1,1,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/07-fm3d_baseline/' ;
        dataLoadDirectories{1,2,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/08-fm3d_baseline/' ;
        dataLoadDirectories{2,1,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/09-fm3d_B11_plus30/' ;
        dataLoadDirectories{2,2,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/10-fm3d_B11_plus30/' ;

        %% ------
        % channel 3 : "A01" (z-gradient)
        dataLoadDirectories{1,1,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/11-fm3d_baseline/' ;
        dataLoadDirectories{1,2,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/12-fm3d_baseline/' ;
        dataLoadDirectories{2,1,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/13-fm3d_A01_plus30/' ;
        dataLoadDirectories{2,2,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/14-fm3d_A01_plus30/' ;

        %% ------
        % channel 4 : "A20"
        dataLoadDirectories{1,1,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/15-fm3d_baseline/' ;
        dataLoadDirectories{1,2,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/16-fm3d_baseline/' ;
        dataLoadDirectories{2,1,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/17-fm3d_A20_plus600/' ;
        dataLoadDirectories{2,2,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/18-fm3d_A20_plus600/' ;
            
        %% ------
        % channel 5 : "A21"
        dataLoadDirectories{1,1,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/19-fm3d_baseline/' ;
        dataLoadDirectories{1,2,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/20-fm3d_baseline/' ;
        dataLoadDirectories{2,1,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/21-fm3d_A21_plus600/' ;
        dataLoadDirectories{2,2,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/22-fm3d_A21_plus600/' ;

        %% ------
        % channel 6 : "B21"
        dataLoadDirectories{1,1,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/23-fm3d_baseline/' ;
        dataLoadDirectories{1,2,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/24-fm3d_baseline/' ;
        dataLoadDirectories{2,1,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/25-fm3d_B21_plus600/' ;
        dataLoadDirectories{2,2,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/26-fm3d_B21_plus600/' ;

        %% ------
        % channel 7 : "A22"
        dataLoadDirectories{1,1,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/27-fm3d_baseline/' ;
        dataLoadDirectories{1,2,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/28-fm3d_baseline/' ;
        dataLoadDirectories{2,1,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/29-fm3d_A22_plus1000/' ;
        dataLoadDirectories{2,2,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/30-fm3d_A22_plus1000/' ;

        %% ------
        % channel 8 : "B22"
        dataLoadDirectories{1,1,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/27-fm3d_baseline/' ;
        dataLoadDirectories{1,2,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/28-fm3d_baseline/' ;
        dataLoadDirectories{2,1,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/33-fm3d_B22_plus1000/' ;
        dataLoadDirectories{2,2,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_39p/34-fm3d_B22_plus1000/' ;
    
    case 'acdc_82p'
        % 2 columns: [ MAG | PHASE ] ;
        dataLoadDirectories = cell( 2, 2, 8 ) ;

        dataLoadDirectories{1,1,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/05-gre_fm_A11_+50_2mmInplane/' ;
        dataLoadDirectories{1,2,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/06-gre_fm_A11_+50_2mmInplane/' ;
        dataLoadDirectories{2,1,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/07-gre_fm_A11_-50/' ;
        dataLoadDirectories{2,2,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/08-gre_fm_A11_-50/' ;

        dataLoadDirectories{1,1,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/09-gre_fm_B11_+50/' ;
        dataLoadDirectories{1,2,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/10-gre_fm_B11_+50/' ;
        dataLoadDirectories{2,1,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/11-gre_fm_B11_-50/' ;
        dataLoadDirectories{2,2,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/12-gre_fm_B11_-50/' ;

        dataLoadDirectories{1,1,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/13-gre_fm_A10_+50/' ;
        dataLoadDirectories{1,2,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/14-gre_fm_A10_+50/' ;
        dataLoadDirectories{2,1,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/15-gre_fm_A10_-50/' ;
        dataLoadDirectories{2,2,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/16-gre_fm_A10_-50/' ;

        dataLoadDirectories{1,1,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/17-gre_fm_A20_+600/' ;
        dataLoadDirectories{1,2,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/18-gre_fm_A20_+600/' ;
        dataLoadDirectories{2,1,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/19-gre_fm_A20_-600/' ;
        dataLoadDirectories{2,2,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/20-gre_fm_A20_-600/' ;

        dataLoadDirectories{1,1,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/21-gre_fm_A21_+600/' ;
        dataLoadDirectories{1,2,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/22-gre_fm_A21_+600/' ;
        dataLoadDirectories{2,1,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/23-gre_fm_A21_-600/' ;
        dataLoadDirectories{2,2,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/24-gre_fm_A21_-600/' ;

        dataLoadDirectories{1,1,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/25-gre_fm_B21_+600/' ;
        dataLoadDirectories{1,2,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/26-gre_fm_B21_+600/' ;
        dataLoadDirectories{2,1,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/27-gre_fm_B21_-600/' ;
        dataLoadDirectories{2,2,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/28-gre_fm_B21_-600/' ;

        dataLoadDirectories{1,1,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/29-gre_fm_A22_+600/' ;
        dataLoadDirectories{1,2,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/30-gre_fm_A22_+600/' ;
        dataLoadDirectories{2,1,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/31-gre_fm_A22_-600/' ;
        dataLoadDirectories{2,2,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/32-gre_fm_A22_-600/' ;

        dataLoadDirectories{1,1,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/33-gre_fm_B22_+600/' ;
        dataLoadDirectories{1,2,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/34-gre_fm_B22_+600/' ;
        dataLoadDirectories{2,1,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/35-gre_fm_B22_-600/' ;
        dataLoadDirectories{2,2,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_82p/36-gre_fm_B22_-600/' ;

    otherwise
        error('series not implemented') ;
end

end
% =========================================================================
function Params = declarecalibrationparameters( Params )
%DECLARECALIBRATIONPARAMETERS
% 
% Initializes parameters for shim reference map construction (aka shim calibration)

%% -----
DEFAULTS.unwrapper = 'AbdulRahman_2007' ;
DEFAULTS.threshold = 0.05 ;

Params.dummy = [] ;
Params = assignifempty( Params, DEFAULTS ) ;
%

Specs             = ShimSpecs_IUGM_Prisma_fit() ;
Params.nChannels  = Specs.Amp.nActiveChannels ;

Params.nCurrents  = 2 ;

% for multi-coil systems, entries of Params.currents need to be manually
% defined.  for the host system, the 'currents' (shim coefficients in
% multi-pole units) can be retrieved directly from the dicoms 
% See ShimOpt_IUGM_Prisma_fit.mapdbdi( )
Params.currents   = zeros( Params.nChannels, Params.nCurrents ) ; 

% image data directories
if ~myisfieldfilled( Params, 'dataLoadDirectories' )
    Params.dataLoadDirectories = ShimOpt_IUGM_Prisma_fit.getdataloaddirectories('acdc_82p') ;
end

end
% =========================================================================
function [ img, Hdr, Interpolant ] = mapdbdi( Params )
%MAPDBDI    map dB/dI : field shift [Hz] per unit-shim
% 
% [ img, Hdr, Interpolant ] = MAPDBDI( Params ) 

DEFAULTS.unwrapper = 'AbdulRahman_2007' ; 
DEFAULTS.threshold = 0.05 ; 

if nargin == 0 || isempty( Params )
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

img   = [] ;
Hdr   = [] ;
dBdI  = [] ;

nChannels = size( Params.currents, 1 ) ;
nCurrents = size( Params.currents, 2 ) ;

assert( nCurrents == 2 )

assert( ( size( Params.dataLoadDirectories, 1 ) == nCurrents ) ...
    && ( size( Params.dataLoadDirectories, 3 ) == nChannels ), ...
    'mapdbdi requires 2 calibration currents be defined for each channel, with corresponding image folders defined for each' )

Mag    = cell( nCurrents, 1, 8 ) ;
Phase  = cell( nCurrents, 1, 8 ) ;

disp( ['Preparing shim calibration...' ] )        
for iChannel = 1 : nChannels 
    disp(['Channel ' num2str(iChannel) ' of ' num2str(nChannels) ] )        

    for iCurrent = 1  : nCurrents
        Mag{ iCurrent, 1, iChannel }   = MaRdI( Params.dataLoadDirectories{ iCurrent, 1, iChannel } ) ;
        Phase{ iCurrent, 1, iChannel } = MaRdI( Params.dataLoadDirectories{ iCurrent, 2, iChannel } ) ;

        %% -----
        % For calibration of Siemens (e.g. Prisma) scanner shims only :
        % Read shim current offsets (relative to baseline values) directly
        % from Siemens DICOM Hdr. 
        [f0,g0,s0] = Mag{ iCurrent, 1, iChannel}.adjvalidateshim( ) ;
        % convert to the 'multipole units' (micro-T/m and micro-T/m^2) of the 3D shim card (Siemens console GUI)
        shimValues = ShimSpecs_IUGM_Prisma_fit.converttomultipole( [g0 ; s0] ) ; 
        Params.currents( iChannel, iCurrent ) = shimValues( iChannel ) ; 
    end
end

%% -----
% assign binary mask for unwrapping based on net magnitude across acquisition series
avgMag = zeros( Mag{1,1,1}.getgridsize() ) ;

for iChannel = 1 : nChannels 
    for iCurrent = 1  : nCurrents
        % add all magnitude images, (averaging over the 2 echoes)    
        avgMag = avgMag + ( Mag{ iCurrent, 1, iChannel }.img(:,:,:,1) + Mag{iCurrent, 1, iChannel}.img(:,:,:,2 ) )/2 ;
    end
end
% normalize
avgMag = avgMag/max(avgMag(:)) ;
mask   = avgMag > Params.threshold ;

%% -----
% compute dB/dI 
Fields = cell( nCurrents, 1, 8 ) ;

for iChannel = 1 : nChannels 
    disp(['Computing field maps: Channel ' num2str(iChannel) ' of ' num2str(nChannels) ] )        
    for iCurrent = 1  : nCurrents
        disp(['Current ' num2str(iCurrent) ' of ' num2str(nCurrents) ] )        
        Phase{iCurrent,1,iChannel}.setmaskingimage( mask ) ;

        Fields{ iCurrent, 1, iChannel } = FieldEval.mapfield( Mag{iCurrent,1,iChannel}, ...
                                                     Phase{iCurrent,1,iChannel}, Params  ) ; 
    end
    
    dB = Fields{2,1,iChannel}.img - Fields{1,1,iChannel}.img ;
    dI = Params.currents(iChannel,2) - Params.currents(iChannel, 1) ;
    dBdI(:,:,:,iChannel) = dB/dI ;
end

%% -----
% for rapid interpolation of dB/dI to a given image, form the interpolant here + save for future use: 
Interpolant = Mag{1,1,1}.getinterpolant() ;

%% -----
% Due to the arbitrariness of the excitation/demodulation frequency, the dB/dI
% maps at this stage are only accurate to within a constant. (For example, if Gx =
% +1Hz/mm, and the imaging frequency is set to the local frequency at x=1mm,
% then the resulting field map will possess zero shift there, rather than at
% isocentre x=0). To correct for this, subtract the measured dB/dI value at
% isocentre:
disp( 'Correcting for demodulation: Recentering measured dB/dI to be 0Hz at isocentre...' )
for iChannel = 1 : nChannels 
    disp(['Channel ' num2str(iChannel) ' of ' num2str(nChannels) ] )        
    dBdI0 = dBdI(:,:,:,iChannel) ;
    Interpolant.Values = dBdI0(:) ;
    dBdI(:,:,:,iChannel) = dBdI0 - Interpolant([0 0 0]) ;
end

img = dBdI ;
img(~mask) = 0 ;
Hdr = Fields{1,1,1}.Hdr ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ shimValues ] = converttomultipole( shimValues )
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
