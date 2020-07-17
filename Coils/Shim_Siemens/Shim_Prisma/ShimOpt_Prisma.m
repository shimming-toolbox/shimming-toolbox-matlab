classdef ShimOpt_Prisma < ShimOpt
%SHIMOPT_Prisma - Shim Optimization for Siemens Prisma scanners 
%     
% ShimOpt_Prisma is a ShimOpt subclass. 
% See ShimOpt documentation for usage.
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt_Prisma( varargin )
%SHIMOPT - Shim Optimization

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;
Shim.System.Specs    = ShimSpecs_Prisma();
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Params = ShimOpt_Prisma.assigndefaultparameters( Params, Shim.System.Specs ) ;

switch Params.shimReferenceMaps
    case 'calibrate'
        [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt_Prisma.calibratereferencemaps( Params ) ;
        Shim.Ref.source = 'data' ;
    case 'model'  
        Shim.Ref.source = 'model' ;
    otherwise 
       [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt.loadshimreferencemaps( Params.shimReferenceMaps ) ; 
        
        Shim.Ref.img = Shim.img ;
        Shim.Ref.Hdr = Shim.Hdr ;
        Shim.Ref.source = 'data' ;
end

if ~isempty( Field ) 
    Shim.setoriginalfield( Field ) ;
end

end
% =========================================================================
function [] = comparemaps( Shim, iSlice )
%COMPAREMAPS    Compare empirical to modeled reference maps

assert( strcmp( Shim.Ref.source, 'data' ), ...
    'Input ShimOpt object should be instantiated with empirical reference maps. See HELP ShimOpt_Prisma' )

if ( nargin == 1 ) || isempty( iSlice )
    iSlice = 63 ;
end

%% -----
Params.shimReferenceMaps = 'model' ;

%% how to get field object?? 
% how to create correct ShimOpt_Prisma object??
% IdealShims = ShimOpt_Prisma( Params, Fields{1,1,1} ) ;
% IdealShims = ShimOpt_IUGM_Prisma_fit( Params, Fields{1,1,1} ) ;

[X,Y,Z] = IdealShims.getvoxelpositions() ;

R = sqrt( X.^2 + Y.^2 + Z.^2 ) ;

deviation = repmat( Shim.getshimsupport(), [1 1 1 8] ) .* abs( IdealShims.img - Shim.img ) ;

%% -----
% 1st order shims:

close all
figure

for iCh = 1:3
    switch iCh
        case 1
            iRow = 0;
        case 2
            iRow = 4;
        case 3 
            iRow = 8 ;
    end

    subplot(3,4,iRow+1) 
    contour( X(:,:,iSlice),Y(:,:,iSlice), R(:,:,iSlice), 'ShowText', 'on' ) ;
    title('Distance to isocentre (mm)') ;
    ylabel(IdealShims.System.Specs.Id.channelNames{iCh})

    subplot(3,4,iRow+2) 
    imagesc( IdealShims.img(:,:,iSlice,iCh) ) ;
    title('Nominal/simulated') ;
    caxis([-1 1] ) ;
    colorbar 
    set(gca,'XTick',[])
    set(gca,'YTick',[])

    subplot(3,4,iRow+3) 
    imagesc( Shim.img(:,:,iSlice,iCh) ) ;
    title('Empirical') ;
    caxis([-1 1] ) ;
    colorbar 
    set(gca,'XTick',[])
    set(gca,'YTick',[])

    subplot(3,4,iRow+4) 
    imagesc( deviation(:,:,iSlice,iCh) ) ;
    title('Abs. deviation') ;
    caxis([0 1] ) ;
    colorbar ;
    set(gca,'XTick',[])
    set(gca,'YTick',[])

end

%% -----
% 2nd order shims:

figure

for iCh = 4:6
    switch iCh
        case 4
            iRow = 0;
        case 5
            iRow = 4;
        case 6 
            iRow = 8 ;
    end

    subplot(3,4,iRow+1) 
    contour( X(:,:,63),Y(:,:,63), R(:,:,63), 'ShowText', 'on' ) ;
    title('Distance to isocentre (mm)') ;
    ylabel(IdealShims.System.Specs.Id.channelNames{iCh})

    subplot(3,4,iRow+2) 
    imagesc( IdealShims.img(:,:,iSlice,iCh) ) ;
    title('Nominal/simulated') ;
    caxis([-0.3 0.3] ) ;
    colorbar 
    set(gca,'XTick',[])
    set(gca,'YTick',[])

    subplot(3,4,iRow+3) 
    imagesc( Shim.img(:,:,iSlice,iCh) ) ;
    title('Empirical') ;
    caxis([-0.3 0.3] ) ;
    colorbar 
    set(gca,'XTick',[])
    set(gca,'YTick',[])

    subplot(3,4,iRow+4) 
    imagesc( deviation(:,:,iSlice,iCh) ) ;
    title('Abs. deviation') ;
    caxis([0 0.1] ) ;
    colorbar ;
    set(gca,'XTick',[])
    set(gca,'YTick',[])

end

%% -----
% 2nd order shims (cont)

figure

for iCh = 7:8 
    switch iCh
        case 7
            iRow = 0;
        case 8
            iRow = 4;
    end

    subplot(2,4,iRow+1) 
    contour( X(:,:,63),Y(:,:,63), R(:,:,63), 'ShowText', 'on' ) ;
    title('Distance to isocentre (mm)') ;
    ylabel(IdealShims.System.Specs.Id.channelNames{iCh})

    subplot(2,4,iRow+2) 
    imagesc( IdealShims.img(:,:,iSlice,iCh) ) ;
    title('Nominal/simulated') ;
    caxis([-0.3 0.3] ) ;
    colorbar 
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    grid on;

    subplot(2,4,iRow+3) 
    imagesc( Shim.img(:,:,iSlice,iCh) ) ;
    title('Empirical') ;
    caxis([-0.3 0.3] ) ;
    colorbar 
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    grid on;

    subplot(2,4,iRow+4) 
    imagesc( deviation(:,:,iSlice,iCh) ) ;
    title('Abs. deviation') ;
    caxis([0 0.1] ) ;
    colorbar ;

    set(gca,'XTick',[])
    set(gca,'YTick',[])
    grid on;
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
%   correspond to isocenter), it is possible that the operator will choose a
%   particular FOV for the following scans which repositions the table by
%   a certain amount ( Field.Hdr.Img.ImaRelTablePosition ), thereby shifting 
%   isocenter (in the patient coordinate system) from Z=0 to Z = Field.Hdr.Img.ImaRelTablePosition.
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

switch Shim.Ref.source
    case 'model'
        Shim.modelreferencemaps( Field ) ;
    case 'data'
        Shim.interpolatetoimggrid( Shim.Field ) ;
end

Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

%% -----
% get the original shim offsets
[f0, g0, s0]                    = Shim.Field.adjvalidateshim() ;
Shim.System.currents            =  [ Shim.System.Specs.converttomultipole( [g0 ; s0] ) ] ; 
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
methods(Hidden)
% =========================================================================
function [ ] = modelreferencemaps( Shim, Field )
%MODELREFERENCEMAPS    Models ideal/nominal Prisma shim reference maps
% 
% Wraps to ShimOpt_SphericalHarmonics.generatebasisfields(), reorders, and
% rescales the basis set to return ideal "shim reference maps" (in units of
% Hz/unit-shim) for the 1st and 2nd order spherical harmonic shims of the
% Siemens Prisma
%
% Usage
%
% [ ] = MODELREFERENCEMAPS( Shim, Field )

assert( nargin == 2 )

[X,Y,Z] = Field.getvoxelpositions() ;

basisFields = ShimOpt_SphericalHarmonics.generatebasisfields( [1:2], X, Y, Z ) ;
% Reorder terms along 4th array dim. in line with Siemens shims: X, Y, Z, Z2, ZX, ZY, X2-Y2, XY
basisFields = reordertosiemens( basisFields ) ; 

scalingFactors = computenormalizationfactors() ;

for iCh = 1 : size( basisFields, 4 ) 
   basisFields(:,:,:,iCh) = scalingFactors(iCh) * basisFields(:,:,:,iCh) ; 
end

Shim.img = basisFields ;
Shim.Hdr = Field.Hdr ;
Shim.Hdr.MaskingImage = true( size( basisFields ) ) ;

Shim.Ref.source = 'model' ;

return;

function [ sh1 ] = reordertosiemens( sh0 )
%REORDERTOSIEMENS
%
%   basisFields1 = REORDERTOSIEMENS( basisFields0 )
%
% basisFields returned from .GENERATEBASISFIELDS() are ordered (along the 4th
% array dimension) like: 
%   Y, Z, X, XY, ZY, Z2, ZX, X2-Y2
%
% REORDERTOSIEMENS orders them in line with Siemens shims: 
%   X, Y, Z, Z2, ZX, ZY, X2-Y2, XY

assert( ( nargin == 1 ) && ( size( sh0, 4 ) == 8 ) )

sh1(:,:,:,1) = sh0(:,:,:,3) ;
sh1(:,:,:,2) = sh0(:,:,:,1) ;
sh1(:,:,:,3) = sh0(:,:,:,2) ;
sh1(:,:,:,4) = sh0(:,:,:,6) ;
sh1(:,:,:,5) = sh0(:,:,:,7) ;
sh1(:,:,:,6) = sh0(:,:,:,5) ;
sh1(:,:,:,7) = sh0(:,:,:,8) ;
sh1(:,:,:,8) = sh0(:,:,:,4) ;

end %reordertosiemens()

function [ scalingFactors ] = computenormalizationfactors()
%COMPUTENORMALIZATIONFACTORS
%
%  scalingFactors = computenormalizationfactors()
%
%  returns a vector of scalingFactors to apply to the (properly reordered)
%  ideal 1st+2nd order spherical harmonic basis returned from GENERATEBASISFIELD
%  to scale the terms as "shim reference maps" in units of Hz/unit-shim
% -----
% Gx, Gy, and Gz should yield 1 micro-T of field shift per metre
% equivalently, 0.042576 Hz/mm
%
% 2nd order terms should yield 1 micro-T of field shift per metre-squared
% equivalently, 0.000042576 Hz/mm^2

%% ------
% create basis on small 3x3x3 mm^3 isotropic grid
[XIso, YIso, ZIso] = meshgrid( [-1:1], [-1:1], [-1:1] ) ;

sh = ShimOpt_SphericalHarmonics.generatebasisfields( [1:2], XIso, YIso, ZIso ) ;
% Reorder terms along 4th array dim. in line with Siemens shims: X, Y, Z, Z2, ZX, ZY, X2-Y2, XY
sh = reordertosiemens( sh ) ; 

nChannels      = size( sh, 4) ; % = 8
scalingFactors = zeros( nChannels, 1 ) ;

% indices of reference positions for normalization: 
iX1   = find( ( XIso == 1 ) & ( YIso == 0 ) & ( ZIso == 0 ) ) ;
iY1   = find( ( XIso == 0 ) & ( YIso == 1 ) & ( ZIso == 0 ) ) ;
iZ1   = find( ( XIso == 0 ) & ( YIso == 0 ) & ( ZIso == 1 ) ) ;

iX1Z1 = find( ( XIso == 1 ) & ( YIso == 0 ) & ( ZIso == 1 ) ) ;
iY1Z1 = find( ( XIso == 0 ) & ( YIso == 1 ) & ( ZIso == 1 ) ) ;
iX1Y1 = find( ( XIso == 1 ) & ( YIso == 1 ) & ( ZIso == 0 ) ) ;

% order the reference indices like the sh field terms 
iRef = [iX1 iY1 iZ1 iZ1 iX1Z1 iY1Z1 iX1 iX1Y1]' ;

%% ------
% scaling:
% 1st order terms yield 1 micro-T of field shift per m (i.e 0.042576 Hz/mm )
% 2nd order terms yield 1 micro-T of field shift per m^2 (i.e 0.000042576 Hz/mm^2 )

% distance from iso/origin to adopted reference point [units: mm]
r = [1 1 1 1 sqrt(2) sqrt(2) 1 sqrt(2)] ;

% invert polarity of certain terms:
sh(:,:,:,[2,3,5,8]) = -sh(:,:,:,[2,3,5,8] ) ;

orders = [1 1 1 2 2 2 2 2] ;

for iCh = 1 : nChannels
    field = sh(:,:,:,iCh) ;
    scalingFactors(iCh) = 42.576*( ( r(iCh) * 0.001 )^orders(iCh) )/field( iRef( iCh ) ) ;
end

end %computenormalizationfactors()

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true, Hidden=true)
% =========================================================================
function [ Params ] = assigndefaultparameters( Params, Specs )
%ASSIGNDEFAULTPARAMETERS  
% 
% Params = ASSIGNDEFAULTPARAMETERS( Params, Specs )
% 
% Add default parameters fields to Params without replacing values (unless empty)

DEFAULTS.shimReferenceMaps       = 'model' ;
DEFAULTS.pathToShimReferenceMaps = [ shimbindir() Specs.Id.systemName ] ;
        
Params = assignifempty( Params, 'shimReferenceMaps', DEFAULTS.shimReferenceMaps ) ;

switch Params.shimReferenceMaps
    case 'calibrate'
        today = datestr( now, 30 ) ;
        today = today(1:8) ; % ignore the time of the day
        Params.pathToShimReferenceMaps = [ shimbindir() Specs.Id.systemName '_' today ] ;
    case 'load'
        Params.shimReferenceMaps = DEFAULTS.pathToShimReferenceMaps ;
    case 'model'
        ; % do nothing
end                        

end
% =========================================================================
function [ img, Hdr, Interpolant ] = calibratereferencemaps( Params )
%CALIBRATEREFERENCEMAPS


% Wraps to .mapdbdi( ) and writes output to disk
% 
% [ img, Hdr ] = CALIBRATEREFERENCEMAPS( Params )

Params = ShimOpt_HGM_Prisma.declarecalibrationparameters( Params ) ;

[ img, Hdr, Interpolant ] = ShimOpt_HGM_Prisma.mapdbdi( Params ) ;

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
% seriesName can be 'sfm_004p'

switch seriesName

    case 'sfm_004p'
        % 2 columns: [ MAG | PHASE ] ;
        dataLoadDirectories    = cell( 2, 2, 8 ) ;

        %% ------
        % channel 1 : "A11" (x-gradient)
        dataLoadDirectories{1,1,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A11+50_0004/' ;
        dataLoadDirectories{1,2,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A11+50_0005/' ;
        dataLoadDirectories{2,1,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A11-50_0006/' ;
        dataLoadDirectories{2,2,1} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A11-50_0007/' ;

        %% ------
        % channel 2 : "B11" (y-gradient)
        dataLoadDirectories{1,1,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B11+50_0008/' ;
        dataLoadDirectories{1,2,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B11+50_0009/' ;
        dataLoadDirectories{2,1,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B11-50_0010/' ;
        dataLoadDirectories{2,2,2} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B11-50_0011/' ;

        %% ------
        % channel 3 : "A01" (z-gradient)
        dataLoadDirectories{1,1,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A10+50_0012/' ;
        dataLoadDirectories{1,2,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A10+50_0013/' ;
        dataLoadDirectories{2,1,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A10-50_0014/' ;
        dataLoadDirectories{2,2,3} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A10-50_0015/' ;

        %% ------
        % channel 4 : "A20"
        dataLoadDirectories{1,1,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A20+600_0016/' ;
        dataLoadDirectories{1,2,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A20+600_0017/' ;
        dataLoadDirectories{2,1,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A20-600_0018/' ;
        dataLoadDirectories{2,2,4} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A20-600_0019/' ;

        %% ------
        % channel 5 : "A21"
        dataLoadDirectories{1,1,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A21+600_0020/' ;
        dataLoadDirectories{1,2,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A21+600_0021/' ;
        dataLoadDirectories{2,1,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A21-600_0022/' ;
        dataLoadDirectories{2,2,5} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A21-600_0023/' ;

        %% ------
        % channel 6 : "B21"
        dataLoadDirectories{1,1,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B21+600_0024/' ;
        dataLoadDirectories{1,2,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B21+600_0025/' ;
        dataLoadDirectories{2,1,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B21-600_0026/' ;
        dataLoadDirectories{2,2,6} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B21-600_0027/' ;

        %% ------
        % channel 7 : "A22"
        dataLoadDirectories{1,1,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A22+600_0028/' ;
        dataLoadDirectories{1,2,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A22+600_0029/' ;
        dataLoadDirectories{2,1,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A22-600_0030/' ;
        dataLoadDirectories{2,2,7} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_A22-600_0031/' ;

        %% ------
        % channel 8 : "B22"
        dataLoadDirectories{1,1,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B22+600_0032/' ;
        dataLoadDirectories{1,2,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B22+600_0033/' ;
        dataLoadDirectories{2,1,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B22-600_0034/' ;
        dataLoadDirectories{2,2,8} = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/sfm_004p/GRE_FM_B22-600_0035/' ;
    
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

Specs             = ShimSpecs_HGM_Prisma() ;
Params.nChannels  = Specs.Amp.nActiveChannels ;

Params.nCurrents  = 2 ;

% for multi-coil systems, entries of Params.currents need to be manually
% defined.  for the host system, the 'currents' (shim coefficients in
% multi-pole units) can be retrieved directly from the dicoms 
% See ShimOpt_HGM_Prisma.mapdbdi( )
Params.currents   = zeros( Params.nChannels, Params.nCurrents ) ; 

% image data directories
if ~myisfieldfilled( Params, 'dataLoadDirectories' )
    Params.dataLoadDirectories = ShimOpt_HGM_Prisma.getdataloaddirectories('sfm_004p') ;
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
        shimValues = ShimSpecs_HGM_Prisma.converttomultipole( [g0 ; s0] ) ; 
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
Hdr.MaskingImage = mask ;

end
% =========================================================================

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
