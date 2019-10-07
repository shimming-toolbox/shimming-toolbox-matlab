classdef ShimOpt_IUGM_Prisma_fit < ShimOpt_Prisma
%SHIMOPT_IUGM_PRISMA_FIT - Shim Optimization for Prisma-fit @ UNF 
%     
% ShimOpt_IUGM_Prisma_fit is a ShimOpt_Prisma subclass. 
% See ShimOpt documentation for usage.
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

Shim.System.Specs    = ShimSpecs_IUGM_Prisma_fit();
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Params = ShimOpt_Prisma.assigndefaultparameters( Params, Shim.System.Specs ) ;

switch Params.shimReferenceMaps
    case 'calibrate'
        [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimSpecs_IUGM_Prisma_fit.calibratereferencemaps( Params ) ;
    case 'model'  
        ; % do nothing until setoriginalfield()
    otherwise 
       [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt.loadshimreferencemaps( Params.shimReferenceMaps ) ; 
        
        Shim.Ref.img = Shim.img ;
        Shim.Ref.Hdr = Shim.Hdr ;
end

if ~isempty( Field ) 
    Shim.setoriginalfield( Field ) ;
end

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Static=true, Hidden=true)
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

end
