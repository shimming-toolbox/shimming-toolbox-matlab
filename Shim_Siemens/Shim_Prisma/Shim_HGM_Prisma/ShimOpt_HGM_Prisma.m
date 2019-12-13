classdef ShimOpt_HGM_Prisma < ShimOpt_Prisma
%SHIMOPT_HGM_PRISMA - Shim Optimization for Prisma @ HGM 
%     
% ShimOpt_HGM_Prisma is a ShimOpt_Prisma subclass. 
% See ShimOpt documentation for usage.
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt_HGM_Prisma( varargin )
%SHIMOPT - Shim Optimization

Shim.System.Specs    = ShimSpecs_HGM_Prisma();
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Params = ShimOpt_Prisma.assigndefaultparameters( Params, Shim.System.Specs ) ;

switch Params.shimReferenceMaps
    case 'calibrate'
        [ Shim.img, Shim.Hdr, Shim.Interpolant ] = ShimOpt_HGM_Prisma.calibratereferencemaps( Params ) ;
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
