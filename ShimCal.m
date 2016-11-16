classdef ShimCal < MaRdI
%SHIMCAL - Shim Calibration 
%
% .......
% 
% Usage
%
% Shim = ShimCal( )
%
%   Shim contains fields
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%
% =========================================================================
% Notes
% 
% ShimCal is a MaRdI subclass [ShimCal < MaRdI]
%
% Part of series of classes pertaining to shimming:
%
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimUse
%    ShimTest 
%     
%
% =========================================================================
% Updated::20160829::ryan.topfer@polymtl.ca
% =========================================================================

properties

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCal( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

if nargin < 1
    % Params = ShimCal.declarecalibrationparameters20160809() ;
    % Params = ShimCal.declarecalibrationparameters20161003() ;
    % Params = ShimCal.declarecalibrationparameters20161004() ;
    Params = ShimCal.declarecalibrationparameters20161007() ;
end

% -------
% Process Zero-current data
Mag   = MaRdI( Params.dataLoadDirectories{1} ) ;
Phase = MaRdI( Params.dataLoadDirectories{2} ) ;

%-----
% Define data support (used for all field maps)
Phase.Hdr.MaskingImage = zeros( Phase.getgridsize( ) ) ;
Phase.Hdr.MaskingImage( Params.limitsForUnwrapping(1,1) : (Params.limitsForUnwrapping(1,2)), ...
                   Params.limitsForUnwrapping(2,1) : (Params.limitsForUnwrapping(2,2)), ...
                   Params.limitsForUnwrapping(3,1) : (Params.limitsForUnwrapping(3,2)) ) = 1 ;

% NB: zero-current acquisition assumed to be the first 2 directories (mag, phase)
BackgroundField = ShimCal.preprocessfielddata( Phase, Params ) ;

Shim.Hdr = BackgroundField.Hdr ;
Shim.img = zeros( BackgroundField.getgridsize ) ;

Params.DataDir.Phase = Params.dataLoadDirectories(2) ; % just to fill field/declare variable

for iChannel = 1 : Params.nChannels 
    
    for iCurrent = 1 : (Params.nCurrents-1)
        % Collect the list of directories corresponding to iChannel
        iImg = 4*iChannel + 2*(iCurrent-1) ;
        Params.DataDir.Phase( iCurrent ) = Params.dataLoadDirectories(iImg) ;
    end

    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
    currents = [0 Params.currents(iChannel, 1) Params.currents(iChannel, 2)] ;
    Shim.img(:,:,:, iChannel) = ...
        ShimCal.mapdbdi( BackgroundField.img, currents, Phase.Hdr.MaskingImage, Params ) ;

end 

if Params.Extension.isExtending
    
    disp(['Performing harmonic extension...']) ;
    disp(['Channel 1 of ' num2str(Params.nChannels) ] ) ;

    gridSizeImg = size( Shim.img(:,:,:,1) ) ;
    Params.Extension.voxelSize = Shim.getvoxelsize() ;
    
    maskReduced = (sum(abs(Shim.img),4))~=0 ; % initial spatial support
    mask        = maskReduced ; % extended spatial support
    mask( Params.Extension.limitsForExtension(1,1) : (Params.Extension.limitsForExtension(1,2)), ...
        Params.Extension.limitsForExtension(2,1) : (Params.Extension.limitsForExtension(2,2)), ...
        Params.Extension.limitsForExtension(3,1) : (Params.Extension.limitsForExtension(3,2)) ) = 1 ;
    
    BackgroundField.Hdr.MaskingImage = mask ;

    [ ~, A, M ] = extendharmonicfield( Shim.img(:,:,:, 1) , mask, maskReduced, Params.Extension ) ;

    for iChannel = 1 : Params.nChannels 
        disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] ) ;
        
        reducedField              = maskReduced .* Shim.img(:,:,:, iChannel) ;
        extendedField             = reshape( M'*A*reducedField(:), gridSizeImg ) ;
        Shim.img(:,:,:, iChannel) = extendedField + reducedField ;
    end
    
end


SpineShim = [] ;
SpineShim.img  = Shim.img ;
SpineShim.Hdr  = Shim.Hdr ;
SpineShim.mask = BackgroundField.Hdr.MaskingImage ;
save( Params.filenameSave ,'SpineShim');

end
% =========================================================================
% =========================================================================
% =========================================================================
end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

function [] = acquirecalibrationdata( )
%ACQUIRECALIBRATIONDATA 
%
%  Re-acquiring calibration data on "waterproof phantom" (28 L water) 
%  Experiment performed: 20160809
%
%  Notes
%
%   **This is a demo script to accompany reference field map acquisition 
%   ( i.e. it's not meant to be run as a standalone function).**
%
%   Halfway through the acquisitions I realised I should record the channel
%   output currents rather than simply the assigned test currents (ideally,
%   these would be identical, but owing to the finite (16 bit) precision of 
%   the DACs, and perhaps other factors, they are close, but not identical).
%
%   Therefore I created a variable called measuredCurrents to keep track of 
%   subsequent recordings, and printed the terminal display output to a file
%   ( ~/Projects/Shimming/Static/Calibration/Data/20160809_currentsLog.pdf)
%   Unfortunately, "infinte scrollback" was disabled on iterm2, so, the 
%   measured currents for the first 7.5 channels were lost. (7.5 because 
%   the *positive* measured current for the 8th channel was recovered.)
%
%   To correct for this, I'm substituting the average measured currents for
%   the remaining channels into the absent values for these first 7.5 
%   channels. (See the last section below.)
%

% -------
% Connect to shim amplifier 
ShimsC = ShimCom;
ShimsC.getsystemheartbeat 

channelsToBankKey = ShimsC.getchanneltobankkey ;

testCurrents = [-0.5 0.5]; % [units : A]

iChannel = 0;

measuredCurrents = zeros(24, length(testCurrents) ) ;
% -------
% Perform routine for iChannel = 1 : 24
iChannel = iChannel + 1 

ShimsC.resetallshims()
iTestCurrent = 1 ;

ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), testCurrents(iTestCurrent) ) ;

pause(0.5)
ChannelOutput = ShimsC.getchanneloutput(channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) )

measuredCurrents( iChannel, iTestCurrent ) = ChannelOutput.current ;


ShimsC.resetallshims()
iTestCurrent = 2 ;

ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), testCurrents(iTestCurrent) ) ;

pause(0.5)
ChannelOutput = ShimsC.getchanneloutput(channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) )

measuredCurrents( iChannel, iTestCurrent ) = ChannelOutput.current ;
% ------- 
% Account for missing measured currents & save as .txt (see 'Notes' above)

measuredCurrents(1:8,1) = mean( measuredCurrents(9:end,1) ) ;
measuredCurrents(1:7,2) = mean( measuredCurrents(8:end,2) ) ;

filenameSave = ...
    '~/Projects/Shimming/Static/Calibration/Data/20160809_measuredCurrents.txt' ;

fid = fopen( filenameSave, 'w' ) ;
fprintf( fid, '%f %f \n', measuredCurrents' ) ;
fclose( fid ) ;

end
% =========================================================================
function [] = comparereferencemaps( )
%COMPAREREFERENCEMAPS 
%

Params.pathToShimReferenceMaps = '~/Projects/Shimming/RRI/data/SpineShimReferenceMaps.mat' ; 
ShimsOld = ShimOpt(Params) ;

Params.pathToShimReferenceMaps = '~/Projects/Shimming/RRI/data/SpineShimReferenceMaps20160809.mat' ;
ShimsNew = ShimOpt(Params) ;

[X,Y,Z] = ShimsOld.getvoxelpositions() ;
ShimsNew.resliceimg( X, Y, Z ) ;

NiiOptions.filename = './old_ref_maps' ;
nii( ShimsOld.img, NiiOptions ) ;

NiiOptions.filename = './new_ref_maps' ;
nii( ShimsNew.img, NiiOptions ) ;
end
% =========================================================================
function Params = declarecalibrationparameters20160809( )
%DECLARECALIBRATIONPARAMETERS20160809 
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)
%
% Notes
%
%   This doesn't need to be a ShimCom method, and could be written as a
%   separate file.

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 7.38] ;
Params.nChannels  = 24 ;
Params.nCurrents  = 3 ;

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_012p/20160809_measuredCurrents.txt') ;
Params.currents = fscanf( fileIn, '%f', [2 Params.nChannels ] )' ;
fclose(fileIn);

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_012p/shimcalibrationdata_dirlist_20160809.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition assumed to be the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20160809' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.Extension.isExtending = false ; % harmonic field extrapolation 
% Define data support (used for all field maps)
Params.limitsForUnwrapping = [2 90; 2 28; 1 80] ; 

end
% =========================================================================
function Params = declarecalibrationparameters20161003( )
    
Params = ShimCal.declarecalibrationparameters20160809( );
Params.Filtering.isFiltering  = false ;

Params.filenameSave = ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161003_nofiltering' ;

end
% =========================================================================
function Params = declarecalibrationparameters20161004( )
    
Params = ShimCal.declarecalibrationparameters20160809( );

% Define data support (used for all field maps)
Params.limitsForUnwrapping = [11 80; 3 26; 12 69] ; 

Params.filenameSave = ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161004' ;

end
% =========================================================================
function Params = declarecalibrationparameters20161007( )
    
Params = ShimCal.declarecalibrationparameters20160809( );

% Define data support (used for all field maps)
Params.limitsForUnwrapping = [11 80; 3 26; 12 69] ; 

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.radius      = 10 ;
Params.Extension.isDisplayingProgress = true ;
Params.Extension.limitsForExtension = [10 80; 1 28; 13 70] ;

Params.filenameSave = ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007' ;

end
% =========================================================================
function [ FitResults ] = fitfieldtocurrent( field, currents, mask, Parameters )
%FITFIELD2CURRENT 
%
% Magnitude-weighted least square regression of field (unwrapped phase) to coil current.
%   
%   b_ij is the field measurement at the i-th pixel for the j-th applied current
%   x_ij is the n-th current applied to the i-th pixel
%   a_i0 is the y-intercept of the line describing the i-th pixel
%   a_i1 is its slope


DEFAULT_TOLERANCE 	  = 1E-6 ;
DEFAULT_MAXITERATIONS = 10000 ;
DEFAULT_MU            = 0 ;

if nargin < 4 || isempty(Parameters)
	Parameters.dummy = [] ;
end

if ~myisfield( Parameters, 'tolerance' ) || isempty(Parameters.tolerance)
	Parameters.tolerance = DEFAULT_TOLERANCE ;
end

if ~myisfield( Parameters, 'maxIterations' ) || isempty(Parameters.maxIterations)
	Parameters.maxIterations = DEFAULT_MAXITERATIONS ;
end

if ~myisfield( Parameters, 'mu' ) || isempty(Parameters.mu)
    Parameters.mu = DEFAULT_MU;
end



gridSize  = size( field ) ;
gridSize  = gridSize(1:3) ;
nVoxels   = prod(gridSize) ;
nCurrents = length( currents ) ;



% linear 'current-to-field' operator in Xa = b
X = [ones( [nCurrents 1] ) currents'] ;
X = sparse(X) ;
while( size(X,2) < 2*nVoxels )
	X = blkdiag(X,X) ;
end
X  = X(1: nCurrents*nVoxels, 1:2*nVoxels) ;


M_0 = spalloc( nVoxels, 2*nVoxels, nVoxels ) ; % extracts fieldOffset a_0
M_1 = spalloc( nVoxels, 2*nVoxels, nVoxels ) ; % extracts dBdI a_1

placeHolder = 1 ;
for voxel = 1 : nVoxels
	M_0(voxel, placeHolder ) 	= 1 ;
	M_1(voxel, placeHolder + 1) = 1 ;
	placeHolder = placeHolder + 2 ;
end 

% vectorized field measurements
b = zeros([nCurrents*nVoxels 1]) ;

% W will be used to form a diagonal matrix consisting of (squared) data-reliability weights
W = zeros([nCurrents*nVoxels 1]) ;

for currentIndex = 1 : nCurrents
	b( currentIndex : nCurrents : end ) = ...
		mask(:) .* reshape( field(:,:,:, currentIndex), [nVoxels 1] )  ;    
	
	W( currentIndex : nCurrents : end ) = mask(:) ; 
end

XtWtW = X' * spdiags(W.^2, 0, nCurrents*nVoxels, nCurrents*nVoxels ) ;

if Parameters.mu ~= 0
    % Regularizing operator R
	[Dx, Dy, Dz] = createdifferenceoperators( gridSize, Parameters.voxelSize, 2 ) ;
	L = Dx + Dy + Dz ; 
	R = spdiags( mask(:), 0, size(L,1), size(L,2) ) * L*(M_0 + M_1);

	fitParameters = cgls( XtWtW*X + Parameters.mu*(R'*R), XtWtW*b, zeros( [size(X,2) 1] ), Parameters ) ;
else

	fitParameters = cgls( XtWtW*X, XtWtW*b, zeros( [size(X,2) 1] ), Parameters ) ;
end

FitResults.fieldOffset = M_0 * fitParameters ;
FitResults.fieldOffset = reshape( FitResults.fieldOffset, gridSize ) ;

FitResults.dBdI = M_1 * fitParameters ;
FitResults.dBdI = reshape( FitResults.dBdI, gridSize ) ;


% -------
% Calculate voxelwise R^2
bestFit = zeros( [gridSize nCurrents] ) ;

for iCurrent = 1 : nCurrents
	bestFit(:,:,:, iCurrent) = FitResults.dBdI * currents(iCurrent) +  FitResults.fieldOffset ;
end

sumOfSquaresResiduals    = sum( (field - bestFit).^2, 4 ) ;
sumOfSquaresMeasurements = ( nCurrents - 1 ) * var( field, 0, 4 ) ;

FitResults.rSquared = 1 - sumOfSquaresResiduals./sumOfSquaresMeasurements ;



end
% =========================================================================
function dBdI = mapdbdi( fieldMaps, currents, mask, Params )
%MAPDBDI
% 
% Maps single channel dB/dI (change in field per unit current [Hz/A])  
% 
% dBdI = MAPDBDI( fields, currents, mask, Params ) 
% 
% fieldMaps
%   field maps to be fitted to current
% 
% currents
%   test currents 
% 
% mask
%   binary mask indicating voi over which dB/dI is to be calculated
%
% Params has fields
%
%   .DataDir
%       .Phase 

ShimUse.display('Fitting dB/dI...' )

dBdI  = [] ; % change in field per change in shim current

nCurrents = length( currents ) ;

if isempty( fieldMaps )
    fieldMaps = [] ; % the field maps to be fitted to current
    nInputFieldMaps = 0 ;
else
    nInputFieldMaps = size( fieldMaps, 4 );
    assert( nCurrents >= nInputFieldMaps )
end

for iCurrent = (nInputFieldMaps + 1) : nCurrents 

    Phase = MaRdI( Params.DataDir.Phase{ (iCurrent - nInputFieldMaps) }  ) ;
    Phase.Hdr.MaskingImage = mask ;    

    Field = ShimCal.preprocessfielddata( Phase, Params ) ;

    fieldMaps( :,:,:, iCurrent ) = Field.img ;

end

Params.isDisplayingProgress = true;

FitResults = ShimCal.fitfieldtocurrent( ...
    fieldMaps, currents, Field.Hdr.MaskingImage, Params ) ;

dBdI = FitResults.dBdI ;

end
% =========================================================================
function Field = preprocessfielddata( Phase, Params ) 
%PREPROCESSFIELDDATA
% 
% Performs masking, phase unwrapping, scaling, and filtering.
%
% Field = PREPROCESSFIELDDATA( Phase, Params )
%
% Phase and Field are MaRdI-type objects 
% 
% Params has fields
%
% .isHardCodingEchoTime
% .EchoTimes
%   uses the difference of .EchoTimes rather than the DICOM header to scale 
%   phase to frequency
%
% .Filtering.filterRadius
%   radius [in mm] of RESHARP kernel. See help MaRdI.extractharmonicfield() for
%   additional parameters

%-----
% 3d path-based unwrapping
Phase = Phase.unwrapphase(  ) ;

%-----
% Scale phase to field
if Params.isHardCodingEchoTime
    Phase.Hdr.EchoTime = Params.EchoTimes(2) - Params.EchoTimes(1) ;
end

Field = Phase.scalephasetofrequency( ) ;

%-----
% Harmonic filtering 
if Params.Filtering.isFiltering
    ShimUse.display('Filtering...' )
    [~, Field] = Field.extractharmonicfield( Params.Filtering ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

