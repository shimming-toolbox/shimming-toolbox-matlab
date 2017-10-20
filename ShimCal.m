classdef (Abstract) ShimCal < MaRdI
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
%           Shim reference maps (4d array) [units: Hz/ampere]
%
%       .Hdr
%           Info re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    Tracking
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimUse
%    ShimTest 
%     
% 
% ShimCal is a MaRdI subclass [ShimCal < MaRdI]
%
% =========================================================================
% Updated::20170206::ryan.topfer@polymtl.ca
% =========================================================================

properties

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCal( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

end
% =========================================================================
function Params = declarecalibrationparameters20160809( )
%DECLARECALIBRATIONPARAMETERS20160809 
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)


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

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ FitResults ] = fitfieldtocurrent( field, currents, mask, Parameters )
%FITFIELD2CURRENT 
%
% Magnitude-weighted least square regression of field (unwrapped phase) to coil current.
%   
% FITFIELDTOCURRENT( field, currents, mask, Parameters )


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
	disp(voxel/nVoxels)
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
% if # fieldMaps == # curents == 2
%   dBdI = mask .* diff(fieldMaps)/diff(currents);
% else
%   performs linear fitting.
%
% Params has fields...

ShimUse.display('Fitting dB/dI...' )

Params.isDisplayingProgress = true;

dBdI  = [] ; % change in field per change in shim current

nCurrents       = length( currents ) ; 
nInputFieldMaps = size( fieldMaps, 4 ) ;

assert( nCurrents == nInputFieldMaps )

if nInputFieldMaps == 2

    dB   = fieldMaps(:,:,:, 2) - fieldMaps(:,:,:, 1) ;
    dI   = currents(2) - currents(1) ;
    dBdI = mask .* (dB/dI) ;

else

    FitResults = ShimCal.fitfieldtocurrent( fieldMaps, currents, mask, Params ) ;

    dBdI = FitResults.dBdI ;
end


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
%   .isHardCodingEchoTime
%   
%   .EchoTimes
%       uses the difference of Params.EchoTimes rather than the DICOM header to scale 
%       phase to frequency
%
%   .Filtering.isFiltering
%       low-pass filtering? [logical]
%
%   .Filtering.filterRadius
%       radius [in mm] of RESHARP kernel. See help MaRdI.extractharmonicfield() for
%       additional parameters

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

