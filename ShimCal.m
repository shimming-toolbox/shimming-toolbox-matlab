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
% Updated::20160826::ryan.topfer@polymtl.ca
% =========================================================================

properties

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCal(  )
%SHIMCAL - Shim calibration (i.e. reference maps) 
%

end
% =========================================================================
function Shim = calibratespineshim(  )
% -------
fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')

isHardCodingEchoTime = true ;

% -------
Params.EchoTimes = [4.92 7.38] ;

Params.projectDirectory  = '~/Projects/Shimming/Static/Calibration/' ;
Params.dataLoadDirectory = [Params.projectDirectory '/Data/shim_012p/'] ;
Params.tmpSaveDirectory  = [Params.projectDirectory 'Tmp/'] ;

Params.nChannels  = 24 ;
Params.nCurrents  = 3 ;

Params.limitsForUnwrapping = [2 90; 2 28; 1 80] ; 

field = [] ; % the field maps to be fitted to current
dBdI  = [] ; % change in field per change in shim current

% -------
% Process Zero-current data
Mag   = MaRdI( [Params.dataLoadDirectory '118-gre_field_mapping_0A/echo_7.38/'] ) ;
Phase = MaRdI( [Params.dataLoadDirectory '119-gre_field_mapping_0A/echo_7.38/'] ) ;


%-----
% Define data support
maskForUnwrapping  = zeros( Phase.getgridsize( ) ) ;
maskForUnwrapping( Params.limitsForUnwrapping(1,1) : (Params.limitsForUnwrapping(1,2)), ...
                   Params.limitsForUnwrapping(2,1) : (Params.limitsForUnwrapping(2,2)), ...
                   Params.limitsForUnwrapping(3,1) : (Params.limitsForUnwrapping(3,2)) ) = 1 ;

%-----
% Unwrap phase
Phase.Hdr.MaskingImage = logical( maskForUnwrapping ) ;
Phase = Phase.unwrapphase(  ) ;

%-----
if isHardCodingEchoTime
    Phase.Hdr.EchoTime = Params.EchoTimes(2) - Params.EchoTimes(1) ;
end



Field = Phase.scalephasetofrequency( ) ;

FilteringParams.filterRadius = Field.getvoxelsize() ;
FilteringParams.filterRadius = FilteringParams.filterRadius(1) ;

[Tmp1,Tmp2] = Field.extractharmonicfield( FilteringParams ) ;

Field1.extractharmonicfield( 
[localField, reducedMask ] = resharp( Field.img, ...
                                 Field.Hdr.MaskingImage, ... 
                                 Field.getvoxelsize(), ...
                                 FilteringParams.filterRadius, ...
                                 0) ;

end
% =========================================================================
function [ FitResults ] = fitfieldtocurrent( field, currentsTested, mask, Parameters )
%FITFIELD2CURRENT Magnitude-weighted least square regression of field (unwrapped phase) to coil current.
%   
%   b_ij is the field measurement at the i-th pixel for the j-th applied current
%   x_ij is the n-th current applied to the i-th pixel
%   a_i0 is the y-intercept of the line describing the i-th pixel
%   a_i1 is its slope

%
%
%     1  x_00                      [   a_00         ] = b_00
%     1  x_01                         a_01          b_01
%     1  x_02                          ...
%                                     a_Np0          b_02
%      ...        0           0           a_Np1          ...
%     1  x_0Nc          
%        0          
%        0                                                     b_Nc
%       ...                                           ...
%              1  x_10                                      b_Np0                  
%              1  x_11                                      b_Np1
%        0     1  x_12        0                              b_Np2
%               ...                                           ...
%              1  x_1Nc                                     b_NpNc
%
%
% argmin_a( || W X a - W b ||^2 )
%
% [coilProfile, res]

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



gridSize               = size( field ) ;
gridSize               = gridSize(1:3) ;
numberOfVoxels         = prod(gridSize) ;
numberOfCurrentsTested = length( currentsTested ) ;



% linear 'current-to-field' operator in Xa = b
X = [ones( [numberOfCurrentsTested 1] ) currentsTested'] ;
X = sparse(X) ;
while( size(X,2) < 2*numberOfVoxels )
	X = blkdiag(X,X) ;
end
X  = X(1: numberOfCurrentsTested*numberOfVoxels, 1:2*numberOfVoxels) ;





M_0 = spalloc( numberOfVoxels, 2*numberOfVoxels, numberOfVoxels ) ; % extracts fieldOffsetParameter a_0
M_1 = spalloc( numberOfVoxels, 2*numberOfVoxels, numberOfVoxels ) ; % extracts current2field Parameter a_1

placeHolder = 1 ;
for voxel = 1 : numberOfVoxels
	M_0(voxel, placeHolder ) 	= 1 ;
	M_1(voxel, placeHolder + 1) = 1 ;
	placeHolder = placeHolder + 2 ;
end 

% vectorized field measurements
b = zeros([numberOfCurrentsTested*numberOfVoxels 1]) ;

% W will be used to form a diagonal matrix consisting of (squared) data-reliability weights
W = zeros([numberOfCurrentsTested*numberOfVoxels 1]) ;

for currentIndex = 1 : numberOfCurrentsTested
	b( currentIndex : numberOfCurrentsTested : end ) = ...
		mask(:) .* reshape( field(:,:,:, currentIndex), [numberOfVoxels 1] )  ;    
	
	W( currentIndex : numberOfCurrentsTested : end ) = mask(:) ; 
end

XtWtW = X' * spdiags(W.^2, 0, numberOfCurrentsTested*numberOfVoxels, numberOfCurrentsTested*numberOfVoxels ) ;

if Parameters.mu ~= 0
% regularizing operator
	[Dx, Dy, Dz] = createdifferenceoperators( gridSize, Parameters.voxelSize, 2 ) ;
	L = Dx + Dy + Dz ; 
	R = spdiags( mask(:), 0, size(L,1), size(L,2) ) * L*(M_0 + M_1);

	fitParameters = cgls( XtWtW*X + Parameters.mu*(R'*R), XtWtW*b, zeros( [size(X,2) 1] ), Parameters ) ;
else

	fitParameters = cgls( XtWtW*X, XtWtW*b, zeros( [size(X,2) 1] ), Parameters ) ;
end

FitResults.fieldOffsetParameter = M_0 * fitParameters ;
FitResults.fieldOffsetParameter = reshape( FitResults.fieldOffsetParameter, gridSize ) ;

FitResults.current2FieldParameter = M_1 * fitParameters ;
FitResults.current2FieldParameter = reshape( FitResults.current2FieldParameter, gridSize ) ;

bestFit = zeros( [gridSize numberOfCurrentsTested] ) ;

for current = 1 : numberOfCurrentsTested
	bestFit(:,:,:, current) = FitResults.current2FieldParameter * currentsTested(current) +  FitResults.fieldOffsetParameter ;
end

sumOfSquaresResiduals    = sum( (field - bestFit).^2, 4 ) ;
sumOfSquaresMeasurements = ( numberOfCurrentsTested - 1 ) * var( field, 0, 4 ) ;

FitResults.rSquared = 1 - sumOfSquaresResiduals./sumOfSquaresMeasurements ;

end

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

function [] = acquirecalibrationdata( )
% Acquire calibration data 
%  Re-acquiring calibration data on "waterproof phantom" (28 L water) 
%  Experiment performed: 20160809
%
%     Notes
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

end
% =========================================================================
% =========================================================================

end

