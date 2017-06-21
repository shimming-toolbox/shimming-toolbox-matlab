classdef (Abstract) ShimOpt < MaRdI 
%SHIMOPT - Shim Optimization
%
% .......
% 
% Usage
%
% Shim = ShimOpt( )
% Shim = ShimOpt( Params )
% Shim = ShimOpt( Params, Field )
% 
% Inputs
%
%   Field [MaRdI-type object]
%       The original field to be shimmed. (Optional.)
%
%   Params can have the following fields
%
%       .pathToShimReferenceMaps
%           file path to .mat containing shim reference maps (ie basis fields) &
%           .Hdr info
%
%       .TrackerSpecs 
%           parameters struct for ProbeTracking(). 
%           See HELP ProbeTracking() for more information.
%       
%       .isInterpolatingReferenceMaps [== true (default) OR false] 
%           if true,
%           when ShimOpt() is called with the Field argument, interpolation of the
%           reference maps onto the grid of the Field image is done automatically.
%
% Outputs
%
%   Shim contains fields
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info Re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%       .Field
%           Object of type MaRdI pertaining to field distribution to be shimmed
%
%       .Model
%           .currents  
%               Optimal shim current vector (i)
%               [units A]
%           .field     
%               Optimal shim field from projection of i onto reference maps (Ai)
%               [units Hz]
%           .couplingCoefficients
%               For realtime shimming, relates vector relating field to tracker
%               measurement (.Tracker.Data.p)
%               [units: Hz/Pa (Probe), or Hz/rad (Nav) ]
%           .dcCurrentsOffsets
%               For realtime shimming, vector of "y-intercept" currents 
%               (i.e. currents when Tracker.Data.p = 0)
%               [units A]
%
%       .Tracker
%           Object of type Tracking (e.g ProbeTracking, NavTracking)
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    Tracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%
% ShimOpt is a MaRdI subclass [ShimOpt < MaRdI]
%     
% =========================================================================
% Updated::20170214::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
%
% .....
%   .Tracker
%
%   Should be optional // ShimOpt constructor should not 
%   necessarily instantiate a ProbeTracking type object, eliciting a msg
%   "Error: Device not found..." if it isn't connected.
%
% .....
% ASSESSSHIMVOLUME()
%   write function to compare spatial support of shim reference maps to the
%   voxel positions of the field map to be shimmed.
%
% .....
% EXTENDHARMONICFIELD()
%   Write function to check if field map is exists + is reasonable over
%   the shim voi. if some portion is missing, a harmonic extension could
%   be performed to fill in some of the missing field.
% 
% .....
% SETORIGINALFIELD()
%   This will overwrite Shim.Field.Hdr.MaskingImage --- may not be desired/
%   expected by user. Avoid or issue warning?
%
% =========================================================================

properties
    Field ; % object of type MaRdI
    Model ;
    Tracker ; % object of type Tracking (e.g. ProbeTracking)
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( Params, Field )
%SHIMOPT - Shim Optimization


if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt.assigndefaultparameters( Params ) ;

Shim.img = [];
Shim.Hdr = [];

% .......
% Load shim basis if provided 
if ~isempty(Params.pathToShimReferenceMaps)
    ShimUse.display(['\n Preparing for shim ...  \n\n'...
            'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

    load( Params.pathToShimReferenceMaps ) ;
    % Loads .mat containing Shim struct
    % which has fields
    % Shim.img              - linear dB/dI 'current-to-field' opterator
    % Shim.Hdr              - defines info like voxel locations 
    % Shim.Hdr.MaskingImage - defines spatial support of reference maps
end

Shim.Model = [ ] ; 
Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 

if (nargin == 2) && (~isempty(Field))
    
    Shim = Shim.setoriginalfield( Field ) ;
    
    if Params.isInterpolatingReferenceMaps
        
        Shim = interpolatetoimggrid( Shim, Field ) ;
        Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

    end

else

    Shim.Field = [ ] ; % user can assign later via Shim.setoriginalfield() 
end

end
% =========================================================================
function [] = delete( Shim )
%DELETE  
% 
% DELETE( Shim )
% 
% Destructor. Calls Probe.deletecomport( ) 

if ~isempty( Shim.Tracker )
    Shim.Tracker.delete();
end

clear Shim ;

end
% =========================================================================
function Shim = calibraterealtimeupdates( Shim, Params )
%CALIBRATEREALTIMEUPDATES
% 
% CALIBRATEREALTIMEUPDATES asks user to select the median Tracker measurement
% from the measurement logs corresponding to the inspired & expired field maps.
% From these, and the associated optimal currents for the 2 respiratory states,
% the following function calls are made:
%
% Shim.Opt.setcouplingcofficients()
% Shim.Opt.setdccurrentoffsets()
% Shim.Opt.setupdateoperator()
%
% Shim = CALIBRATEREALTIMEUPDATES( Shim, Params ) 
%
%   Params.
%       .Inspired
%           .currents
%           .measurementLog
%
%       .Expired
%           .currents
%           .measurementLog
close all; 
ShimUse.display( ['\n ------- \n ' ...
    'Determine median measurement over inspired apnea : \n \n '] ) ;

Params.Inspired.medianP = ...
    Shim.Tracker.userselectmedianmeasurement( Params.Inspired.measurementLog ) ; 

ShimUse.display( ...
    ['Median measurement : ' num2str( Params.Inspired.medianP )] ) ;

ShimUse.display( ['\n ------- \n ' ...
    'Determine median measurement over expired apnea : \n \n '] ) ;

Params.Expired.medianP = ...
    Shim.Tracker.userselectmedianmeasurement( Params.Expired.measurementLog ) ; 

ShimUse.display( ...
    ['Median measurement : ' num2str( Params.Expired.medianP )] ) ;

Shim.setcouplingcoefficients( ...
    Params.Inspired.currents, Params.Expired.currents, ...
    Params.Inspired.medianP, Params.Expired.medianP ) ;

Shim.setdccurrentoffsets( ...
    Params.Inspired.currents, Params.Expired.currents, ...
    Params.Inspired.medianP, Params.Expired.medianP ) ;

ShimUse.display( ['\n ------- \n ' ...
    'Optimal DC current offsets (in amperes): ' num2str(Shim.Model.dcCurrentOffsets') '\n \n '] ) ;

Shim.setupdateoperator() ;

end
% =========================================================================
function Shim = setcouplingcoefficients( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pInspired, pExpired )
%SETCOUPLINGCOEFFICIENTS
%
% Shim = SETCOUPLINGCOEFFICIENTS( Shim, iIn, iEx, pIn, pEx )
%
% iIn & iEx are vectors of optimal currents for inspired and expired fields.
% pIn & pEx the associated pressures (scalars)
%
% Sets Shim.Model.couplingCoefficients

A = Shim.getshimoperator ;

Shim.Model.couplingCoefficients = ...
    A*(currentsInspired - currentsExpired)/(pInspired - pExpired) ;

end
% =========================================================================
function Shim = setdccurrentoffsets( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pInspired, pExpired )
%SETDCCURRENTOFFSETS
% 
% Compute and set optimal shim DC current offsets (bias)
%
% Shim = SETDCCURRENTOFFSET( Shim, iIn, iEx, pIn, pEx  )
%
% Sets Shim.Model.dcCurrentsOffsets

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

shimFieldOffset = ...
    A*(currentsExpired*pInspired - currentsInspired*pExpired)/(pInspired - pExpired);

X  = (A'*M'*M*A)\A'*M'*M ;
Shim.Model.dcCurrentOffsets = X*shimFieldOffset ;

end
% =========================================================================
function Shim = interpolatetoimggrid( Shim, Field )
%INTERPOLATETOIMGGRID 
%
% Shim = INTERPOLATETOIMGGRID( Shim, Field )
%
% Interpolates Shim.img (reference maps) to the grid (voxel positions) of
% MaRdI-type Img
% 
% i.e.
%
%   [X,Y,Z] = Field.getvoxelpositions ;
%   Shim.resliceimg( X, Y, Z ) ;

[X,Y,Z]      = Field.getvoxelpositions ;
[X0, Y0, Z0] = Shim.getvoxelpositions ;

% check if voxel positions already happen to coincide
% if all( size(X) == size(X0) )
%     if ~all( X(:)==X0(:)) || ~
Shim.resliceimg( X, Y, Z ) ;

end
% =========================================================================
function Shim = setoriginalfield( Shim, Field )
%SETORIGINALFIELD 
%
% Shim = SETORIGINALFIELD( Shim, Field )
%
% Sets Shim.Field
%
% Field is a MaRdI type object with .img in Hz

Shim.Field = Field ;

end
% =========================================================================
function Shim = setshimvolumeofinterest( Shim, mask )
%SETSHIMVOLUMEOFINTEREST 
% 
% Shim = SETSHIMVOLUMEOFINTEREST( Shim, mask )
%
% Sets Shim.Field.Hdr.MaskingImage
%
% mask is a binary image (with the same dimensions as Shim.Field.img) of 
% the desired shim region.

Shim.Field.Hdr.MaskingImage = mask ;

end
% =========================================================================
function Shim = setforwardmodelfield( Shim )
% SETFORWARDMODELFIELD
%
% Shim = SETFORWARDMODELFIELD( Shim ) ;
%
% Predicts shim field (output: Shim.Model.field) for given set of currents 
% (input: Shim.Model.currents)
    
A = Shim.getshimoperator() ;

Shim.Model.field = reshape( A*Shim.Model.currents, size( Shim.Field.img ) ) ;

end
% =========================================================================
function UO = getupdateoperator( Shim )
% GETUPDATEOPERATOR
%
% UO = GETUPDATEOPERATOR( Shim ) ;
%
%   where UO * Shim.Tracker.Data.p = currentsUpdate

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;
c = Shim.Model.couplingCoefficients ;

X  = (A'*M'*M*A)\A'*M'*M ;
UO = X*c ;

end
% =========================================================================
function mask = getvaliditymask( Shim, Params, Field1, Field2 )
%GETVALIDITYMASK 
%
% Returns binary mask - TRUE where field values are well defined and within 
% the expected range and where (interpolated) shim reference maps are well 
% defined.
%
% mask = GETVALIDITYMASK( Shim )
% mask = GETVALIDITYMASK( Shim, Params )
% mask = GETVALIDITYMASK( Shim, Params, Field1 )
% mask = GETVALIDITYMASK( Shim, Params, Field1, Field2 )
% 
% Field1/2 are MaRdI-type objects and may correspond to 'Inspired' and 
% 'Expired' fields.
%
% .......................
%   
% The following Params.fields are supported
%
% .maxAbsField 
%   maximum absolute voxel value assumed to represent an accurate field
%   measurement. Voxels with abs-values greater than this might stem from
%   errors in the unwrapping.  [default: 500 Hz]
%
% .maxFieldDifference
%   maximum absolute voxel-wise difference assumed to be valid between Field1 &
%   Field2 (e.g. 'inspired field' vs. 'expired field') [default: 150 Hz] (See
%   Verma T, Magn Reson Med, 2014)

DEFAULT_MAXABSFIELD        = 500 ;
DEFAULT_MAXFIELDDIFFERENCE = 150 ;

if nargin < 2
    Params.dummy = []
end

if ~myisfield( Params, 'maxAbsField' ) || isempty( Params.maxAbsField ) 
    Params.maxAbsField = DEFAULT_MAXABSFIELD ;
end

if ~myisfield( Params, 'maxFieldDifference' ) || isempty( Params.maxFieldDifference ) 
    Params.maxFieldDifference = DEFAULT_MAXFIELDDIFFERENCE ;
end

mask = Shim.getshimsupport() ;

if nargin >= 3

    mask = mask & logical(Field1.Hdr.MaskingImage) ;
    mask = mask & ( abs(Field1.img) <= Params.maxAbsField ) ;

end

if nargin >= 4

    mask = mask & logical(Field2.Hdr.MaskingImage) ;
    mask = mask & ( abs(Field2.img) <= Params.maxAbsField ) ;

    mask = mask & ( abs( Field1.img - Field2.img ) <= Params.maxFieldDifference ) ;
end

end
% =========================================================================
function Shim = setupdateoperator( Shim )
% SETUPDATEOPERATOR
%
% Shim = SETUPDATEOPERATOR( Shim ) ;
%
% Calls Shim.getupdateoperator() to set field Shim.Model.updateOperator

Shim.Model.updateOperator = Shim.getupdateoperator() ; 

end
% =========================================================================
function M = gettruncationoperator( Shim )
% GETTRUNCATIONOPERATOR
%
% M = GETTRUNCATIONOPERATOR( Shim ) ;
%
% Truncation Shim.Field.Hdr.MaskingImage operator (e.g. M*b, 'picks out' the
% VOI voxels from vector b)

nVoxelsImg = numel( Shim.Field.Hdr.MaskingImage ) ;
nVoxelsVoi = nnz( Shim.Field.Hdr.MaskingImage ) ;

indicesVoi = find( Shim.Field.Hdr.MaskingImage(:) ) ;

M = sparse( [1:nVoxelsVoi], indicesVoi, ones([nVoxelsVoi 1]), nVoxelsVoi, nVoxelsImg ) ;

end
% =========================================================================
function A = getshimoperator( Shim )
% GETSHIMOPERATOR
%
% A = GETSHIMOPERATOR( Shim ) ;
%
%   where A * vectorOfShimCurrents = shimField

nVoxelsImg      = Shim.getnumberofvoxels() ;
nActiveChannels = Shim.getnactivechannels() ;

A = zeros( nVoxelsImg, nActiveChannels ) ; 

for channel = 1 : nActiveChannels
    A(:, channel) = reshape( Shim.img(:,:,:, channel), [nVoxelsImg 1] ) ;
end

end
% =========================================================================
function nActiveChannels = getnactivechannels( Shim )
%GETNACTIVECHANNELS 
%
% Returns number of active shim channels
%
% nActiveChannels = GETNACTIVECHANNELS( Shim ) ;
%
% nActiveChannels = size( Shim.img, 4 ) ;

nActiveChannels = size( Shim.img, 4 ) ;

end
% =========================================================================
function shimSupport = getshimsupport( Shim )
% GETSHIMSUPPORT
%
% shimSupport = GETSHIMSUPPORT( Shim ) ;
%
%   shimSupport is a logical map over the grid (voxel positions) defined by
%   Shim.img of where the shim reference maps have well defined values.

shimSupport = sum(abs(Shim.img),4) > Shim.getnactivechannels()*eps  ;

end
% =========================================================================
function currents = computerealtimeupdate( Shim )
% COMPUTEREALTIMEUPDATE
% 
% Usage
%
% currents = COMPUTEREALTIMEUPDATE( Shim )
%
%   currents = Shim.Model.dcCurrentOffsets + Shim.Model.updateOperator * Shim.Tracker.Data.p(end) ; 

currents = Shim.Model.dcCurrentOffsets + ...
        Shim.Model.updateOperator * Shim.Tracker.Data.p(end) ; 

end
% =========================================================================
function [predictedField, Stats] = predictshim( Shim )
%PREDICTSHIM 
%
% [predictedField, Stats] = PREDICTSHIM( Shim ) ;


mask = logical( Shim.Field.Hdr.MaskingImage ) ;

predictedField = mask .* ( Shim.Field.img + Shim.Model.field ) ;

Stats.Mean = [] ; 
Stats.Mean.predicted          = mean( predictedField( mask ) ) ;
Stats.Mean.original           = mean( Shim.Field.img( mask ) ) ;
Stats.Mean.percentImprovement = ...
    100*(1-Stats.Mean.predicted/Stats.Mean.original) ;

Stats.Deviation = [] ;
Stats.Deviation.predicted          = std( predictedField( mask ) ) ;
Stats.Deviation.original           = std( Shim.Field.img( mask ) ) ;
Stats.Deviation.percentImprovement = ...
    100*(1-Stats.Deviation.predicted/Stats.Deviation.original) ;


end
% =========================================================================
function [f0, f0Voi, f0VoiShimmed] = optimizelarmor( Shim, voi )
%OPTIMIZELARMOR 
%

voi = logical( voi ) ;

predictedField = Shim.predictshim() ;

f0           = Shim.Field.Hdr.ImagingFrequency *(1E6) ; % [units: Hz]

f0Voi        = f0 + mean( Shim.Field.img( voi ) ) ;
f0VoiShimmed = f0 + mean( predictedField( voi ) ) ;

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
methods(Static=true)

% =========================================================================
function [Field, Extras] = mapfield( ImgArray, Params, ObjectiveImg )
% MAPFIELD
%
% Field = MAPFIELD( ImgArray ) 
% Field = MAPFIELD( ImgArray, Params ) 
% Field = MAPFIELD( ImgArray, Params, ObjectiveImg ) 
% 
% ImgArray --- a cell array where the 1st dimension (i.e. row) corresponds to
% echo number and the index along the second dimension (i=1,2) corresponds to
% Magnitude & Wrapped Phase respectively (each a MaRdI-type object) 
%
% if nEchoes == 1 
%   (i.e. ImgArray = { Mag, Phase } )
%   It is assumed that Phase is in fact a phase-difference image.
%
% if nEchoes == 2
%   (i.e. ImgArray = { MagEcho1, PhaseEcho1 ; MagEcho2, PhaseEcho2 } )
%
% ObjectiveImg --- a MaRdI-type object *example* image (e.g. EPI volume!) 
% that, most importantly, has its voxels positioned precisely where the field
% information/shim is desired. That is, if a spherical harmonic fitting of
% the field is performed (see below) then the fitted-field will be interpolated
% at the voxel positions from ObjectiveImg. 
%
% What is the immediate purpose of this? To enabled gapped/partial gre_field_map
% acquisitions in order to reduce the acquisition time for real-time shim training
% i.e. these generally involve the subject holding their breath, so, for
% feasibility + patient comfort, the acquisitions should be as brief as
% possible.
%
% Params --- *must* contain the following fields
%
% if nEchoes == 1
%   .echoTimeDifference [units: ms]
%       
%
% Params --- may contain the following fields
%
%   .mask
%       binary array indicating phase region to be unwrapped 
%       [default: formed by thresholding magnitude images > Params.threshold)
%
%   .threshold  
%   (as a fraction (<1) of max measured magnitude intensity)
%   Determines the phase region to be unwrapped (i.e. areas of low signal are
%   ignored) 
%       [default: 0.01] 
%
%   .isFilteringField 
%       Applies harmonic ('RESHARP', i.e. low-pass) filtering & returns filtered field map.
%       [default: false]
%   
%   .isFittingSphericalHarmonics 
%       Fits field map to spherical harmonic basis set & returns the fitted field map.
%       See doc MaRdI.extractharmonicfield() for relevant Params.
%       [default: false]
%
%   .ordersToGenerate
%       Orders of SH incorporated into fitting 
%       [default: [0:1:8]]

DEFAULT_ISCORRECTINGPHASEOFFSET        = true ; % deprecated
DEFAULT_ISUNWRAPPINGECHOESINDIVIDUALLY = false ; % deprecated
DEFAULT_ISFILTERINGFIELD               = false ;
DEFAULT_THRESHOLD                      = 0.01 ;

DEFAULT_ISFITTINGSPHERICALHARMONICS    = false ;
DEFAULT_ORDERSTOGENERATE               = [0:8] ;

assert( (nargin >= 1) && ~isempty( ImgArray ) ) ;

if ~myisfield( Params, 'isCorrectingPhaseOffset' ) || isempty( Params.isCorrectingPhaseOffset ) 
    Params.isCorrectingPhaseOffset = DEFAULT_ISCORRECTINGPHASEOFFSET ;
end

if ~myisfield( Params, 'isUnwrappingEchoesIndividually' ) || isempty( Params.isUnwrappingEchoesIndividually ) 
    Params.isUnwrappingEchoesIndividually = DEFAULT_ISUNWRAPPINGECHOESINDIVIDUALLY ;
end

if ~myisfield( Params, 'isFilteringField' ) || isempty( Params.isFilteringField ) 
    Params.isFilteringField = DEFAULT_ISFILTERINGFIELD ;
end

if ~myisfield( Params, 'threshold' ) || isempty( Params.threshold ) 
    Params.threshold = DEFAULT_THRESHOLD ;
end

if ~myisfield( Params, 'isFittingSphericalHarmonics' ) || isempty( Params.isFittingSphericalHarmonics ) 
    Params.isFittingSphericalHarmonics = DEFAULT_ISFITTINGSPHERICALHARMONICS ;
end

if ~myisfield( Params, 'ordersToGenerate' ) || isempty( Params.ordersToGenerate ) 
    Params.ordersToGenerate = DEFAULT_ORDERSTOGENERATE ;
end

Extras = [] ;

% .......
nEchoes = size( ImgArray, 1 ) ;

% -------
% define spatial support for unwrapping
if ~myisfield( Params, 'mask' ) || isempty( Params.mask )

    Params.mask = ones( ImgArray{1,1}.getgridsize ) ;

    for iEcho = 1 : nEchoes 

        Params.mask = Params.mask .* ( ImgArray{ iEcho, 1 }.img > Params.threshold ) ;

    end

end


if nEchoes == 1
    
    PhaseDiff = ImgArray{ 1, 2 } ;
    PhaseDiff.Hdr.EchoTime = Params.echoTimeDifference ;

else
    
    % -------
    % phase difference image via complex division
    PhaseDiff      = ImgArray{ 1, 2}.copy() ;

    img            = ImgArray{ 1, 1 }.img .* exp(-i*ImgArray{ 1, 2 }.img) ;
    img(:,:,:,2)   = ImgArray{ 2, 1 }.img .* exp(-i*ImgArray{ 2, 2 }.img) ;

    PhaseDiff.img = angle( img(:,:,:,2) ./ img(:,:,:,1) ) ;

    PhaseDiff.Hdr.EchoTime = ( ImgArray{ 1, 2 }.Hdr.EchoTime - ImgArray{ 2, 2 }.Hdr.EchoTime );
end

% -------
% 3d path-based unwrapping
PhaseDiff.Hdr.MaskingImage = Params.mask ;

PhaseDiff = PhaseDiff.unwrapphase(  ) ;

Field     = PhaseDiff.scalephasetofrequency( ) ;


if Params.isFilteringField

[~,Field] = Field.extractharmonicfield( Params ) ;

end

    if Params.isFittingSphericalHarmonics % UNTESTED
    % -------
    % fit spherical harmonic basis set to input Field 

    % generates basis set, with field positions same as those of input Field 
    Shims = ShimOptSHarmonics( Params, Field ) ;

    % calculate fitting coefficients ('currents')
    Shims = Shims.optimizeshimcurrents( Params ) ;
    
    Extras.FieldResidual = Field.img + Shims.Model.field ;

    Field.img = -Shims.Model.field ;
    
    if (nargin == 3) & ~isempty( ObjectiveImg )
        % Interpolate the field @ VoxelPositions
        %
        % Main purpose: to enable gapped slices in the field map acquisitions
        % for real-time shim training --- by reducing nSlices, acq. time is
        % reduced, & therefore, the necessary duration of the breath hold.

        [X0, Y0, Z0] = Field.getvoxelpositions( ) ; % original
        [X, Y, Z]    = ObjectiveImg.getvoxelpositions( ) ; % final

        % recalculate the set of harmonics at the given voxel positions      
        basisFields = ShimOptSHarmonics.generatebasisfields( Params.ordersToGenerate, X, Y, Z ) ;
        
        % scale each harmonic by the fitted 'currents' (coefficients)
        for iHarmonic = 1 : size( basisFields, 4 ) 
            basisFields(:,:,:, iHarmonic) = Shims.Model.currents(iHarmonic) * basisFields(:,:,:, iHarmonic) ;
        end
        
        Field.img = sum( -basisFields, 4 ) ;
        

        disp( ['Interpolating phase/field mask...' ]) ;
            
        Field.Hdr.MaskingImage = griddata( X0, Y0, Z0, Field.Hdr.MaskingImage, X, Y, Z, 'nearest' ) ;

        % if new positions are outside the range of the original, 
        % interp3/griddata replaces array entries with NaN
        Field.Hdr.MaskingImage( isnan( Field.Hdr.MaskingImage ) ) = 0 ; 

        % -------
        % Update Hdr 
        %
        % Note: the Hdr could probably simply be copied from ObjectiveImg but recalculating the entries 
        % is more general ('extensible') should the future user not have a
        % fully-formed 'ObjectiveImg' set of dicoms, but merely the target
        % voxel positions [X,Y,Z]
        % 
        % That said, the way SliceLocation is updated below may not always be correct.
        % (borrowed from MaRdI.resliceimg() )
        
        Field.Hdr.ImagePositionPatient( 1 ) = X(1) ; 
        Field.Hdr.ImagePositionPatient( 2 ) = Y(1) ;
        Field.Hdr.ImagePositionPatient( 3 ) = Z(1) ;

        %-------
        % Rows 
        Field.Hdr.Rows = size(Field.img, 1) ;

        dx = X(2,1,1) - X(1,1,1) ;
        dy = Y(2,1,1) - Y(1,1,1) ;
        dz = Z(2,1,1) - Z(1,1,1) ;  

        % vertical (row) spacing
        Field.Hdr.PixelSpacing(1) = ( dx^2 + dy^2 + dz^2 )^0.5 ; 

        % column direction cosine (expressing angle btw column direction and X,Y,Z axes)
        Field.Hdr.ImageOrientationPatient(4) = dx/Field.Hdr.PixelSpacing(1) ;
        Field.Hdr.ImageOrientationPatient(5) = dy/Field.Hdr.PixelSpacing(1) ;
        Field.Hdr.ImageOrientationPatient(6) = dz/Field.Hdr.PixelSpacing(1) ;

        %-------
        % Columns 
        Field.Hdr.Columns = size(Field.img, 2) ;       

        dx = X(1,2,1) - X(1,1,1) ;
        dy = Y(1,2,1) - Y(1,1,1) ;
        dz = Z(1,2,1) - Z(1,1,1) ;  

        % horizontal (column) spacing
        Field.Hdr.PixelSpacing(2) = ( dx^2 + dy^2 + dz^2 )^0.5 ;

        % row direction cosine (expressing angle btw column direction and X,Y,Z axes)
        Field.Hdr.ImageOrientationPatient(1) = dx/Field.Hdr.PixelSpacing(2) ;
        Field.Hdr.ImageOrientationPatient(2) = dy/Field.Hdr.PixelSpacing(2) ;
        Field.Hdr.ImageOrientationPatient(3) = dz/Field.Hdr.PixelSpacing(2) ;

        %-------
        % Slices
        Field.Hdr.NumberOfSlices       = size(Field.img, 3) ;
        Field.Hdr.SpacingBetweenSlices = ( (X(1,1,2) - X(1,1,1))^2 + ...
                                           (Y(1,1,2) - Y(1,1,1))^2 + ...
                                           (Z(1,1,2) - Z(1,1,1))^2 ) ^(0.5) ;

        [~, ~, sHat] = Field.getdirectioncosines( ) ;  
        Field.Hdr.SliceLocation = dot( Field.Hdr.ImagePositionPatient, sHat ) ;
    end

end

end
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
% DEFAULT_TRACKERSPECS = [] ;
%
% DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;

DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
DEFAULT_TRACKERSPECS            = [] ;

DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_TRACKERSPECS ;
end

if ~myisfield( Params, 'isInterpolatingReferenceMaps' ) || isempty(Params.isInterpolatingReferenceMaps)
   Params.isInterpolatingReferenceMaps = DEFAULT_ISINTERPOLATINGREFERENCEMAPS ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Abstract)
% =========================================================================
[Shim] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Shim = OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% OPTIMIZESHIMCURRENTS should return the optimal shim currents in the field
% Shim.Model.currents

% =========================================================================
% =========================================================================
end

end
