classdef (Abstract) ShimOpt < FieldEval 
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
% ShimOpt is a FieldEval subclass [ShimOpt < FieldEval < MaRdI]
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
%   --> implement   
%
%       bool ProbeTracking.checkconnection( )
%
%   
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
    Field ; % object of type FieldEval
    Model ; % Modeled quantities for shimming
    Params ; % Parameters struct 
    Tracker ; % object of type Tracking (e.g. ProbeTracking)
    ShimmedField; % object of type FieldEval 
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

Shim.Params = ShimOpt.assigndefaultparameters( Params ) ;

Shim.img = [];
Shim.Hdr = [];
Shim.Aux = []; % ShimOpt inherits this from FieldEval... not sure yet what it will be used for

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
function [Results] = assessshim( Shim, voi )
%ASSESSSHIM  
% 
% ASSESSSHIM( Shim )

DEFAULT_ISASSESSINGSHIMMEDFIELD = false;

if nargin == 1
    voi = Shim.Field.Hdr.MaskingImage ;
end

PredictedField = Shim.predictshimmedfield();

Original  = Shim.Field.assessfielddistribution( voi ) ;
Predicted = PredictedField.assessfielddistribution( voi ) ;

Results.volume = Original.volume ; 

Results.Std.Original     = Original.std ;
Results.Std.Predicted    = Predicted.std ;

Results.Std.PredictedImprovement = 100 * (1 - Predicted.std  ./ Original.std  )  ;

Results.Mean.Original    = Original.mean ;
Results.Mean.Predicted   = Predicted.mean ;

Results.Median.Original  = Original.median ;
Results.Median.Predicted = Predicted.median ;

Results.Norm.Original    = Original.norm ;
Results.Norm.Predicted   = Predicted.norm ;

if ~myisfield( Shim, 'ShimmedField' ) || isempty( Shim.ShimmedField )
    
   isAssessingShimmedField = DEFAULT_ISASSESSINGSHIMMEDFIELD ;

else

    isAssessingShimmedField = true ;
    
    Shimmed = Shim.ShimmedField.assessfielddistribution( voi ) ;

    Results.Std.Shimmed      = Shimmed.std ;
    Results.Std.ShimmedImprovement = 100 * (1 - Shimmed.std  ./ Original.std  )  ;
    Results.Mean.Shimmed     = Shimmed.mean ;
    Results.Median.Shimmed   = Shimmed.median ;
    Results.Norm.Shimmed     = Shimmed.norm ;
end

fprintf('\n')
display( '----- Shim Results -----' )
fprintf('\n')
display( '---- St. dev. (Hz) & Improvement (%) -----' )
display(Results.Std)
display( '-----  Mean (Hz)   -----' )
display(Results.Mean)
display( '----- Median (Hz)  -----' )
display(Results.Median)
display( '-----    Norm (Hz) -----' )
display(Results.Norm)

end
% =========================================================================
function [] = delete( Shim )
%DELETE  
% 
% DELETE( Shim )
% 
% Destructor. Calls Probe.deletecomport( ) followed by clear Shim

if ~isempty( Shim.Tracker )
    Shim.Tracker.delete();
end

clear Shim ;

end
% =========================================================================
function Params = calibraterealtimeupdates( Shim, Params )
%CALIBRATEREALTIMEUPDATES
% 
% CALIBRATEREALTIMEUPDATES asks user to select the median tracker measurement
% from the measurement logs corresponding to the inspired & expired field maps.
% From these, and the associated optimal currents for the 2 respiratory states,
% the following function calls are made:
%
% --> Shim.Opt.setcouplingcofficients()
% --> Shim.Opt.setdccurrentoffsets()
% --> Shim.Opt.setupdateoperator()
% 
%
% Usage
%
% Params = CALIBRATEREALTIMEUPDATES( Shim, Params ) 
%
%   Params.
%       .Inspired
%           .currents
%           .measurementLog
%
%       .Expired
%           .currents
%           .measurementLog
% 
% The returned Params struct has additional fields (Params.Inspired.medianP,
% and Params.Expired.medianP) corresponding to the user-selected medians (e.g.
% pressures)

if ~myisfield( Params, 'Inspired' ) || ~myisfield( Params, 'Expired' )
    error( 'See HELP ShimOpt.calibraterealtimeupdates' ) ;
end

if ~myisfield( Params.Inspired, 'medianP' ) || ~myisfield( Params.Expired, 'medianP' ) ...
    || isempty( Params.Inspired.medianP ) || isempty( Params.Expired.medianP ) 
% User selects median p-measurement

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
else

    ShimUse.display( 'Using user-supplied median p-measurements...' ) ;

    ShimUse.display( ...
        ['Inspired : ' ] ) ;
    ShimUse.display( ...
        ['Median measurement : ' num2str( Params.Inspired.medianP )] ) ;

    ShimUse.display( ...
        ['Expired : ' ] ) ;
    ShimUse.display( ...
        ['Median measurement : ' num2str( Params.Expired.medianP )] ) ;
end


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
function [] = setcouplingcoefficients( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pInspired, pExpired )
%SETCOUPLINGCOEFFICIENTS
%
% [] = SETCOUPLINGCOEFFICIENTS( Shim, iIn, iEx, pIn, pEx )
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
function [] = setdccurrentoffsets( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pInspired, pExpired )
%SETDCCURRENTOFFSETS
% 
% Compute and set optimal shim DC current offsets (bias)
%
% [] = SETDCCURRENTOFFSET( Shim, iIn, iEx, pIn, pEx  )
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

[X,Y,Z]      = Field.getvoxelpositions ;
[X0, Y0, Z0] = Shim.getvoxelpositions ;

% -------
% translate original coordinates to shifted ref. frame
dR = Field.Hdr.Private_0019_1014 ; % Private header field indicating shifted position of table w.r.t ISO
assert( dR(1) == 0, 'Table shifted in L/R direction?' ) ;
assert( dR(2) == 0, 'Table shifted in A/P direction?' ) ;

if ( dR(3) ~= 0 ) 
% field positions originally at Z0 have been shifted
    display('Correcting for table shift with respect to shim reference images')
    Shim.Hdr.ImagePositionPatient(3) = Shim.Hdr.ImagePositionPatient(3) - dR(3) ;   
end

% TODO check if voxel positions already happen to coincide
% if all( size(X) == size(X0) )
%     if ~all( X(:)==X0(:)) || ~
Shim.resliceimg( X, Y, Z ) ;

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

end
% =========================================================================
function [] = setshimmedfield( Shim, Field )
%SETSHIMMEDFIELD 
%
% [] = SETSHIMMEDFIELD( Shim, Field )
%
% Sets Shim.ShimmedField
%
% Field is a FieldEval type object with .img in Hz

Shim.ShimmedField = Field.copy() ;

end
% =========================================================================
function [] = setshimvolumeofinterest( Shim, mask )
%SETSHIMVOLUMEOFINTEREST 
% 
% [] = SETSHIMVOLUMEOFINTEREST( Shim, mask )
%
% Sets Shim.Field.Hdr.MaskingImage
%
% mask is a binary image (with the same dimensions as Shim.Field.img) of 
% the desired shim region.

assert( all( size(mask) == size( Shim.Field.img ) ), ...
    'mask (shim VOI) and target field (Shim.Field.img) must be the same size' ) ; 
Shim.Field.Hdr.MaskingImage = mask ;

end
% =========================================================================
function [] = setforwardmodelfield( Shim )
% SETFORWARDMODELFIELD
%
% [] = SETFORWARDMODELFIELD( Shim ) ;
%
% Sets Shim.Model.field --- the predicted shim field for the *current* set of currents 
% (held in Shim.Model.currents)
    
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
dbstop in ShimOpt at 532;
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
function [] = setupdateoperator( Shim )
% SETUPDATEOPERATOR
%
% [] = SETUPDATEOPERATOR( Shim ) ;
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
function [PredictedField] = predictshimmedfield( Shim, currents )
%PREDICTSHIMMEDFIELD
%
% [PredictedField] = PREDICTSHIMMEDFIELD( Shim ) ;
% 
% Returns FieldEval-type object PredictedField where
%
% PredictedField.img = ( Shim.Field.img + Shim.Model.field ) ;
%
% NOTE
%   The regions of spatial support for Shim.Model.field and Shim.Field.img 
%   are likely somewhat different (tho hopefully overlapping!).
%   PredictedField.img does not account for the finite spatial support of 
%   either field term!

% voi = logical( Shim.Field.Hdr.MaskingImage ) ;

if nargin == 2
    currentsOriginal = Shim.Model.currents ; 
    Shim.Model.currents = currents ;
else
    assert( ~isempty( Shim.Model.currents ), ...
    'Requires valid set of shim currents in Shim.Model.currents to predict shim field.' ) ;
end
    
% set Shim.Model.field according to Shim.Model.currents
Shim.setforwardmodelfield() ; 

PredictedField     = Shim.Field.copy() ;
PredictedField.img = ( Shim.Field.img + Shim.Model.field ) ;

if nargin == 2
% revert entry
Shim.Model.currents = currentsOriginal ;
end

end
% =========================================================================
function [f0, f0Voi, f0VoiShimmed] = optimizelarmor( Shim, voi )
%OPTIMIZELARMOR 
%
% [f0, f0Voi, f0VoiShimmed] = OPTIMIZELARMOR( Shim, voi ) 


voi = logical( voi ) ;

PredictedField = Shim.predictshimmedfield() ;

f0           = Shim.Field.Hdr.ImagingFrequency *(1E6) ; % [units: Hz]

f0Voi        = f0 + mean( Shim.Field.img( voi ) ) ;
f0VoiShimmed = f0 + mean( PredictedField.img( voi ) ) ;

end
% =========================================================================
function [currents] = optimizeshimcurrents( Shim, Specs, Params, FieldShift, checknonlinearconstraints )
%OPTIMIZESHIMCURRENTS 
%
% currents = OPTIMIZESHIMCURRENTS( Shim, Specs, Params )
%
% currents = OPTIMIZESHIMCURRENTS( Shim, Specs, Params, FieldShift )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]
%
%   .dataWeights
%       Array (size of Shim.Field.img) of data reliability weights [0:1].
%       Field entries corresponding to lower weighting coefficients receive
%       less consideration in the shim optimization. 
%       [default: ones(size(Shim.Field.img))]
%   
%   .Inspired.medianP
%   .Expired.medianP
%       Scalar values
%

% TODO add Scanner shim feature
%
% if ~myisfield(Params, 'controlOffsetAndLinear')
%     Params.controlOffsetAndLinear = false;
% end
%
% if Params.controlOffsetAndLinear
%     
%     % compute the shimming offset (Hz)
%     
%     offset = sum(b)/sum(M*ones(size(b)));
%     
%     %compute the desired linear shimming (Hz/mm)
%     
%     [X,Y,Z] = Shim.Field.getvoxelpositions();
%     
%     xx = M*X(:);
%     yy = M*Y(:);
%     zz = M*Z(:);
%     
%     dBdx = xx'*b./(xx'*xx);
%     dBdy = yy'*b./(yy'*yy);
%     dBdz = zz'*b./(zz'*zz);
%     
%     % Create and write on test file
%     
%     fid = fopen([Params.dataLoadDir datestr(now, 30) '-values_for_offline_shimming.txt'], 'w+');
%     fprintf(fid, '%f\n', offset);
%     fprintf(fid, '%f\n', dBdx);
%     fprintf(fid, '%f\n', dBdy);
%     fprintf(fid, '%f\n', dBdz);
%     fclose(fid);
% end

% DEFAULT_REGULARIZATIONPARAMETER     = 0;
DEFAULT_ISRETURNINGPSEUDOINVERSE    = 0;

if nargin < 2
    error('Function requires at least 2 arguments: a ShimOpt instance, and a ShimSpecs instance')
elseif nargin == 3
    Params.dummy = [];
    DEFAULT_ISSOLVINGAUGMENTEDSYSTEM    = false;
    nInputFields = 1;
elseif nargin == 4
    DEFAULT_ISSOLVINGAUGMENTEDSYSTEM    = true;
    nInptFields = 2;
end

if ~myisfield(Params, 'isReturningPseudoInverse') || isempty( Params.isReturningPseudoInverse ) 
    Params.isReturningPseudoInverse = DEFAULT_ISRETURNINGPSEUDOINVERSE ; 
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end

if ~myisfield(Params, 'isSolvingAugmentedSystem') || isempty( Params.isSolvingAugmentedSystem )
    Params.isSolvingAugmentedSystem = DEFAULT_ISSOLVINGAUGMENTEDSYSTEM ;
    
end

if Params.isSolvingAugmentedSystem
    % b0 = ( pIn* FieldShift.img(:) - pEx * Shims.Field.img(:) )/dp ;
    assert( myisfield( Params, 'Inspired') ) 
    pIn = Params.Inspired.medianP ;
    pEx = Params.Expired.medianP ;
    dp  = pIn - pEx ;
end

% if ~myisfield(Params, 'regularizationParameter') || isempty( Params.regularizationParameter ) 
%     Params.regularizationParameter = DEFAULT_REGULARIZATIONPARAMETER ;
% end


% Params for conjugate-gradient optimization
CgParams.tolerance     = 1E-10 ;
CgParams.maxIterations = 100000 ;    


nImg = numel( Shim.Field.img(:) ) ; % number of voxels



% -------
% define matrix of data-weighting coefficients : W
if ~myisfield( Params, 'dataWeights' ) || isempty( Params.dataWeights ) 

    W = speye( nImg, nImg ) ;

else

    assert( numel( Params.dataWeights ) == nImg ) 

    if ( size( Params.dataWeights, 1 ) ~= nImg ) || ( size( Params.dataWeights, 2) ~= nImg )
        
        W = spdiags( Params.dataWeights(:), 0, nImg, nImg ) ;

    end

end

M  = Shim.gettruncationoperator*W ;


A  = M*Shim.getshimoperator ; % masked current-to-field operator
MA = A;

if ~Params.isSolvingAugmentedSystem
    
    solutionVectorLength = Specs.Amp.nActiveChannels ;

    b = M*(-Shim.Field.img(:)) ;

else 
   
    % stacked/augmented solution vector length (top half: dc currents, bottom: insp-exp current shifts)
    solutionVectorLength = 2*Specs.Amp.nActiveChannels ;

    % the 'dc' field, assuming a linear model of the field shift from respiration:
    b0 = Shim.Field.img(:) ;

    % field shift from RIRO  
    db = FieldShift.img(:) ;

    % augmented data vector : dc field, field shift -- vertically concatenated
    b = [ M*-b0 ; M*-db ] ;
    
    % augmented matrix operator
    AA = [A ;  A ] ;

end


% ------- 
% Least-squares solution via conjugate gradients
x = cgls( A'*A, ... % least squares operator
          A'*b, ... % effective solution vector
          zeros( [solutionVectorLength 1] ), ... % initial model (currents) guess
          CgParams ) ;

if ~Params.isSolvingAugmentedSystem
    
    currents = x ;
    Shim.Model.currents = x ;
    
else

    [i0, di] = splitsolutionvector( x ) ;
    
    currents = i0 ;
    Shim.Model.currents = i0 ;
    Shim.Model.couplingCoefficients = di/dp ;

end


% ------- 
% TODO: CHECK SOL. VECTOR FOR CONDITIONS 
% + allow more optional params for optimization 
isCurrentSolutionOk = false ;

if ~Params.isReturningPseudoInverse && ~isCurrentSolutionOk

    Options = optimset(...
        'DerivativeCheck','off',...
        'GradObj','on',...
        'Display', 'off',... %'iter-detailed',...
        'MaxFunEvals',36000,...
        'TolX',1e-11,...
        'TolCon',1E-8);

    tic
    if ~Params.isSolvingAugmentedSystem
   
        [currents] = fmincon( ...
            @shimcost,...
            zeros( solutionVectorLength, 1),...
            [],...
            [],...
            [],...
            [],...
            -Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            @checknonlinearconstraints,...
            Options);

        Shim.Model.currents = currents ;

    else

        [currents] = fmincon( ...
            @shimcost,...
            zeros( solutionVectorLength, 1),...
            [],...
            [],...
            [],...
            [],...
            -Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            @checknonlinearconstraints_augmented,...
            Options);
        
        [currents, currentsExpired] = splitsolutionvector( currents ) ;


    end

    toc
    
end

function [f, df] = shimcost( currents )
    
     y = A*currents - b;
     f = y'*y;
     df = 2*A'*y;
    
end

function [C, Ceq] = checknonlinearconstraints_riro( solutionVector )
%CHECKNONLINEARCONSTRAINTS_RIRO
%
% Checks solution satisfies nonlinear system constraints for augmented (real-time) problem
% 
% C(x) <= 0
% (e.g. x = currents)

Ceq = [];

[i0, di] = splitsolutionvector( solutionVector ) ;

didp = di/(Params.Inspired.medianP - Params.Expired.medianP) ;

% dcCurrents ok?
[dcCurrentsC, ~] = checknonlinearconstraints( i0 ) ;
% inspired Currents ok?
[inCurrentsC, ~] = checknonlinearconstraints( i0 + ( didp * pIn ) ) ;
% expired Currents ok?
[exCurrentsC, ~] = checknonlinearconstraints( i0 + ( didp * pEx ) ) ;

C = [dcCurrentsC; inCurrentsC; exCurrentsC];

end

function [i0, di] = splitsolutionvector( solutionVector ) 
%SPLITSOLUTIONVECTOR
%
% De-concatenates vertically stacked/agumented solution vector into 
%
% [i0, di]   
% 
% i0 
%   vector of DC static shim currents
% 
% di
%   vector of current differences (inspired - expired), such that the
%   time-varying shim current i(t) 
%       = i0 (di/dp)*p(t)
%
%   where di/dp are the relational constants that scale the respiratory
%   tracking data p(t) into real-time shim current corrections

    assert( length(solutionVector) == 2*Specs.Amp.nActiveChannels ) ;
    i0 = solutionVector(1:Specs.Amp.nActiveChannels) ;
    di = solutionVector((Specs.Amp.nActiveChannels+1):end) ;

end


end
%end
% =========================================================================
% =========================================================================
methods(Access=protected)

% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true)

% =========================================================================
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

DEFAULT_SHIMMODE            = 'static' ;

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_TRACKERSPECS ;
end

if ~myisfield( Params, 'isInterpolatingReferenceMaps' ) || isempty(Params.isInterpolatingReferenceMaps)
   Params.isInterpolatingReferenceMaps = DEFAULT_ISINTERPOLATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'shimMode' ) || isempty(Params.shimMode)
   Params.shimMode = DEFAULT_SHIMMODE ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Abstract)
% =========================================================================

% =========================================================================
% =========================================================================
end

end
