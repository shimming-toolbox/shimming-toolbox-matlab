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
%       .Aux
%           [default : isempty(Shim.Aux)]
%           When Shim does not itself refer to the scanner shims, then .Aux 
%           is a ShimOpt object corresponding to the MRI host system.
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
% Updated::20170306::ryan.topfer@polymtl.ca
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
    % Model ; % Modeled quantities for shimming
    Tracker ; % object of type Tracking (e.g. ProbeTracking)
    ShimmedField; % object of type FieldEval 
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( Params, Field )
%SHIMOPT - Shim Optimization

Shim.img = [];
Shim.Hdr = [];
Shim.Model = [ ] ; 
Shim.Aux = []; 

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt.assigndefaultparameters( Params ) ;

% .......
% Load shim basis if provided 
if ~isempty(Params.pathToShimReferenceMaps)
    ShimUse.customdisplay(['\n Preparing for shim ...  \n\n'...
            'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

    load( Params.pathToShimReferenceMaps ) ;
    % Loads .mat containing Shim struct
    % which has fields
    % Shim.img              - linear dB/dI 'current-to-field' opterator
    % Shim.Hdr              - defines info like voxel locations 
    % Shim.Hdr.MaskingImage - defines spatial support of reference maps
end

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
    ShimUse.customdisplay( ['\n ------- \n ' ...
        'Determine median measurement over inspired apnea : \n \n '] ) ;

    Params.Inspired.medianP = ...
        Shim.Tracker.userselectmedianmeasurement( Params.Inspired.measurementLog ) ; 

    ShimUse.customdisplay( ...
        ['Median measurement : ' num2str( Params.Inspired.medianP )] ) ;

    ShimUse.customdisplay( ['\n ------- \n ' ...
        'Determine median measurement over expired apnea : \n \n '] ) ;

    Params.Expired.medianP = ...
        Shim.Tracker.userselectmedianmeasurement( Params.Expired.measurementLog ) ; 

    ShimUse.customdisplay( ...
        ['Median measurement : ' num2str( Params.Expired.medianP )] ) ;
else

    ShimUse.customdisplay( 'Using user-supplied median p-measurements...' ) ;

    ShimUse.customdisplay( ...
        ['Inspired : ' ] ) ;
    ShimUse.customdisplay( ...
        ['Median measurement : ' num2str( Params.Inspired.medianP )] ) ;

    ShimUse.customdisplay( ...
        ['Expired : ' ] ) ;
    ShimUse.customdisplay( ...
        ['Median measurement : ' num2str( Params.Expired.medianP )] ) ;
end


Shim.setcouplingcoefficients( ...
    Params.Inspired.currents, Params.Expired.currents, ...
    Params.Inspired.medianP, Params.Expired.medianP ) ;

Shim.setdccurrentoffsets( ...
    Params.Inspired.currents, Params.Expired.currents, ...
    Params.Inspired.medianP, Params.Expired.medianP ) ;

ShimUse.customdisplay( ['\n ------- \n ' ...
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

% Private header field indicating absolute position of table [units: mm?]
if ~myisfield( Field.Hdr, 'Private_0019_1013' ) 
    Field.Hdr.Private_0019_1013 = input( 'Enter table position in units of mm (e.g. -1638): ' ) ;
end

% --------------tablePos0-----
tablePos0 = double( Shim.Hdr.Private_0019_1013 )
% -----tablePos1--------------
tablePos1 = double( Field.Hdr.Private_0019_1013 ) 


% tablePos0 empty if Shim refers to the scanner shims (spatially fixed)
if ~isempty( tablePos0 )  
    
    % -------
    % translate original coordinates to shifted ref. frame
    warning('Adjusting for variable table position...  ')
    
    dR = tablePos1 - tablePos0 

    assert( dR(1) == 0, 'Table shifted in L/R direction?' ) ;
    assert( dR(2) == 0, 'Table shifted in A/P direction?' ) ;

    if ( dR(3) ~= 0 ) 
    % field positions originally at Z0 have been shifted
        warning('Correcting for table shift with respect to shim reference images')
        Shim.Hdr.ImagePositionPatient(3) = Shim.Hdr.ImagePositionPatient(3) - dR(3) ;   
    end

end
% -------
% check if voxel positions already happen to coincide. if they do, don't interpolate (time consuming).
if any( size(X) ~= size(X0) ) || any( X0(:) ~= X(:) ) || any( Y0(:) ~= Y(:) ) || any( Z0(:) ~= Z(:) )
    Shim.resliceimg( X, Y, Z ) ;
else
    % voxel positions already coincide,
    % i.e.
    assert( all(X0(:) == X(:) ) && all( Y0(:) == Y(:) ) && all( Z0(:) == Z(:) ) ) ;
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

if ~isempty( Shim.Aux )
    Shim.Aux.Field = Shim.Field ;
    Shim.Aux.interpolatetoimggrid( Shim.Field ) ;
end

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
function mask = getvaliditymask( Shim, Fields, Params )
%GETVALIDITYMASK 
%
% Returns binary mask - TRUE where field values are well defined and within 
% the expected range and where (interpolated) shim reference maps are well 
% defined.
%
% mask = GETVALIDITYMASK( Shim, Fields )
% mask = GETVALIDITYMASK( Shim, Fields, Params )
% 
% Fields : a cell array of FieldEval-type objects 
% (e.g. may correspond to 'Inspired' and/or 'Expired' fields.)
%
% .......................
%   
% The following Params.fields are supported
%
% .maxAbsField 
%   maximum absolute voxel value assumed to represent an accurate field
%   measurement. Voxels with abs-values greater than this might stem from
%   errors in the unwrapping.  [default: 500 Hz]
%   (To ignore this criterion, set value to Inf)

%TODO
%
% .maxFieldDifference
%   in terms of st. dev?
%
%   maximum absolute voxel-wise difference assumed to be valid between Field1 &
%   Field2 (e.g. 'inspired field' vs. 'expired field') [default: 150 Hz] (See
%   Verma T, Magn Reson Med, 2014)
%
% (Set either to Inf to ignore the criterion)

DEFAULT_MAXABSFIELD        = 500 ;
DEFAULT_MAXFIELDDIFFERENCE = 150 ;

mask = Shim.getshimsupport() ;

if nargin < 2
    warning('No input Fields provided? Returned mask is simply the shim field support.')
    return;
elseif nargin == 2 || isempty(Params)
    Params.dummy = [] ;
end

if ~myisfield( Params, 'maxAbsField' ) || isempty( Params.maxAbsField ) 
    Params.maxAbsField = DEFAULT_MAXABSFIELD ;
end

if ~myisfield( Params, 'maxFieldDifference' ) || isempty( Params.maxFieldDifference ) 
    Params.maxFieldDifference = DEFAULT_MAXFIELDDIFFERENCE ;
end


nFields = length( Fields ) ;

for iField = 1 : nFields

    mask = mask & logical(Fields{iField}.Hdr.MaskingImage) ;

    mask = mask & ( abs(Fields{iField}.img) <= Params.maxAbsField ) ;

end

% if nargin == 3 
% 
%     nFields = length( Fields ) ;
%
%     if nFields > 1 % max field deviation?
%         
%         sumField = 0;
%
%         for iField = 1 : nFields
%             
%             sumField = sumField + Fields{iField}.img ;
%
%         end
%
%         meanField = sumField/nFields ;
%
%     end
%
%     for iField = 1 : nFields
%
%         mask = mask & logical(Field.Hdr.MaskingImage) ;
%
%         mask = mask & ( abs(Field{iField}.img) <= Params.maxAbsField ) ;
%
%         if iField > 1
%             
%             mask = mask & ( abs( Fields{iField}.img - meanField ) <= Params.maxFieldDifference ) ;
%
%         end
%     end
%
% end

if nnz(mask(:)) == 0
    warning( 'Calculated region of validity is nul. Returning mask array of zeros.' ) ;
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

currents = Shim.Model.currents + ...
        Shim.Model.couplingCoefficients * Shim.Tracker.Data.p(end) ; 

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
function [currents] = optimizeshimcurrents( Shim, Specs, Params, checknonlinearconstraints )
%OPTIMIZESHIMCURRENTS 
%
% currents = OPTIMIZESHIMCURRENTS( Shim, Specs, Params, @checknonlinearconstraints )
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
    assert( myisfield( Shim.Field.Model, 'Shift') && ~isempty( Shim.Field.Model.Shift.img ) ) 

    % change to Params.minP and Params.maxP
    pIn = Params.pMax 
    pEx = Params.pMin 
    dp  = Shim.Field.Model.Shift.Aux.Tracker.Data.p

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
    db = Shim.Field.Model.Shift.img(:) ;

    % augmented data vector : dc field, field shift -- vertically concatenated
    b = [ M*-b0 ; M*-db ] ;
    
    % augmented matrix operator
    A = [MA zeros(size(A)); zeros(size(A)) MA] ;

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
isCurrentSolutionOk = all( checknonlinearconstraints( x ) < 0 ) ;

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
            x,...
            [],...
            [],...
            [],...
            [],...
            -Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            checknonlinearconstraints,...
            Options);

        Shim.Model.currents = currents ;

    else

        [currents] = fmincon( ...
            @shimcost,...
            x,...
            [],...
            [],...
            [],...
            [],...
            -Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            Params.maxCurrentPerChannel * ones(solutionVectorLength,1),...
            @checknonlinearconstraints_riro,...
            Options);
        
        [currents, di] = splitsolutionvector( currents ) ;
        
        Shim.Model.currents = currents ;
        Shim.Model.couplingCoefficients = di/dp ;


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

didp = di/dp ;

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
function  [ dataWeights ] = derivedataweights( Mag, targetEchoTime, targetMask )
%DERIVEDATAWEIGHTS  
% 
% dataWeights = DERIVEDATAWEIGHTS( Mag, targetEchoTime )
% dataWeights = DERIVEDATAWEIGHTS( Mag, targetEchoTime, targetMask )
% 
% Assumes basic model of T2* signal decay using the magnitude images of the
% 2 echoes in the dual-echo GRE
%
% if: Mag(TE) = M0 * exp(-TE/T2star)
% then:
% t2star  = -deltaTe ./ log( Mag2.img ./ Mag1.img ) ; 
% 
% dataWeights is the predicted + normalized signal intensity at the targetEchoTime
% (e.g. a long TE for a GRE-EPI sequence).
%
% targetMask, if provided, is a logical array specifying a priority region
% (e.g. spinal canal) that will receive maximal (unity) weighting in dataWeights,
% irrespective of the t2*-forecasted signal.
% (i.e. dataWeights( targetMask ) = 1 ;)


[t2star, reliabilityMask] = computet2star( Mag ) ;
% if: Mag(TE) = M0 * exp(-TE/T2star)
M0 = Mag{1}.img ./ exp(-Mag{1}.Hdr.EchoTime/t2star) ;
M0( ~reliabilityMask ) = 0 ; 
M0 = M0 ./ max( M0(:) ) ;

if ( nargin < 2 ) || isempty( targetEchoTime ) 
    targetEchoTime = Mag{2}.Hdr.EchoTime ;
end

magAtEcho = M0 .* exp( -targetEchoTime/t2star ) ; % forecasted undistorted magnitude at EchoTime 
magAtEcho = magAtEcho./max( magAtEcho(:) ) ;

dataWeights = magAtEcho ;

if ( nargin == 3 ) & ~isempty( targetMask )
    assert( all( size( targetMask ) == Mag{1}.getgridsize ) ) 
    dataWeights( targetMask ) = 1 ;
end

function [ t2star, reliabilityMask ] = computet2star( Mag )
%COMPUTET2STAR
% Returns poor man's estimate of t2star [units: ms] using Mag images from 2 echoes
% 
% [ t2star, reliabilityMask ] = computet2star( Mag )
% 
% if: Mag(TE) = M0 * exp(-TE/T2star)
% then: 
% t2star  = -deltaTe ./ log( Mag{2}.img ./ Mag{1}.img ) ; 
% 
% reliabilityMask = t2star > 0 ;
    deltaTe = (Mag{2}.Hdr.EchoTime - Mag{1}.Hdr.EchoTime) ; % [units: ms]
    
    t2star  = -deltaTe ./ log( Mag{2}.img ./ Mag{1}.img ) ; 

    t2star(t2star == -Inf) = Inf ; % corresponding to saturated voxels ;
    reliabilityMask = t2star > 0 ;
    t2star( ~reliabilityMask ) = 0 ;% unreliable voxels (i.e. where Mag(TE2)>Mag(TE1) )

end

end
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

DEFAULT_INSTITUTIONNAME = 'IUGM' ;
DEFAULT_STATIONNAME     = 'MRC35049' ;

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

if ~myisfield( Params, 'InstitutionName' ) || isempty(Params.InstitutionName)
   Params.InstitutionName = DEFAULT_INSTITUTIONNAME ;
end

if ~myisfield( Params, 'StationName' ) || isempty(Params.StationName)
   Params.StationName = DEFAULT_STATIONNAME ;
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
