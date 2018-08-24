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
%           Object of type FieldEval pertaining to field distribution to be shimmed
%
%       .Aux
%           .Shim
%               When Shim does not itself refer to the scanner shims, then .Aux 
%               is a ShimOpt object corresponding to the MRI host system.
%
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
%       .System
%           .currents
%               
%           .Specs
%               Object of a type of ShimSpecs sub-class (e.g. ShimSpecs_IUGM_Prisma_fit) 
%
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimUse
%
% ShimOpt is a FieldEval subclass [ShimOpt < FieldEval < MaRdI]
%     
% =========================================================================
% Updated::20180723::ryan.topfer@polymtl.ca
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
    ShimmedField; % object of type FieldEval 
    Tracker;
    System; %
end
properties( Hidden = true )
    Ref ; % original shim reference maps before interpolation
end
% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( Params, Field )
%SHIMOPT - Shim Optimization

Shim.img    = [] ;
Shim.Hdr    = [] ;
Shim.Field  = [] ;  
Shim.Model  = [] ;
Shim.Aux    = [] ;
Shim.System = [] ;

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOpt.assigndefaultparameters( Params ) ;

% .......
% Load shim basis if provided 
if ~isempty(Params.pathToShimReferenceMaps)
    
ShimUse.customdisplay(['\n Preparing for shim ...  \n\n'...
        'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

    % Loads .mat: struct with fields
    %
    % .img 
    %       linear dB/dI 'current-to-field' operator
    % .Hdr
    %       dicom Hdr from the reference map acquisition 
    RefMaps = load( Params.pathToShimReferenceMaps ) ; 

    %%-----
    % dB/dI linear 'Current-to-Field' operator
    Shim.img              = RefMaps.Shim.img ;
    Shim.Hdr              = RefMaps.Shim.Hdr ;

    Shim.Ref.img = Shim.img ;
    Shim.Ref.Hdr = Shim.Hdr ;

end

Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 

if (nargin == 2) && (~isempty(Field))
    
    Shim.setoriginalfield( Field ) ;

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

% function [ShimCopy] = copy( Shim )
% %COPY  
% % 
% % ShimCopy = COPY( Shim )
% ShimCopy = ShimOpt()
% ShimCopy.img    = Shim.img ;
% ShimCopy.Hdr    = Shim.Hdr ;
% ShimCopy.Field  = Shim.Field.copy() ;
%
% ShimCopy.Model  = Shim.Model ;
% ShimCopy.System = Shim.System ;
%
% ShimCopy.Tracker = Shim.Tracker.copyinert() ;
% ShimCopy.Aux     = Shim.Aux.copy() ;
%
% end
% =========================================================================
function [] = resettoreference( Shim )
%RESETTOREFERENCE  
% 
% RESETTOREFERENCE( Shim )
% 
% Returns Shim.img and Shim.Hdr to their original (reference) values (before
% interpolations, cropping etc.) 

Shim.img = Shim.Ref.img ;
Shim.Hdr = Shim.Ref.Hdr ;

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

% A = Shim.getshimoperator ;
%
% Shim.Model.couplingCoefficients = ...
%     A*(currentsInspired - currentsExpired)/(pInspired - pExpired) ;

Shim.Model.couplingCoefficients = ...
     (currentsInspired - currentsExpired)/(pInspired - pExpired) ;

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
    % TODO
    % none of this table pos. stuff matters for the ac/dc array and is all deprecated anyway: use siemens shadow hdr fields for table pos.
    % Tmp fix, merely copy in case i've planted an error flag somewhere else for an empty entry:
    Field.Hdr.Private_0019_1013 = Shim.Hdr.Private_0019_1013 ;
end

% --------------tablePos0-----
tablePos0 = double( Shim.Hdr.Private_0019_1013 ) ;

% -----tablePos1--------------
tablePos1 = double( Field.Hdr.Private_0019_1013 ) ;


% tablePos0 empty if Shim refers to the scanner shims (spatially fixed)
% if ~isempty( tablePos0 )  
%     
%     % -------
%     % translate original coordinates to shifted ref. frame
%     warning('Adjusting for variable table position...  ')
%     
%     dR = tablePos1 - tablePos0 
%
%     assert( dR(1) == 0, 'Table shifted in L/R direction?' ) ;
%     assert( dR(2) == 0, 'Table shifted in A/P direction?' ) ;
%
%     if ( dR(3) ~= 0 ) 
%     % field positions originally at Z0 have been shifted
%     % tablePosition is increasinly negative the more it is into the scanner.
%     % the opposite is true for the z-coordinate of a voxel in the dicom reference system.
%         warning('Correcting for table shift with respect to shim reference images')
%     % dbstop in ShimOpt at 455
%         % Shim.Hdr.ImagePositionPatient(3) = Shim.Hdr.ImagePositionPatient(3) + dR(3) ;   
%     end
%
% end
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
function [] = setoriginalfield( Shim, Field, currents )
%SETORIGINALFIELD 
%
% [] = SETORIGINALFIELD( Shim, Field )
% [] = SETORIGINALFIELD( Shim, Field, currents )
%
% Sets Shim.Field
%
% Field is a FieldEval type object with .img in Hz

if nargin < 2
    error('Not enough input arguments.') ;
elseif nargin == 2
    currents = 0;
    warning('Assuming field map was acquired with all shim channels at 0 A.');
end

Shim.Model.currents = currents ;
Shim.Field = Field.copy() ;

Shim.interpolatetoimggrid( Shim.Field ) ;
Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

if ~isempty( Shim.Aux )  
    Shim.Aux.setoriginalfield( Shim.Field ) ;
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
function [] = setshimvolumeofinterestriro( Shim, mask )
%SETSHIMVOLUMEOFINTERESTRIRO 
% 
% [] = SETSHIMVOLUMEOFINTERESTRIRO( Shim, mask )
%
% Sets Shim.Field.Model.Riro.Hdr.MaskingImage
%
% mask is a binary image (with the same dimensions as Shim.Field.Model.Riro.img) of 
% the desired shim region.

assert( myisfield( Shim.Field.Model, 'Riro' ) && ~isempty( Shim.Field.Model.Riro ) )

assert( all( size(mask) == size( Shim.Field.Model.Riro.img ) ), ...
    'mask (shim VOI) and target field (Shim.Field.Model.Riro.img) must be the same size' ) ; 

Shim.Field.Model.Riro.Hdr.MaskingImage = mask ;

end
% =========================================================================
function [shimCorrection] = forwardmodelshimcorrection( Shim, correctionCoefficients )
%FORWARDMODELSHIMCORRECTION
% 
% shimCorrection = FORWARDMODELSHIMCORRECTION( Shim, correctionCoefficients ) ;
% 
% Forward projection of the shim correction :
%
% shimCorrection = reshape( Shim.getshimoperator()*correctionCoefficients, Shim.Field.getgridsize() ) ;

assert( nargin == 2 )

shimCorrection = reshape( Shim.getshimoperator()*correctionCoefficients, Shim.Field.getgridsize() ) ;

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
c = Shim.forwardmodelshimcorrection( Shim.Model.couplingCoefficients ) ;
c = c(:) ;

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
% Returns sparse linear truncation operator M 
% i.e. M*b, 'picks out' the VOI voxels from vector b 
% where the VOI is defined by the full 3d array Shim.Field.Hdr.MaskingImage
%
% if b has length nImg, and nVoi is the # of non-zero VOI voxels, then length(M*b)=nVoi.

nVoxelsImg = numel( Shim.Field.Hdr.MaskingImage ) ;
nVoxelsVoi = nnz( Shim.Field.Hdr.MaskingImage ) ;

indicesVoi = find( Shim.Field.Hdr.MaskingImage(:) ) ;

M = sparse( [1:nVoxelsVoi], indicesVoi, ones([nVoxelsVoi 1]), nVoxelsVoi, nVoxelsImg ) ;

end
% =========================================================================
function M = gettruncationoperatorriro( Shim )
% GETTRUNCATIONOPERATORRIRO
%
% M = GETTRUNCATIONOPERATORRIRO( Shim ) ;
%
% Returns sparse linear truncation operator M 
% i.e. M*b, 'picks out' the VOI voxels from vector b 
% where the VOI is defined by the full 3d array Shim.Field.Model.Riro.Hdr.MaskingImage
%
% if b has length nImg, and nVoi is the # of non-zero VOI voxels, then length(M*b)=nVoi.

nVoxelsImg = numel( Shim.Field.Model.Riro.Hdr.MaskingImage ) ;
nVoxelsVoi = nnz( Shim.Field.Model.Riro.Hdr.MaskingImage ) ;

indicesVoi = find( Shim.Field.Model.Riro.Hdr.MaskingImage(:) ) ;

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
function currents = computerealtimeupdate( Shim, p )
% COMPUTEREALTIMEUPDATE
% 
% Usage
%
% currents = COMPUTEREALTIMEUPDATE( Shim, p )
%
%   currents = Shim.Model.currents + p*Shim.Model.couplingCoefficients ;

currents = Shim.Model.currents + p*Shim.Model.couplingCoefficients ;

end
% % =========================================================================
% function [ShimmedField] = GETINSTANTANEOUSFIELDSHIMMED( Shim, p )
% %GETINSTANTANEOUSFIELDSHIMMED
%
% Field = Shim.Field.getinstantaneousfield( p ) ;
%
% ShimmedField = Field.copy ;
%
% ShimmedField.img = ShimmedField.img + Shim.forwardmodelshimcorrection( ...
%
% end
% =========================================================================
function [PredictedField] = predictshimmedfield( Shim )
%PREDICTSHIMMEDFIELD
%
% [PredictedField] = PREDICTSHIMMEDFIELD( Shim ) ;
% 
% Returns FieldEval-type object(s) 
%
% PredictedField.img = ( Shim.Field.img + Shim.Model.field ) ;
%
% NOTE
%   The regions of spatial support for Shim.Model.field and Shim.Field.img 
%   are likely somewhat different (though ideally overlapping!).
%   The predictions do not account for the finite spatial support of 
%   either field term!

% voi = logical( Shim.Field.Hdr.MaskingImage ) ;

assert( ~isempty( Shim.Model.currents ), ...
    'Requires valid set of shim currents in Shim.Model.currents to predict shim field.' ) ;

PredictedField     = Shim.Field.copy() ;
shimCurrentsUpdate = Shim.Model.currents - Shim.System.currents ;
Shim.Model.field   = Shim.forwardmodelshimcorrection( shimCurrentsUpdate ) ;

PredictedField.img = PredictedField.img + Shim.Model.field ;

if ~isempty(Shim.Model.Tx.imagingFrequency)
    
    % frequency shift
    df0 = Shim.Model.Tx.imagingFrequency - (Shim.Field.Hdr.ImagingFrequency*1E6) ; % [units: Hz]
    
    if ( abs( df0 ) > 0.1 )
        PredictedField.Hdr.ImagingFrequency = Shim.Model.Tx.imagingFrequency*1E-6 ;  % [units: MHz]
        PredictedField.img = PredictedField.img + df0 ;
    end

end

if ~isempty(Shim.Aux.Model.currents)

    auxCurrentsUpdate = Shim.Aux.Model.currents - Shim.Aux.System.currents ;

    Shim.Aux.Model.field = Shim.Aux.forwardmodelshimcorrection( auxCurrentsUpdate ) ;

    PredictedField.img = PredictedField.img + Shim.Aux.Model.field ;

end

end
% =========================================================================
function [ PredictedRiro ] = predictshimmedriro( Shim, p )
%PREDICTSHIMMEDRIRO
%
% [PredictedRiro] = PREDICTSHIMMEDRIRO( Shim ) ;
% 
% Returns FieldEval-type object(s) 
%
% PredictedRiro.img  = ( Shim.Field.Model.Riro.img + Shim.Model.riro ) ;
%
% NOTE
%   The regions of spatial support for Shim.Model.field and Shim.Field.img 
%   are likely somewhat different (though ideally overlapping!).
%   The predictions do not account for the finite spatial support of 
%   either field term!

% voi = logical( Shim.Field.Hdr.MaskingImage ) ;
PredictedRiro = [] ;

assert( ~isempty( Shim.Field.Model.Riro ) && ~isempty( Shim.Field.Model.Riro.img ), ...
    'Requires input model of respiration-induced reference offset (RIRO) in Shim.Field.Model.Riro' ) ;


if ~myisfield( Shim.Model, 'couplingCoefficients' ) || isempty( Shim.Model.couplingCoefficients ) 
    
    warning('No RIRO correction: Returning input model RIRO')

else    
    
    % all field + correction terms scaled by dp :
    if nargin < 2
        PredictedRiro = Shim.Field.getriro( ) ;
        % (dp: recorded inspired - expired pressure difference)
        dp = Shim.Field.Model.Riro.Aux.Tracker.Data.p ; 
    else
        PredictedRiro = Shim.Field.getriro( p ) ;
        dp = Shim.Field.Aux.Tracker.debias( p ) ; 
    end
    
    Shim.Model.riro   = dp*Shim.forwardmodelshimcorrection( Shim.Model.couplingCoefficients ) ;
    PredictedRiro.img = PredictedRiro.img + Shim.Model.riro ;

    if ~isempty(Shim.Model.Tx.couplingCoefficients) && ( Shim.Model.Tx.couplingCoefficients ~= 0 )
       
       warning('Untested: TX RIRO adjustment') 

       PredictedRiro.img = PredictedRiro.img + dp*Shim.Model.Tx.couplingCoefficients ;

    end

    if myisfield( Shim.Aux.Model, 'couplingCoefficients' ) && ~isempty( Shim.Aux.Model.couplingCoefficients )

        Shim.Aux.Model.riro = dp*Shim.Aux.forwardmodelshimcorrection( Shim.Aux.Model.couplingCoefficients ) ;
        PredictedRiro.img   = PredictedRiro.img + Shim.Aux.Model.riro ;

    end

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
function [Corrections] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Corrections = OPTIMIZESHIMCURRENTS( Shim, Params )
%
% Corrections.static
%
% Corrections.realtime

% TODO: write usage
% Params can have the following fields 
%
%   .isOptimizingStaticTxFrequency
%       [default: TRUE]
%
%   .isOptimizingDynamicTxFrequency
%       [default: FALSE]
%   
%   .maxCurrentPerChannel
%       [default defined as the .Amp.maxCurrentPerChannel property in the ShimSpecs_ [sub]class definition]
%
%   .dataWeights
%       Array (size of Shim.Field.img) of data reliability weights [0:1].
%       Field entries corresponding to lower weighting coefficients receive
%       less consideration in the shim optimization. 
%       [default: ones(size(Shim.Field.img))]

% DEFAULT_REGULARIZATIONPARAMETER     = 0; % deprecated
DEFAULT_ISRETURNINGPSEUDOINVERSE       = false ;
DEFAULT_ISOPTIMIZINGSTATICTXFREQUENCY  = true ;
DEFAULT_ISOPTIMIZINGDYNAMICTXFREQUENCY = false ;
DEFAULT_ISOPTIMIZINGAUX                = false ;
DEFAULT_ISREALTIMESHIMMING             = false ;

DEFAULT_MEDIASAVEDIR                   = [ './shimopt_results_' datestr(now,30) '/' ] ;
DEFAULT_ISDISPLAYINGRESULTS            = true ;
DEFAULT_ISSAVINGRESULTSTABLE           = true ;

% TODO: Create new class to handle all Tx related aspects
DEFAULT_MINTXFREQUENCY = 123100100 ; % [units: Hz]
DEFAULT_MAXTXFREQUENCY = 123265000 ; % [units: Hz]


if nargin < 1
    error('Function requires at least 1 argument: a ShimOpt instance')
elseif nargin == 1 || isempty( Params ) 
    Params.dummy = [] ; % Using default parameters
end

if ~myisfield(Params, 'isReturningPseudoInverse') || isempty( Params.isReturningPseudoInverse ) 
    Params.isReturningPseudoInverse = DEFAULT_ISRETURNINGPSEUDOINVERSE ; 
end



if ~myisfield(Params, 'isOptimizingStaticTxFrequency') || isempty( Params.isOptimizingStaticTxFrequency )
    Params.isOptimizingStaticTxFrequency = DEFAULT_ISOPTIMIZINGSTATICTXFREQUENCY ;
end

if ~myisfield(Params, 'isOptimizingDynamicTxFrequency') || isempty( Params.isOptimizingDynamicTxFrequency )
    Params.isOptimizingDynamicTxFrequency = DEFAULT_ISOPTIMIZINGDYNAMICTXFREQUENCY ;
end



if ~myisfield(Params, 'maxPositiveTxFrequencyShift') || isempty( Params.maxPositiveTxFrequencyShift )
    % max + freq. shift relative to original Larmor
    Params.maxPositiveTxFrequencyShift = DEFAULT_MAXTXFREQUENCY - double(Shim.Field.Hdr.ImagingFrequency)*1E6 ;
end

if ~myisfield(Params, 'maxNegativeTxFrequencyShift') || isempty( Params.maxNegativeTxFrequencyShift )
    % max - freq. shift relative to original Larmor
    Params.maxNegativeTxFrequencyShift = DEFAULT_MINTXFREQUENCY - double(Shim.Field.Hdr.ImagingFrequency)*1E6 ;
end



assert( ~isempty( Shim.Aux ) )

% total active channels across systems
Params.nActiveChannels = 1 + Shim.System.Specs.Amp.nActiveChannels + Shim.Aux.System.Specs.Amp.nActiveChannels ; % +1 for freq. adjust

if ~myisfield(Params, 'activeStaticChannelsMask') || isempty( Params.activeStaticChannelsMask )
    Params.activeStaticChannelsMask = [ Params.isOptimizingStaticTxFrequency ; ...
                                        Shim.System.Specs.Amp.staticChannels(:) ; ...
                                        Shim.Aux.System.Specs.Amp.staticChannels(:) ; ] ;
end
    
if ~myisfield(Params, 'activeDynamicChannelsMask') || isempty( Params.activeDynamicChannelsMask )
    Params.activeDynamicChannelsMask = [ Params.isOptimizingDynamicTxFrequency ; ...
                                         Shim.System.Specs.Amp.dynamicChannels(:) ; ...
                                         Shim.Aux.System.Specs.Amp.dynamicChannels(:) ; ] ;
end


if ~myisfield(Params, 'maxCorrectionPerChannel') || isempty( Params.maxCorrectionPerChannel ) 
   
    % 'tmp' because this is a TODO: 
    % since the optimization solves for the Tx frequency *shift*, it can shift
    % up or down relative to Larmor.  For now, I'll limit the max shift equally
    % in both up/down directions, but in the future both should be specifically
    % accounted for (only comes into play in checknonlinearconstraints_static() )
    tmpMaxFreqShift = min(abs( [ Params.maxPositiveTxFrequencyShift Params.maxNegativeTxFrequencyShift ] )) ;

    Params.maxCorrectionPerChannel = [ tmpMaxFreqShift ; Shim.System.Specs.Amp.maxCurrentPerChannel ; Shim.Aux.System.Specs.Amp.maxCurrentPerChannel ]; 

end

if ~myisfield(Params, 'minCorrectionPerChannel') || isempty( Params.minCorrectionPerChannel ) 
    Params.minCorrectionPerChannel = -Params.maxCorrectionPerChannel ; 
end

assert( ( length( Params.maxCorrectionPerChannel ) == length( Params.minCorrectionPerChannel ) ) && ...
        ( length( Params.maxCorrectionPerChannel ) == Params.nActiveChannels ), ...
    'Shim system limits (Params.maxCorrectionPerChannel and Params.minCorrectionPerChannel) must possess an entry for each shim channel (primary and, if in use, auxiliary shims).' ) ;

% set min/max static correction to 0 for inactive static channels
Params.minStaticCorrectionPerChannel = Params.activeStaticChannelsMask .* Params.minCorrectionPerChannel ; 
Params.maxStaticCorrectionPerChannel = Params.activeStaticChannelsMask .* Params.maxCorrectionPerChannel ;

if ~myisfield(Params, 'assessmentVoi') || isempty( Params.assessmentVoi )
    Params.assessmentVoi = Shim.Field.Hdr.MaskingImage ;
end

if ~myisfield(Params, 'isRealtimeShimming') || isempty( Params.isRealtimeShimming )
    Params.isRealtimeShimming = DEFAULT_ISREALTIMESHIMMING ;
end

% if ~myisfield(Params, 'regularizationParameter') || isempty( Params.regularizationParameter ) 
%     Params.regularizationParameter = DEFAULT_REGULARIZATIONPARAMETER ;
% end

if ~myisfield(Params, 'mediaSaveDir') || isempty( Params.mediaSaveDir ) 
    Params.mediaSaveDir = DEFAULT_MEDIASAVEDIR ; 
end
mkdir( Params.mediaSaveDir ) ;

if ~myisfield(Params, 'isDisplayingResults') || isempty( Params.isDisplayingResults ) 
    Params.isDisplayingResults = DEFAULT_ISDISPLAYINGRESULTS ; 
end

if ~myisfield(Params, 'isSavingResultsTable') || isempty( Params.isSavingResultsTable ) 
    Params.isSavingResultsTable = DEFAULT_ISSAVINGRESULTSTABLE ; 
end

if Params.isRealtimeShimming

    assert( myisfield( Shim.Field.Model, 'Riro') && ~isempty( Shim.Field.Model.Riro.img ) ) 
    
    pMax = Params.pMax ; % max relevant/recorded pressure (e.g. ~inspired state)
    pMin = Params.pMin ; % min relevant/recorded pressure (e.g. ~expired state) 
    dp   = Shim.Field.Model.Riro.Aux.Tracker.Data.p ; % delta pressure shift to which input Riro has been scaled
    % i.e. Riro = dp * ( Inspired_Field - Expired_Field )/( inspired_pressure - expired_pressure )
    
end



nImg = numel( Shim.Field.img(:) ) ; % number of voxels

% -------
% define matrix of data-weighting coefficients for the static optimization: W0
if ~myisfield( Params, 'dataWeights' ) || isempty( Params.dataWeights ) 

    W0 = speye( nImg, nImg ) ;

else

    assert( numel( Params.dataWeights ) == nImg ) 

    if ( size( Params.dataWeights, 1 ) ~= nImg ) || ( size( Params.dataWeights, 2) ~= nImg )
        
        W0 = spdiags( Params.dataWeights(:), 0, nImg, nImg ) ;

    end

end

% truncated (VOI-masked) weighting operator for static shim: MW0
MW0 = Shim.gettruncationoperator()*W0 ; 

if Params.isRealtimeShimming 
    % define matrix of data-weighting coefficients for real-time (RIRO) correction: W1
    if ~myisfield( Params, 'dataWeightsRiro' ) || isempty( Params.dataWeightsRiro ) 

        W1 = W0 ;

    else

        assert( numel( Params.dataWeightsRiro ) == nImg ) 

        if ( size( Params.dataWeightsRiro, 1 ) ~= nImg ) || ( size( Params.dataWeightsRiro, 2) ~= nImg )
            
            W1 = spdiags( Params.dataWeightsRiro(:), 0, nImg, nImg ) ;

        end

    end

    % truncated (VOI-masked) weighting operator for static shim: MW0
    MW1 = Shim.gettruncationoperatorriro()*W1 ; 
end

% -------
% define off-resonance correction operator : A

% global resonance offset operator (Tx frequency) : A_tx
A_tx = ones( nImg, 1 ) ;

% primary shim (i.e. current-to-field) operator ('mc' since it's likely a multi-coil array) : A_mc
A_mc = Shim.getshimoperator() ; 

% auxiliary shim operator (e.g. scanner's spherical harmonic coils)
A_aux = Shim.Aux.getshimoperator() ;

% combined: static off-resonance correction operator: A0
A0 = [ MW0*A_tx MW0*A_mc MW0*A_aux ] ;


if Params.isRealtimeShimming
    
    activeChannelsMask = [ Params.activeStaticChannelsMask ; Params.activeDynamicChannelsMask ] ;
    
    A1 = [ MW1*A_tx MW1*A_mc MW1*A_aux ] ;
    % left half applies to the static field, right half to the resp.-induced dynamic component 
    A  = [ A0 zeros(size(A0)) ; zeros(size(A1)) A1 ] ; 

    Params.minDynamicCorrectionPerChannel = zeros(size(Params.activeDynamicChannelsMask)) ;
    Params.minDynamicCorrectionPerChannel( Params.activeDynamicChannelsMask ) = -Inf ;
    
    Params.maxDynamicCorrectionPerChannel = zeros(size(Params.activeDynamicChannelsMask)) ;
    Params.maxDynamicCorrectionPerChannel( Params.activeDynamicChannelsMask ) = Inf ;

    Params.minCorrectionPerChannel = [ Params.minStaticCorrectionPerChannel Params.minDynamicCorrectionPerChannel ] ;
    Params.maxCorrectionPerChannel = [ Params.maxStaticCorrectionPerChannel Params.maxDynamicCorrectionPerChannel ] ;

else
    
    A = A0 ;
    
    activeChannelsMask = [ Params.activeStaticChannelsMask ] ;

    Params.minCorrectionPerChannel = [ Params.minStaticCorrectionPerChannel ] ;
    Params.maxCorrectionPerChannel = [ Params.maxStaticCorrectionPerChannel ] ;

end

% possible TODO ?
%   Remove inactive channels from the operator (and solution vector)
%   rather than setting them to zero. e.g. A = A(:, activeChannelsMask) ; 
A(:, ~activeChannelsMask) = 0;

solutionVectorLength = size( A, 2 ) ;






% -------
% Define the solution (data) vector: b

% Solve for residual field offset (bx) *without* existing shim fields (bs) 
% *from shim channels that will be optimized*
% (i.e. if shim settings are to remain constant anyway, there is no need to
% subtract out their contribution from the field)
% 
% this way, optimizeshimcurrents() solves for the absolute shim currents,
% rather than a shift, which simplifies handling of the min/max constraints on
% the correction terms.

% retain original corrections for 'initial guess' of conjugate gradient solver
% 1st element: initial guess of zero frequency shift relative to the imaging frequency of the input field map
x0 = [ 0 ; Shim.System.currents ; Shim.Aux.System.currents ] ;

x0( ~Params.activeStaticChannelsMask ) = 0 ;

% Residual static off-resonance field without existing shim fields from the
% active static channels: bx0
%
% i.e. the respiration-independent static field, assuming a linear model of the
% time variation (van Gelderen et. al, MRM 2007):
bx0 = Shim.Field.img(:) - [A_tx A_mc A_aux]*x0 ;

if ~Params.isRealtimeShimming

    b  = MW0*-bx0 ;

else 

    % field shift from RIRO  
    db = Shim.Field.Model.Riro.img(:) ;

    %  stacked data vector : dc field, RIRO field shift -- vertically concatenated
    b = [ MW0*-bx0 ; MW1*-db ] ;

end

if Params.isRealtimeShimming 
% initial guess of zero shift for the all respiratory/dynamic terms
    x0 = [ x0 ; zeros( size(x0) ) ] ; 
end
    


% ------- 
% linear optimization 

% Params for conjugate-gradient optimization
CgParams.tolerance     = 1E-10 ;
CgParams.maxIterations = 100000 ;    

x = cgls( A'*A, ... % least squares operator
          A'*b, ... % effective solution vector
          x0, ... % initial 'guess' solution vector
          CgParams ) ;

% check linear solution
isCurrentSolutionOk = all( checknonlinearconstraints_default( x ) <= 0 ) ;

% ------- 
% nonlinear optimization
if ~Params.isReturningPseudoInverse && ~isCurrentSolutionOk

    Options = optimset(...
        'DerivativeCheck','off',...
        'GradObj','on',...
        'Display', 'off',... %'iter-detailed',...
        'MaxFunEvals',36000,...
        'TolX',1e-11,...
        'TolCon',1E-8);

    tic
   
    x = fmincon( ...
            @shimcost,...
            x,...
            [],...
            [],...
            [],...
            [],...
            Params.minCorrectionPerChannel,...
            Params.maxCorrectionPerChannel,...
            @checknonlinearconstraints_default,...
            Options);

    toc
    
end

[ fo0, i0, i0Aux, dfo, di, diAux] = splitsolutionvector( x ) ;


Shim.Model.Tx.imagingFrequency      = fo0 + (Shim.Field.Hdr.ImagingFrequency*1E6) ; % DC avg. [units: Hz]

Shim.Model.currents                 = i0 ;
% replace inactive terms with their initial values: TODO reorganize
Params.activeStaticChannelsMaskMc   = Params.activeStaticChannelsMask(2:2+length(i0)) ;
Shim.Model.currents( ~Params.activeStaticChannelsMaskMc ) = Shim.System.currents( ~Params.activeStaticChannelsMaskMc ) ;

Shim.Aux.Model.currents             = i0Aux ;
% replace inactive terms with their initial values: TODO ditto
Params.activeStaticChannelsMaskAux  = Params.activeStaticChannelsMask(2+length(i0):end) ;
Shim.Aux.Model.currents( ~Params.activeStaticChannelsMaskAux ) = Shim.Aux.System.currents( ~Params.activeStaticChannelsMaskAux ) ;

Corrections.static = [ Shim.Model.Tx.imagingFrequency; Shim.Model.currents; Shim.Aux.Model.currents ] ;

if Params.isRealtimeShimming
    Shim.Model.Tx.couplingCoefficients  = dfo/dp ; % resp-induced shift [units: Hz/unit-pressure]
    Shim.Model.couplingCoefficients     = di/dp ;
    Shim.Aux.Model.couplingCoefficients = diAux/dp ;

    Corrections.realtime = [Shim.Model.Tx.couplingCoefficients; Shim.Model.couplingCoefficients; Shim.Aux.Model.couplingCoefficients] ; 
end

% -------
% results summary 
%   TODO: clean up

T = table ;

fprintf(['\n' '-----\t' 'Optimization Results' '\t-----' '\n\n']) ;

% [table units: mA]
i0 = i0*1000 ;     

if Params.isRealtimeShimming
    di = di*1000 ;
else
    dfo = 0 ;
    di  = zeros( size( i0 ) ) ;
    diAux = zeros( size( i0Aux ) ) ;
end

x0 = [ fo0; i0; i0Aux; ] ;
dx = [ dfo; di; diAux; ] ;

Correction_Term = {'Tx Freq. (Hz)'} ;

for iCh = 1 : length(i0)
    Correction_Term(end+1) = { [ Shim.System.Specs.Id.systemName ' ' Shim.System.Specs.Id.channelNames{iCh} ' (mA)'] } ;
end

for iCh = 1 : length(i0Aux)  
    Correction_Term(end+1) = { [ Shim.Aux.System.Specs.Id.systemName ' ' Shim.Aux.System.Specs.Id.channelNames{iCh} ' (D.U.)'] } ;
end

Original = [ round( Shim.Field.Hdr.ImagingFrequency*1E6, 9, 'significant' ) ; ... 
             round( 1000*Shim.System.currents ) ; ... % [units: mA] 
             round(Shim.Aux.System.currents, 4, 'significant') ] ; 

Optimal = [ round( Shim.Model.Tx.imagingFrequency, 9, 'significant' ) ; ...
            round(i0) ; ...
            round(i0Aux, 4, 'significant') ] ; 

Optimal( ~Params.activeStaticChannelsMask ) = NaN ;

Update  = [ round( fo0, 6, 'significant' ) ; ...
            round(i0 - 1000*Shim.System.currents) ; ...
            round((i0Aux - Shim.Aux.System.currents), 4, 'significant') ] ;

Update( ~Params.activeStaticChannelsMask )  = NaN ;

Realtime       = round( dx, 4, 'significant' ) ;

% relative power: RIRO / static
relativePower = 100*(dx./x0).^2 ;
relativePower( isnan( relativePower ) ) = 0 ;

Relative_Power = round( relativePower, 4, 'significant' ) ;


T.Correction_Term = char( Correction_Term ) ;
T.Original        = Original ;
T.Optimal         = Optimal ;
T.Update          = Update ;
T.Realtime        = Realtime ;
T.Relative_Power  = Relative_Power ;

if Params.isDisplayingResults
    T    
end

if Params.isSavingResultsTable
    
    writetable( T, [ Params.mediaSaveDir '/shimSystemStats' ]  ) ;

    filename = [ Params.mediaSaveDir '/fieldStats_original' ] ;
    Shim.Field.assessfielddistribution( Params.assessmentVoi, filename ) ;
    
    PredictedShimmedField = Shim.predictshimmedfield() ;
    filename = [ Params.mediaSaveDir '/fieldStats_shimmedPrediction' ] ;
    PredictedShimmedField.assessfielddistribution( Params.assessmentVoi, filename ) ;

    % Shim.Field.write( [Params.mediaSaveDir '/field'], 'nii' ) ; 

    Params.imgSlice     = 7 ; 
    Params.scaling      = [ -100 100 ] ;
    Params.colormap     = 'default' ;

    Params.filename     = [Params.mediaSaveDir '/field_Dc_s' num2str(Params.imgSlice)] ;
    MaRdI.writeimg( Shim.Field.img(:,:,Params.imgSlice), Params ) ;
    
    Params.filename     = [Params.mediaSaveDir '/fieldShimmed_Dc_s' num2str(Params.imgSlice)] ;
    MaRdI.writeimg( PredictedShimmedField.img(:,:,Params.imgSlice), Params ) ;
    
    if Params.isRealtimeShimming
        Params.scaling      = [ -10 10 ] ;

        filename = [ Params.mediaSaveDir '/riroStats_original' ] ;
        Shim.Field.Model.Riro.assessfielddistribution( Params.assessmentVoi, filename ) ;

        PredictedShimmedRiro = Shim.predictshimmedriro() ;

        filename = [ Params.mediaSaveDir '/riroStats_shimmedPrediction' ] ;
        PredictedShimmedRiro.assessfielddistribution( Params.assessmentVoi, filename ) ;

        Params.filename     = [Params.mediaSaveDir '/riro_s' num2str(Params.imgSlice)] ;
        MaRdI.writeimg( Shim.Field.Model.Riro.img(:,:,Params.imgSlice), Params ) ;

        Params.filename     = [Params.mediaSaveDir '/riroShimmed_s' num2str(Params.imgSlice)] ;
        MaRdI.writeimg( PredictedShimmedRiro.img(:,:,Params.imgSlice), Params ) ;

    end

end

function [f, df] = shimcost( x )
    
y = A*x - b;
f = y'*y;
df = 2*A'*y;
    
end

function [C, Ceq] = checknonlinearconstraints_default( solutionVector )
%CHECKNONLINEARCONSTRAINTS_DEFAULT
%
% Checks solution satisfies nonlinear system constraints 
%
% (e.g. for augmented (real-time) problem)
% 
% C(x) <= 0
% (e.g. x = solution vector of shim currents)
% 
% NOTE
%   "Default" since it may be necessary later to introduce
%   system-specific functions.

Ceq = [];

[fo0, i0, i0Aux, dfo, di, diAux] = splitsolutionvector( solutionVector ) ;

% NB: auxiliary shim currents may be empty or all zeros 
x0 = [ fo0; i0; i0Aux] ; 

% dcCurrents ok?
[C, ~] = checknonlinearconstraints_static( x0 ) ;

if Params.isRealtimeShimming
    
    dx = [ dfo; di; diAux] ; 
    
    dxdp = dx/dp ;
    
    % currents for max recorded pressure ok? (e.g. ~inspired state)
    [C_pMax, ~] = checknonlinearconstraints_static( x0 + ( dxdp * pMax ) ) ;
    % currents for min recorded pressure ok? (e.g. ~expired state)
    [C_pMin, ~] = checknonlinearconstraints_static( x0 + ( dxdp * pMin ) ) ;
    
    C = [C ; C_pMax ; C_pMin] ;
    
end

function [C, Ceq] = checknonlinearconstraints_static( x0 )
%CHECKNONLINEARCONSTRAINTS_STATIC
    Ceq = [];

    % checks on :
    %   global frequency shift ( C(1) )
    %   abs current per channel ( C(2:end) )
    C = abs( x0 ) - Params.maxStaticCorrectionPerChannel ;

end

end

function [fo0, i0, i0Aux, dfo, di, diAux] = splitsolutionvector( x ) 
%SPLITSOLUTIONVECTOR
%
% De-concatenates vertically stacked solution vector into 
%
% [fo0, i0, i0Aux, df0, di, diAux]
% 
% ---- Static terms:
%
% fo0
%   scalar Tx frequency offset (relative to original Larmor!)
%
% i0 
%   vector of DC static shim currents for the primary (e.g. multicoil) shim
% 
% i0Aux
%   vector of DC static shim currents for the auxiliary (e.g. scanner) shim
% 
% ---- Dynamic terms:
%
% dfo
%   scalar Tx frequency /global offset shift for respiratory correction
%
% di
%   vector of current differences (inspired - expired) for the primary shim, 
%   such that the time-varying shim current i(t) 
%       = i0 (di/dp)*p(t)
%
%   where di/dp are the relational constants that scale the respiratory
%   tracking data p(t) into real-time shim current corrections
%
% diAux 
%   vector of current differences for the auxiliary (e.g. scanner) shim

% static terms
fo0   = [];
i0    = [];
i0Aux = [];

% dynamic terms
dfo   = [];
di    = [];
diAux = [];

nTerms = length( x ) ;

if Params.isRealtimeShimming
    assert( mod(nTerms,2) == 0, 'Expected solution vector of even length.' ) ; % should be 1 static + 1 dynamic term for every channel
    x0 = x(1:nTerms/2) ;   % static terms
    x1 = x(nTerms/2+1:end) ; % dynamic terms
else
    x0 = x ;
end

fo0 = x0(1) ;
i0  = x0(2:Shim.System.Specs.Amp.nActiveChannels+1) ;

i0Aux = x0(2+Shim.System.Specs.Amp.nActiveChannels:end) ;
assert( length(i0Aux) == Shim.Aux.System.Specs.Amp.nActiveChannels ) ;


if Params.isRealtimeShimming 

    dfo = x1(1) ;
    di  = x1(2:Shim.System.Specs.Amp.nActiveChannels+1) ;

    diAux = x1(2+Shim.System.Specs.Amp.nActiveChannels:end) ;
        
end


end %splitsolutionvector


end %optimizeshimcurrents()
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
function [ img, Hdr ] = mapdbdi( Params )
%MAPDBDI
% 
% [ img, Hdr ] = MAPDBDI( Params ) 
%
% map dB/dI --- change in field [Hz] per unit current (A)
% 
% Params
%   .reliabilityMask 
%       binary array indicating where phase unwrapping should take place 
%
% TODO
%   Clean-up + Documentation

img = [] ;
Hdr = [] ;

Mag         = MaRdI( Params.dataLoadDirectories{1} ) ;
dBdIRaw     = zeros( [Mag.getgridsize Params.nChannels] ) ;

Params.mask = Params.reliabilityMask ;

for iChannel = 1 : Params.nChannels 
    
    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
    
    fieldMaps = zeros( [Mag.getgridsize Params.nCurrents] ) ;

    for iCurrent = 1  : Params.nCurrents

        ImgArray = cell( Params.nEchoes, 2 ) ;
        
        for iEcho = 1 : Params.nEchoes
            % -------
            % PROCESS GRE FIELD MAPS

            ImgArray{iEcho, 1}  = MaRdI( Params.dataLoadDirectories{ iEcho, 1, iCurrent, iChannel +1 }  ) ; % mag
            ImgArray{iEcho, 2}  = MaRdI( Params.dataLoadDirectories{ iEcho, 2, iCurrent, iChannel +1 }  ) ; % phase
        end 

        [Field, Extras] = FieldEval.mapfield( ImgArray, Params ) ;
        
        Field.Hdr.MaskingImage = Field.Hdr.MaskingImage & Params.mask ;
        
        Params.mask = Field.Hdr.MaskingImage ;

        fieldMaps( :,:,:, iCurrent ) = Field.img ;

    end
    
    dBdIRaw(:,:,:, iChannel) = mapdbdi_singlechannel( fieldMaps, Params.currents(iChannel,:), Params.mask, Params ) ;  

end 

% -------
% filter the dB/dI maps
if Params.Filtering.isFiltering

    disp(['Filtering dB/dI maps...'] )
    Params.filterRadius = Params.Filtering.filterRadius ;

    dBdIFiltered = zeros( size( dBdIRaw ) ) ;

    Tmp = Field.copy();
    Tmp.Hdr.MaskingImage = shaver( Params.mask, 1 ) ;

    for iChannel = 1 : Params.nChannels 

        disp(['iChannel ' num2str(iChannel)] )

        tmp = dBdIRaw(:,:,:, iChannel) ;
        tmp( ~Tmp.Hdr.MaskingImage ) = NaN ; % exclude from medfilt3( )
        tmp = medfilt3( tmp, round( Field.getvoxelsize()./Params.Filtering.filterRadius ) ) ;
        tmp( ~Tmp.Hdr.MaskingImage ) = 0 ;

        dBdIFiltered(:,:,:,iChannel) = tmp ;

    end

    Hdr = Field.Hdr ;
    Hdr.MaskingImage = Tmp.Hdr.MaskingImage ;
    img = dBdIFiltered ;

    % -------
    % Extrapolate the dB/dI maps (EXPERIMENTAL)
    if Params.Extension.isExtending

        disp(['Extrapolating fields ...']) ;
        disp(['Channel 1 of ' num2str(Params.nChannels) ] ) ;

        gridSizeImg = size( Params.mask ) ;
        maskFov     = ones( gridSizeImg )  ; % extended spatial support

        [ ~, A, M ] = extendharmonicfield( img(:,:,:, 1) , maskFov, Hdr.MaskingImage, Params.Extension ) ;

        for iChannel = 1 : Params.nChannels 
            disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] ) ;

            reducedField         = Hdr.MaskingImage .* img(:,:,:, iChannel) ;
            extendedField        = reshape( M'*A*reducedField(:), gridSizeImg ) ;
            img(:,:,:, iChannel) = extendedField + reducedField ;
        end

        Hdr.MaskingImage = (sum(abs(img),4))~=0 ; % extended spatial support

    end

else
    Hdr = Field.Hdr ;
    img = dBdIRaw ;
end

function dBdI = mapdbdi_singlechannel( fieldMaps, currents, mask, Params )
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

ShimUse.customdisplay('Fitting dB/dI...' )

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
    
    error('Unimplemented feature');
    % FitResults = ShimCal.fitfieldtocurrent( fieldMaps, currents, mask, Params ) ;
    % dBdI = FitResults.dBdI ;
end


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
function [ img, Hdr ] = calibratereferencemaps( Params )
%CALIBRATEREFERENCEMAPS
% 
% Wraps to ShimOpt.mapdbdi( ) and writes output to disk
% 
% [ img, Hdr ] = CALIBRATEREFERENCEMAPS( Params )

[ img, Hdr ] = ShimOpt.mapdbdi( Params ) ;

disp(['Saving shim reference maps for future use: '])
disp( Params.pathToShimReferenceMaps ) ;

save( Params.pathToShimReferenceMaps, 'img', 'Hdr' ) ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true)
% =========================================================================
function [ img, Hdr ] = loadshimreferencemaps( pathToShimReferenceMaps )
%LOADSHIMREFERENCEMAPS
%
% [ img, Hdr ] = LOADSHIMREFERENCEMAPS( filename )

ShimUse.customdisplay(['\n Preparing for shim ...  \n\n'...
        'Loading shim reference maps from ' pathToShimReferenceMaps '\n\n']) ;

% Loads .mat: struct with fields
%
% .img 
%       linear dB/dI 'current-to-field' operator
% .Hdr
%       dicom Hdr from the reference map acquisition 
RefMaps = load( pathToShimReferenceMaps ) ; 

assert( myisfield( RefMaps, 'img' ) && myisfield( RefMaps, 'Hdr' ), ...
    [ 'Deprecated format. Contents of ' pathToShimReferenceMaps ' are invalid.'] ) ;

img = RefMaps.img ;
Hdr = RefMaps.Hdr ;

end
% =========================================================================

% =========================================================================
% =========================================================================
end

end
