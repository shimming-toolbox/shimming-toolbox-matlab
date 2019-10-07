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
% Inputs (Optional)
%
%   Field 
%       A FieldEval-type object pertaining to the field map to be targetted by
%       shimming.
%
%   Params 
%       A struct of parameters that can possess the following fields:
%
%       .pathToShimReferenceMaps
%           File path to .mat containing shim reference maps (ie basis fields) &
%           .Hdr info
%
% Outputs
%
%   Shim possesses the following properties:
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info Re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%       .Field
%           As defined above in Inputs.
%           Field can be supplied as an input during ShimOpt instantiation,
%           or, at later point by calling Shim.setoriginalfield( Field ) ;
%
%       .Aux
%           .Shim
%               When Shim does not itself refer to the scanner shims, then .Aux 
%               is a ShimOpt object corresponding to the MRI host system.
%       
%       .Model
%           .currents  
%               Optimal shim current vector (i)
%               [units A]
%
%           .field     
%               Optimal shim field from projection of i onto reference maps (Ai)
%               [units Hz]
%
%           .couplingCoefficients
%               For realtime shimming, vector relating shim correction to
%               respiratory state measurement (e.g. ProbeTracking.Data.p)
%               [units: Hz/Pa (Probe), or Hz/rad (Nav) ]
%
%           .dcCurrentsOffsets
%               For realtime shimming, vector of "y-intercept" currents 
%               (e.g currents when ProbeTracking.Data.p = 0)
%               [units A]
%
%       .System
%           .currents
%               
%           .Specs
%               Object of a type of ShimSpecs sub-class (e.g. ShimSpecs_IUGM_Prisma_fit) 
%
% .......
%
% NOTE
%
% ShimOpt is a MaRdI subclass [ShimOpt < MaRdI]
%     
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
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
    ShimmedField; % object of type FieldEval 
    System; %
end
properties( Hidden = true )
    Interpolant ;
    Ref ; % original shim reference maps before interpolation
end
% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( varargin )
%SHIMOPT - Shim Optimization

Shim.img    = [] ;
Shim.Hdr    = [] ;
Shim.Field  = [] ;  
Shim.Model  = [] ;
Shim.Aux    = [] ;
Shim.System = [] ;
Shim.System.Specs = [] ;

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Params = ShimOpt.assigndefaultparameters( Params, Shim.System.Specs ) ;

%% .......
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
    Shim.img     = RefMaps.Shim.img ;
    Shim.Hdr     = RefMaps.Shim.Hdr ;

    Shim.Ref.img = Shim.img ;
    Shim.Ref.Hdr = Shim.Hdr ;

end
    
Shim.Ref.source = 'data' ;

if ~isempty( Field ) 
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
% clear Shim 
clear Shim
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
% CALIBRATEREALTIMEUPDATES asks user to select the intervals of a pair of
% respiratory recordings (e.g. ProbeTracking.Data.t ) corresponding to inspired
% and expired field maps acquisitions.
%
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
%
% NOTE : Possibly deprecated
warning('Possibly deprecated function call: ShimOpt.calibraterealtimeupdates') 

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
        ProbeTracking.selectmedianmeasurement( Params.Inspired.measurementLog ) ; 

    ShimUse.customdisplay( ...
        ['Median measurement : ' num2str( Params.Inspired.medianP )] ) ;

    ShimUse.customdisplay( ['\n ------- \n ' ...
        'Determine median measurement over expired apnea : \n \n '] ) ;

    Params.Expired.medianP = ...
        ProbeTracking.selectmedianmeasurement( Params.Expired.measurementLog ) ; 

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

% Interpolates Shim.img (reference maps) to the grid (voxel positions) of Field
[X,Y,Z] = Shim.Field.getvoxelpositions() ;
Shim.resliceimg( X, Y, Z ) ;

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
%   where UO * respiratory_measurement = currentsUpdate

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
%
% .isAuxIncluded
%   true || false, account for spatial support of the Shim.Aux system?

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
if ~isempty( Shim.Aux )
    DEFAULT_ISAUXINCLUDED      = true ;
else
    DEFAULT_ISAUXINCLUDED      = false ;
end

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

if ~myisfield( Params, 'isAuxIncluded' ) || isempty( Params.isAuxIncluded ) 
    Params.isAuxIncluded = DEFAULT_ISAUXINCLUDED ;
end

if Params.isAuxIncluded
    assert( ~isempty( Shim.Aux ), 'Auxiliary shim system unavailable' )   
    mask = mask & Shim.Aux.getshimsupport() ;
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

% if myisfieldfilled( Shim.Hdr, 'MaskingImage' )
%     shimSupport = shimSupport & Shim.Hdr.MaskingImage ;
% end

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
% =========================================================================
function [] = predictslicewiseshim( Shim, Params )
%PREDICTSLICEWISESHIM
% 
% TEMPORARY FUNCTION FOR SAGITTAL FIELD MAPS:
%   --> Decompose global shim voi into axial segments and shim segments individually

Params.isSavingResultsTable = false ;
Params.isDisplayingResults  = false ;
Params.dataWeights   = [] ;
   
shimVoi = Shim.Field.Hdr.MaskingImage ;

PredictedField = Shim.Field.copy() ;
PredictedRiro  = Shim.Field.Model.Riro.copy() ;

nZ = size( shimVoi, 1 ) ;% i.e. nImageRows, or nAxialSlices

for iZ = 1 : nZ

    axialShimVoi = false( size( shimVoi ) ) ;

    if nnz( shimVoi( iZ, :, : ) ) > 0
    % shim axial slice:
        axialShimVoi(iZ,:,: ) = shimVoi(iZ,:,: ) ;
        Shim.setshimvolumeofinterest( axialShimVoi ) ;
        Shim.setshimvolumeofinterestriro( axialShimVoi ) ;

        [currents] = Shim.optimizeshimcurrents( Params ) ;

        AxialPredictedField = Shim.predictshimmedfield() ;

        display( [ 'Axial slice ' num2str(iZ) 'of ' num2str(nZ) ] )
        display('Field maxima: 1) unshimmed,2) shimmed:') 
        max( abs( Shim.Field.img( axialShimVoi ) ) )
        max( abs( AxialPredictedField.img( axialShimVoi ) ) )

        PredictedField.img(iZ,:,: ) = AxialPredictedField.img(iZ,:,: ) ;
        
        if Params.isRealtimeShimming
            AxialPredictedRiro  = Shim.predictshimmedriro() ;
            PredictedRiro.img(iZ,:,: )  = AxialPredictedRiro.img(iZ,:,: ) ;
            
           display('Riro maxima: 1) unshimmed,2) shimmed:') 
            max( abs( Shim.Field.Model.Riro.img( axialShimVoi ) ) )
            max( abs( AxialPredictedRiro.img( axialShimVoi ) ) )
        end

    end

end

filename = [ Params.mediaSaveDir '/fieldStats_shimmedPrediction_dynamic' ] ;
PredictedField.assessfielddistribution( Params.assessmentVoi, filename ) ;

Params.colormap     = 'default' ;
Params.imgSlice     = 14 ;

Params.scaling      = [ -300 300 ] ;
Params.filename     = [Params.mediaSaveDir '/fieldShimmed_Dc_s_dynamic' num2str(Params.imgSlice)] ;
MaRdI.savefigure( Params.validityMask(:,:, Params.imgSlice ).*PredictedField.img(:,:,Params.imgSlice), Params ) ;
        
if Params.isRealtimeShimming
    filename = [ Params.mediaSaveDir '/riroStats_shimmedPrediction_dynamic' ] ;
    PredictedRiro.assessfielddistribution( Params.assessmentVoi, filename ) ;

    Params.scaling      = [ -20 20 ] ;
    Params.filename     = [Params.mediaSaveDir '/riro_s_dynamic' num2str(Params.imgSlice)] ;
    MaRdI.savefigure( Params.validityMask(:,:, Params.imgSlice ).*PredictedRiro.img(:,:, Params.imgSlice), Params ) ;
end

% reset shim voi
Shim.setshimvolumeofinterest( shimVoi ) ;
Shim.setshimvolumeofinterestriro( shimVoi ) ;
pause(5)
close all

end
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

if myisfield( Shim.Model, 'Tx' ) && myisfieldfilled( Shim.Model.Tx, 'imagingFrequency' ) 
    
    df0 = Shim.Model.Tx.imagingFrequency - Shim.Field.getimagingfrequency() ; % frequency shift [units: Hz]
    
    if ( abs( df0 ) > 0.1 )
        PredictedField.Hdr.ImagingFrequency = Shim.Model.Tx.imagingFrequency*1E-6 ;  % [units: MHz]
        PredictedField.img = PredictedField.img + df0 ;
    end
end

if isa( Shim.Aux, 'ShimOpt' ) && myisfield( Shim.Aux, 'Model' ) && myisfieldfilled( Shim.Aux.Model, 'currents' )

    auxCurrentsUpdate    = Shim.Aux.Model.currents - Shim.Aux.System.currents ;
    Shim.Aux.Model.field = Shim.Aux.forwardmodelshimcorrection( auxCurrentsUpdate ) ;
    PredictedField.img   = PredictedField.img + Shim.Aux.Model.field ;
end

end
% =========================================================================
function [ PredictedRiro ] = predictshimmedriro( Shim, p )
%PREDICTSHIMMEDRIRO     
%
% [PredictedRiro] = PREDICTSHIMMEDRIRO( Shim ) ;
% [PredictedRiro] = PREDICTSHIMMEDRIRO( Shim, dp ) ;
% 
% Returns a FieldEval-type object

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

if ~myisfieldfilled( Shim.Model, 'couplingCoefficients' )  
    warning('No RIRO correction: Returning input model RIRO') ;
    PredictedRiro = Shim.Field.Model.Riro.copy() ;
else    
    % all field + correction terms scaled by dp :
    if nargin < 2
        PredictedRiro = Shim.Field.getriro( ) ;
        % (dp: recorded inspired - expired pressure difference)
        dp = Shim.Field.Model.Riro.Aux.Data.p ; 
    else
        PredictedRiro = Shim.Field.getriro( p ) ;
        dp = Shim.Field.Aux.debias( p ) ; 
    end
    
    Shim.Model.riro   = dp*Shim.forwardmodelshimcorrection( Shim.Model.couplingCoefficients ) ;
    PredictedRiro.img = PredictedRiro.img + Shim.Model.riro ;

    if myisfield( Shim.Model, 'Tx' ) && myisfieldfilled( Shim.Model.Tx, 'couplingCoefficients' ) 
        if ( Shim.Model.Tx.couplingCoefficients ~= 0 )  
           warning('Untested: TX RIRO adjustment') 
        end
       PredictedRiro.img = PredictedRiro.img + dp*Shim.Model.Tx.couplingCoefficients ;
    end

    if isa( Shim.Aux, 'ShimOpt' ) && myisfieldfilled( Shim.Aux.Model, 'couplingCoefficients' ) 
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
DEFAULTS.isReturningPseudoInverse       = false ;
DEFAULTS.isRealtimeShimming             = false ;
DEFAULTS.isDisplayingResults            = false ;
DEFAULTS.isSavingResultsTable           = true ;

DEFAULTS.isOptimizingAux                = false ;
DEFAULTS.mediaSaveDir                   = [ './shimopt_results_' datestr(now,30) '/' ] ;

DEFAULTS.isOptimizingStaticTxFrequency  = false ;
DEFAULTS.isOptimizingDynamicTxFrequency = false ;

% TODO: Create new class to handle all Tx related aspects. These defaults are for the IUGM Prisma fit
DEFAULTS.minTxFrequency = 123100100 ; % [units: Hz]
DEFAULTS.maxTxFrequency = 123265000 ; % [units: Hz]

%% -----
% Check inputs + assign parameters 
if nargin == 1 || isempty( Params ) 
    Params.dummy = [] ; % Using default parameters
end

Params = assignifempty( Params, 'isReturningPseudoInverse', DEFAULTS.isReturningPseudoInverse ) ;
Params = assignifempty( Params, 'isRealtimeShimming', DEFAULTS.isRealtimeShimming ) ;
Params = assignifempty( Params, 'isDisplayingResults', DEFAULTS.isDisplayingResults ) ;
Params = assignifempty( Params, 'isSavingResultsTable', DEFAULTS.isSavingResultsTable ) ;

Params = assignifempty( Params, 'assessmentVoi', Shim.Field.Hdr.MaskingImage ) ;

if isa( Shim.Aux, 'ShimOpt' )
    if ~myisfieldfilled(Params, 'isJointOptimization')  
        Params.isJointOptimization = true ; 
    end
else
    Params.isJointOptimization = false ;
end

% check imaging frequency
Params = assignifempty( Params, 'maxTxFrequency', DEFAULTS.maxTxFrequency ) ;
Params = assignifempty( Params, 'minTxFrequency', DEFAULTS.minTxFrequency ) ;

f0 = Shim.Field.getimagingfrequency() ;

if ( f0 > DEFAULTS.maxTxFrequency ) || ( f0 < DEFAULTS.minTxFrequency )
    error('Imaging frequency of the input Field is outside the expected range. See input Params.minTxFrequency and Params.maxTxFrequency' ) ;
end 

% max +/- freq. shift relative to original Larmor
Params = assignifempty( Params, 'maxPositiveTxFrequencyShift', DEFAULTS.maxTxFrequency - f0 ) ;
Params = assignifempty( Params, 'maxNegativeTxFrequencyShift', DEFAULTS.minTxFrequency - f0 ) ;


if Params.isJointOptimization 
    % total active channels across systems
    Params.nActiveChannels = 1 + Shim.System.Specs.Amp.nActiveChannels + Shim.Aux.System.Specs.Amp.nActiveChannels ; % +1 for freq. adjust
    
    if ~myisfieldfilled(Params, 'activeStaticChannelsMask') 
        Params.activeStaticChannelsMask = [ DEFAULTS.isOptimizingStaticTxFrequency ; ...
                                            Shim.System.Specs.Amp.staticChannels(:) ; ...
                                            Shim.Aux.System.Specs.Amp.staticChannels(:) ; ] ;
    end
        
    if ~myisfieldfilled(Params, 'activeDynamicChannelsMask') 
        Params.activeDynamicChannelsMask = [ DEFAULTS.isOptimizingDynamicTxFrequency ; ...
                                             Shim.System.Specs.Amp.dynamicChannels(:) ; ...
                                             Shim.Aux.System.Specs.Amp.dynamicChannels(:) ; ] ;
    end
    
    if ~myisfieldfilled(Params, 'maxCorrectionPerChannel')  
       
        % 'tmp' because this is technically a TODO: 
        % since the optimization solves for the Tx frequency *shift*, it can shift
        % up or down relative to Larmor.  For now, limit the max shift equally
        % in both up/down directions, but in the future both should be specifically
        % accounted for (only comes into play in checknonlinearconstraints_static() )
        % (in practice the limits of this tmp solution seem unlikely to arise, i.e. adjusting f0 by many kHz?)
        tmpMaxFreqShift = min(abs( [ Params.maxPositiveTxFrequencyShift Params.maxNegativeTxFrequencyShift ] )) ;

        Params.maxCorrectionPerChannel = [ tmpMaxFreqShift ; Shim.System.Specs.Amp.maxCurrentPerChannel ; Shim.Aux.System.Specs.Amp.maxCurrentPerChannel ]; 

    end

else
    % total active channels across systems
    Params.nActiveChannels = 1 + Shim.System.Specs.Amp.nActiveChannels ; % +1 for freq. adjust
    
    if ~myisfieldfilled(Params, 'activeStaticChannelsMask') 
        Params.activeStaticChannelsMask = [ DEFAULTS.isOptimizingStaticTxFrequency ; ...
                                            Shim.System.Specs.Amp.staticChannels(:) ; ] ;
    end
        
    if ~myisfieldfilled(Params, 'activeDynamicChannelsMask') 
        Params.activeDynamicChannelsMask = [ DEFAULTS.isOptimizingDynamicTxFrequency ; ...
                                             Shim.System.Specs.Amp.dynamicChannels(:) ; ] ;
    end
    
    if ~myisfieldfilled(Params, 'maxCorrectionPerChannel')  
       
        % see todo note above for ~isempty(Shim.Aux) case:
        tmpMaxFreqShift = min(abs( [ Params.maxPositiveTxFrequencyShift Params.maxNegativeTxFrequencyShift ] )) ;

        Params.maxCorrectionPerChannel = [ tmpMaxFreqShift ; Shim.System.Specs.Amp.maxCurrentPerChannel ; ]; 

    end
end

Params = assignifempty( Params, 'minCorrectionPerChannel', -Params.maxCorrectionPerChannel ) ;

assert( ( length( Params.maxCorrectionPerChannel ) == length( Params.minCorrectionPerChannel ) ) && ...
        ( length( Params.maxCorrectionPerChannel ) == Params.nActiveChannels ), ...
    'Shim system limits (Params.maxCorrectionPerChannel and Params.minCorrectionPerChannel) must possess an entry for the TX frequency followed by entries for each shim channel (primary and, if in use, auxiliary shims).' ) ;

% Set min/max static correction to 0 for inactive static channels
% NB: some elements of Params.min/maxCorrectionPerChannel may be +/-Inf (so
% don't just multiply by the logical array activeChannelsMask --> yields NaN!)
Params.minStaticCorrectionPerChannel = zeros( size( Params.activeStaticChannelsMask ) ) ;
Params.minStaticCorrectionPerChannel( Params.activeStaticChannelsMask ) = Params.minCorrectionPerChannel( Params.activeStaticChannelsMask ) ;

Params.maxStaticCorrectionPerChannel = zeros( size( Params.activeStaticChannelsMask ) ) ;
Params.maxStaticCorrectionPerChannel( Params.activeStaticChannelsMask ) = Params.maxCorrectionPerChannel( Params.activeStaticChannelsMask ) ;

% if ~myisfieldfilled(Params, 'regularizationParameter')  
%     Params.regularizationParameter = DEFAULT_REGULARIZATIONPARAMETER ;
% end

Params = assignifempty( Params, 'mediaSaveDir', DEFAULTS.mediaSaveDir ) ;
mkdir( Params.mediaSaveDir ) ;


if Params.isRealtimeShimming

    assert( myisfield( Shim.Field.Model, 'Riro') && ~isempty( Shim.Field.Model.Riro.img ), ...
        'Real-time shim optimization requires a valid FieldEval object pertaining to respiration-induced resonance offsets (RIRO) in Shim.Field.Model.Riro.img' ) ;
   
    if ~myisfield(Params, 'pMax') || isempty( Params.pMax ) 
        pMax = Shim.Field.Model.Riro.Aux.Specs.limits(2) ;
    else
        pMax = Params.pMax ; % max relevant/recorded pressure (e.g. ~inspired state)
    end
    
    if ~myisfield(Params, 'pMin') || isempty( Params.pMin ) 
        pMin = Shim.Field.Model.Riro.Aux.Specs.limits(2) ;
    else
        pMin = Params.pMin ; % min relevant/recorded pressure (e.g. ~expired state) 
    end

    % dp : delta pressure shift to which input Riro.img has been scaled
    % e.g. Riro = dp * ( Inspired_Field - Expired_Field )/( inspired_pressure - expired_pressure )
    dp = Shim.Field.Model.Riro.Aux.Data.p ; 
    
end

%%-------------------------------------------------------------------------------
% define optimization terms:
%--------------------------------------------------------------------------------

nVoxels = Shim.Field.getnumberofvoxels() ; 

%% -----
% define matrix of data-weighting coefficients for the static optimization: W0
if ~myisfield( Params, 'dataWeights' ) || isempty( Params.dataWeights ) 
    W0 = speye( nVoxels, nVoxels ) ;
else
    assert( numel( Params.dataWeights ) == nVoxels, ...
        'Params.dataWeights must have the same number of elements (voxels) as Shim.Field.img' )

    if ( size( Params.dataWeights, 1 ) ~= nVoxels ) || ( size( Params.dataWeights, 2) ~= nVoxels )
        
        W0 = spdiags( Params.dataWeights(:), 0, nVoxels, nVoxels ) ;

    end
end

% truncated (VOI-masked) weighting operator for static shim: MW0
MW0 = Shim.gettruncationoperator()*W0 ; 

%% -----
% define matrix of data-weighting coefficients for real-time (RIRO) correction: W1
if Params.isRealtimeShimming 
    if ~myisfield( Params, 'dataWeightsRiro' ) || isempty( Params.dataWeightsRiro ) 
        W1 = W0 ;
    else
        assert( numel( Params.dataWeightsRiro ) == nVoxels, ... 
            'Params.dataWeightsRiro must have the same number of elements (voxels) as Shim.Field.img' )

        if ( size( Params.dataWeightsRiro, 1 ) ~= nVoxels ) || ( size( Params.dataWeightsRiro, 2) ~= nVoxels )
            
            W1 = spdiags( Params.dataWeightsRiro(:), 0, nVoxels, nVoxels ) ;

        end
    end

    % truncated (VOI-masked) weighting operator for riro shim: MW1
    MW1 = Shim.gettruncationoperatorriro()*W1 ; 
end

%% -----
% define off-resonance correction operator : A

% global resonance offset operator (Tx frequency) : A_tx
A_tx = ones( nVoxels, 1 ) ;

% primary shim (e.g. scanner's spherical harmonic coils) (i.e. current-to-field) operator : A_p
A_p = Shim.getshimoperator() ; 

% combined: static off-resonance correction operator: A0

if Params.isJointOptimization
    % auxiliary shim operator (e.g. multi-coil array)
    A_aux = Shim.Aux.getshimoperator() ;
    A0    = [ MW0*A_tx MW0*A_p  MW0*A_aux] ;
else
    A0    = [ MW0*A_tx MW0*A_p ] ;
end

if Params.isRealtimeShimming
    
    activeChannelsMask = [ Params.activeStaticChannelsMask ; Params.activeDynamicChannelsMask ] ;

    if Params.isJointOptimization
        A1 = [ MW1*A_tx MW1*A_p MW1*A_aux ] ;
    else
        A1 = [ MW1*A_tx MW1*A_p ] ;
    end

    % left half applies to the static field, right half to RIRO 
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





%% -------
% Define the solution (data) vector: b

% Solve for residual field offset (bx0) *without* existing shim fields (bs) 
% *from shim channels that will be optimized* (i.e. if shim settings are to
% remain constant anyway, there is no need to subtract out their contribution
% from the field).
% 
% This way, optimizeshimcurrents() solves for the absolute shim currents,
% rather than a shift, which simplifies handling of the min/max constraints on
% the correction terms.

% Retain original corrections for 'initial guess' x0 of conjugate gradient solver:
%
% e.g. 1st element: corresponds to an initial guess of zero for the optimal frequency 
% shift relative to the input field map

if Params.isJointOptimization
    x0  = [ 0 ; Shim.System.currents ; Shim.Aux.System.currents ] ;
    bx0 = Shim.Field.img(:) - [A_tx A_p A_aux]*x0 ;
else
    x0  = [ 0 ; Shim.System.currents ] ;
    bx0 = Shim.Field.img(:) - [A_tx A_p]*x0 ;
end

x0( ~Params.activeStaticChannelsMask ) = 0 ;

%% -----
% concatenate/stack bx0 onto the RIRO vector db to define the complete solution vector b
if Params.isRealtimeShimming

    % field shift from RIRO  
    db = Shim.Field.Model.Riro.img(:) ;

    %  stacked data vector : dc field, RIRO field shift -- vertically concatenated
    b = [ MW0*-bx0 ; MW1*-db ] ;

    % initial guess of zero for the all RIRO corrections 
    x0 = [ x0 ; zeros( size(x0) ) ] ; 
else 
    b  = MW0*-bx0 ;
end
    
%%-------------------------------------------------------------------------------
% perform optimization 
%--------------------------------------------------------------------------------

%% -----
% linear optimization:

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

%%-------------------------------------------------------------------------------
% organize results
%--------------------------------------------------------------------------------
[ fo0, i0, i0Aux, dfo, di, diAux] = splitsolutionvector( x ) ;

Shim.Model.Tx.imagingFrequency    = fo0 + (Shim.Field.Hdr.ImagingFrequency*1E6) ; % DC avg. [units: Hz]

Shim.Model.currents               = i0 ;

Corrections.activeStaticChannelsMask  = Params.activeStaticChannelsMask ;
Corrections.activeDynamicChannelsMask = Params.activeDynamicChannelsMask ;

%% -----
% replace inactive static terms with their initial values:

if Params.isJointOptimization
    Shim.Aux.Model.currents             = i0Aux ;
    Params.activeStaticChannelsMaskAux  = Params.activeStaticChannelsMask(2+length(i0):end) ;
    Shim.Aux.Model.currents( ~Params.activeStaticChannelsMaskAux ) = Shim.Aux.System.currents( ~Params.activeStaticChannelsMaskAux ) ;

    Corrections.static = [ Shim.Model.Tx.imagingFrequency; Shim.Model.currents; Shim.Aux.Model.currents ] ;
else

    Params.activeStaticChannelsMaskP   = Params.activeStaticChannelsMask(2:(length(i0)+1)) ;
    Shim.Model.currents( ~Params.activeStaticChannelsMaskP ) = Shim.System.currents( ~Params.activeStaticChannelsMaskP ) ;
    
    Corrections.static = [ Shim.Model.Tx.imagingFrequency; Shim.Model.currents ] ;
end

%% -----
if Params.isRealtimeShimming
    Shim.Model.Tx.couplingCoefficients  = dfo/dp ; % resp-induced shift [units: Hz/unit-pressure]
    Shim.Model.couplingCoefficients     = di/dp ;

    if Params.isJointOptimization
        Shim.Aux.Model.couplingCoefficients = diAux/dp ;
        Corrections.realtime = [Shim.Model.Tx.couplingCoefficients; Shim.Model.couplingCoefficients; Shim.Aux.Model.couplingCoefficients] ; 
    else
        Corrections.realtime = [Shim.Model.Tx.couplingCoefficients; Shim.Model.couplingCoefficients] ; 
    end
end

PredictedShimmedField = Shim.predictshimmedfield() ;

if Params.isRealtimeShimming
    PredictedShimmedRiro = Shim.predictshimmedriro() ;
end

%%-------------------------------------------------------------------------------
% results summary 
%--------------------------------------------------------------------------------
T = tableoptimizationresults( Corrections ) ;

if Params.isSavingResultsTable
    writetable( T, [ Params.mediaSaveDir '/shimSystemStats.csv' ]  ) ;

    Shim.Field.assessfielddistribution( Params.assessmentVoi, [ Params.mediaSaveDir '/fieldStats_original.csv' ] ) ;

    PredictedShimmedField.assessfielddistribution( Params.assessmentVoi, [ Params.mediaSaveDir '/fieldStats_shimmedPrediction.csv' ] ) ;

    if Params.isRealtimeShimming
        Shim.Field.Model.Riro.assessfielddistribution( Params.assessmentVoi, [ Params.mediaSaveDir '/riroStats_original' ] ) ;

        PredictedShimmedRiro.assessfielddistribution( Params.assessmentVoi, [ Params.mediaSaveDir '/riroStats_shimmedPrediction' ] ) ;
    end
end

%%-------------------------------------------------------------------------------

%--------------------------------------------------------------------------------

if Params.isDisplayingResults

    fprintf(['\n' '-----\t' 'Optimization Results' '\t-----' '\n\n']) ;
    display(T)

    Params.imgSlice     = round( size( Shim.Field.img, 3 )/2 ) ; 
    
    Params.scaling      = [ 0 1 ] ;
    Params.colormap     = 'gray' ;

    Params.filename     = [Params.mediaSaveDir '/shimVoi_s' num2str(Params.imgSlice)] ;
    MaRdI.savefigure( Params.assessmentVoi(:,:, Params.imgSlice ), Params ) ;
    
    Params.scaling      = [ -100 100 ] ;
    Params.colormap     = 'default' ;
    
    Params.filename     = [Params.mediaSaveDir '/field_Dc_s' num2str(Params.imgSlice)] ;
    MaRdI.savefigure( Params.assessmentVoi(:,:, Params.imgSlice ).* Shim.Field.img(:,:,Params.imgSlice), Params ) ;
    
    Params.filename     = [Params.mediaSaveDir '/fieldShimmed_Dc_s' num2str(Params.imgSlice)] ;
    MaRdI.savefigure( Params.assessmentVoi(:,:, Params.imgSlice ).*PredictedShimmedField.img(:,:,Params.imgSlice), Params ) ;
        
    NiiOptions.filename = [Params.mediaSaveDir 'fieldB0'] ;
    nii( Shim.Field.img, NiiOptions ) ;
    
    NiiOptions.filename = [Params.mediaSaveDir 'fieldB0_shimmed'] ;
    nii( PredictedShimmedField.img, NiiOptions ) ;
    
    if Params.isRealtimeShimming
       
        Params.scaling      = [ -dp dp ]/3 ;
        Params.filename     = [Params.mediaSaveDir '/riro_s' num2str(Params.imgSlice)] ;
        MaRdI.savefigure( Params.assessmentVoi(:,:, Params.imgSlice ).*Shim.Field.Model.Riro.img(:,:,Params.imgSlice), Params ) ;

        Params.filename     = [Params.mediaSaveDir '/riroShimmed_s' num2str(Params.imgSlice)] ;
        MaRdI.savefigure( Params.assessmentVoi(:,:, Params.imgSlice ).*PredictedShimmedRiro.img(:,:,Params.imgSlice), Params ) ;
        
        NiiOptions.filename = [Params.mediaSaveDir 'riro'] ;
        nii( Shim.Field.Model.Riro.img, NiiOptions ) ;
        
        NiiOptions.filename = [Params.mediaSaveDir 'riro_shimmed'] ;
        nii( PredictedShimmedRiro.img, NiiOptions ) ;

    end
    pause(5)
    close all
end

%%-------------------------------------------------------------------------------
% Local functions
%--------------------------------------------------------------------------------

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

C( C == -Inf ) = 0 ;

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
    assert( mod(nTerms,2) == 0, 'Expected shim currents vector of even length.' ) ; % should be 1 static + 1 dynamic term for every channel
    x0 = x(1:nTerms/2) ;   % static terms
    x1 = x(nTerms/2+1:end) ; % dynamic terms
else
    x0 = x ;
end

fo0 = x0(1) ;
i0  = x0(2:Shim.System.Specs.Amp.nActiveChannels+1) ;

if Params.isJointOptimization
    i0Aux = x0(2+Shim.System.Specs.Amp.nActiveChannels:end) ;
    assert( length(i0Aux) == Shim.Aux.System.Specs.Amp.nActiveChannels ) ;
end

if Params.isRealtimeShimming 

    dfo = x1(1) ;
    di  = x1(2:Shim.System.Specs.Amp.nActiveChannels+1) ;

    if Params.isJointOptimization
        diAux = x1(2+Shim.System.Specs.Amp.nActiveChannels:end) ;
    end
end


end %splitsolutionvector

function [T] = tableoptimizationresults( Corrections )
%TABLEOPTIMIZATIONRESULTS      Table of shim results
%
% T = TABLEOPTIMIZATIONRESULTS( Corrections )
%
% T is a table where the 1st column contains the shim terms (channel names [+
% units]), 2nd column contains their original settings, 3rd is the optimal
% (static) shim terms, and 4th the optimal - original difference.
%
% Terms excluded from the optimization are assigned entries of NaN.
%
% If realtime shimming was enabled, the 5th column contains entries with its
% corresponding correction term's units, per unit (respiratory) probe
% intensity.  (e.g. 0.01 A/probe-unit-intensity). The 6th column contains the
% squared ratio of a channel's real-time correction (scaled by the probe RMS
% intensity) to its static correction, multiplied by 100%.
%
% e.g. T =
% 
% Correction_Term  | Original  | Optimal | Update | Realtime | Relative_Power
%   Tx Freq. [Hz]  | 123259218 |   NaN   |  NaN   |    NaN   |       0
%    Ch.1    [A]   |    0      |  1.23   |  1.23  |   0.01   |     0.05
%         .        |    .      |    .    |   .    |    .     |      .
%         .        |    .      |    .    |   .    |    .     |      .
%         .        |    .      |    .    |   .    |    .     |      .

T = table ;

Correction_Term = {'Tx Freq. [Hz]'} ;

for iCh = 1 : Shim.System.Specs.Amp.nActiveChannels 
    Correction_Term(end+1) = { [ Shim.System.Specs.Id.systemName ' ' ...
                                 Shim.System.Specs.Id.channelNames{iCh} ' ' ...
                                 Shim.System.Specs.Id.channelUnits{iCh}  ] } ;
end
    
Original = [ Shim.Field.getimagingfrequency ; Shim.System.currents ; ] ;

if isa( Shim.Aux, 'ShimOpt' )

    for iCh = 1 : Shim.Aux.System.Specs.Amp.nActiveChannels 
        Correction_Term(end+1) = { [ Shim.Aux.System.Specs.Id.systemName ' ' ...
                                     Shim.Aux.System.Specs.Id.channelNames{iCh} ' ' ...
                                     Shim.Aux.System.Specs.Id.channelUnits{iCh}  ] } ;
    end

    Original = [ Original ;  Shim.Aux.System.currents ;  ] ;
end

nChannelsTotal = length( Original ) ;

if length( Corrections.static ) < nChannelsTotal 
   % append NaNs to indicate Aux channels were excluded from the optimization
   Optimal = [ Corrections.static ; nan(nChannelsTotal - length(Corrections.static), 1) ] ;
   activeStaticChannelsMask  = [ Corrections.activeStaticChannelsMask; false(nChannelsTotal - length(Corrections.static), 1) ] ;
   activeDynamicChannelsMask = [ Corrections.activeDynamicChannelsMask ; false(nChannelsTotal - length(Corrections.static), 1) ] ;
else
   Optimal = [ Corrections.static ] ;
   activeStaticChannelsMask  = Corrections.activeStaticChannelsMask ;
   activeDynamicChannelsMask = Corrections.activeDynamicChannelsMask ;
end

Optimal( ~activeStaticChannelsMask ) = NaN ;

Update = Optimal - Original ;

if myisfield( Corrections, 'realtime' ) 
    Realtime = Corrections.realtime ;
    Realtime( ~activeDynamicChannelsMask ) = NaN ;

    Relative_Power = 100*( ( ( Shim.Field.Model.Riro.Aux.Data.p * Realtime )./Optimal ).^2 ) ;
    Relative_Power( isnan( Relative_Power ) ) = 0 ;
    
    T.Correction_Term = char( Correction_Term ) ;
    T.Original        = round( Original, 3 ) ;
    T.Optimal         = round( Optimal, 3 ) ;
    T.Update          = round( Update, 3 ) ;
    T.Realtime        = round( Realtime, 5, 'significant' ) ;
    T.Relative_Power  = round( Relative_Power, 5, 'significant' ) ;
else
    T.Correction_Term = char( Correction_Term ) ;
    T.Original        = round( Original, 3 ) ;
    T.Optimal         = round( Optimal, 3 ) ;
    T.Update          = round( Update, 3 ) ;
end

end %tableoptimizationresults()


end %optimizeshimcurrents()
% =========================================================================
function [T] = tableshim( Shim, Corrections )
%TABLESHIM      Table of shim terms
%
% T = TABLESHIM( Shim )
% T = TABLESHIM( Shim, Corrections )
%
% When the only input argument is a ShimOpt-type object, T is a table where
% the first column contains the shim terms (channel names), and the second column 
% contains their current setting. 
%
% e.g. T =
%
% Correction_Term  | Original 
%   Tx Freq. [Hz]  | 123259218
%    Ch.1    [A]   |    0
%         .        |    .
%         .        |    .
%         .        |    .
%    Ch.N    [A]   |    0
 
Correction_Term = {'Tx Freq. [Hz]'} ;

for iCh = 1 : Shim.System.Specs.Amp.nActiveChannels 
    Correction_Term(end+1) = { [ Shim.System.Specs.Id.systemName ' ' ...
                                 Shim.System.Specs.Id.channelNames{iCh} ' ' ...
                                 Shim.System.Specs.Id.channelUnits{iCh}  ] } ;
end
    
Original = [ Shim.Field.getimagingfrequency ; Shim.System.currents ; ] ;

if isa( Shim.Aux, 'ShimOpt' )

    for iCh = 1 : Shim.Aux.System.Specs.Amp.nActiveChannels 
        Correction_Term(end+1) = { [ Shim.Aux.System.Specs.Id.systemName ' ' ...
                                     Shim.Aux.System.Specs.Id.channelNames{iCh} ' ' ...
                                     Shim.Aux.System.Specs.Id.channelUnits{iCh}  ] } ;
    end

    Original = [ Original ;  Shim.Aux.System.currents ;  ] ;
end

T = table ;
    
T.Correction_Term = char( Correction_Term ) ;    
T.Original        = round( Original, 3 ) ;

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

magAtEcho = M0 .* exp( -targetEchoTime/t2star ) ; % forecasted undistorted magnitude at targetEchoTime 
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
%MAPDBDI    map dB/dI : field shift [Hz] per unit current (A)
% 
% [ img, Hdr ] = MAPDBDI( Params ) 

%
% 
% Params
%
%   For unwrapping (see MaRdI.unwrapphase for details):
%
%       .threshold 
%       .unwrapper
dbstop in ShimOpt at 1797
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

disp( ['Performing shim calibration...' ] )        


    Fields = cell( nCurrents, 1, 8 ) ;
for iChannel = 1 : nChannels 
    
    disp(['Channel ' num2str(iChannel) ' of ' num2str(nChannels) ] )        

    % for iCurrent = 1  : nCurrents
    %     Fields{ iCurrent, 1, iChannel } = FieldEval( Params.dataLoadDirectories{ iCurrent, 1, iChannel }, ...
    %                                                  Params.dataLoadDirectories{ iCurrent, 2, iChannel }, Params  ) ; 
    % end

    dB = Fields{2,1,iChannel}.img - Fields{1,1,iChannel}.img ;
    dI = Params.currents(iChannel,2) - Params.currents(iChannel, 1) ;

    dBdI(:,:,:,iChannel) = dB/dI ;
end

for iChannel = 1 : nChannels 
    dBdI0 = dBdI(:,:,:,iChannel) ;
    dBdI2(:,:,:,iChannel) = dBdI0 - dBdI0(iR) ;
end

for iChannel = 1 : nChannels 
    
    disp(['Channel ' num2str(iChannel) ' of ' num2str(nChannels) ] )        

    Img = cell( nCurrents, 2 ) ;

    for iCurrent = 1  : nCurrents
        Img{ iCurrent, 1 } = MaRdI( Params.dataLoadDirectories{ iCurrent, 1, iChannel }  ) ; % mag
        Img{ iCurrent, 2 } = MaRdI( Params.dataLoadDirectories{ iCurrent, 2, iChannel }  ) ; % phase
    end
    
    % ASSUMES :
    %   NO TX FREQUENCY SHIFT BETWEEN ACQUISITIONS
    % TODO
    %   adjust Field in the event Tx Freq. changed between scans/shim currents
    
    assert( all( Img{1,1}.getechotime() ) == all( Img{2,1}.getechotime() ), ...
        'Expected identical echo times across acquisitions' )
    assert( MaRdI.compareimggrids( Img{1,1}, Img{2,1} ), ...
         'Expected identical voxel positions across acquisitions' )

    if iChannel == 1
        dBdI    = zeros( [ Img{1,1}.getgridsize Params.nChannels ] ) ;
        [X,Y,Z] = Img{1,1}.getvoxelpositions() ;
    else
        [X1,Y1,Z1] = Img{1,1}.getvoxelpositions ;
        if ~MaRdI.compareimggrids( X,Y,Z, X1,Y1,Z1)
            error('Voxel positions should remain constant across all acquisitions' ) ;
        end
    end

    TE        = Img{1,1}.getechotime() ;
    nEchoes   = length( TE ) ;
    
    Mag       = Img{ 1, 1 }.copy() ;
    Phase     = Img{ 1, 2 }.copy() ;

    % complex difference between the 2 acquisitions:
    Phase.img = angle( ( Img{ 2, 1 }.img .*exp(i*Img{ 2, 2 }.img) ) .* conj( Img{ 1, 1 }.img .*exp(i*Img{ 1, 2 }.img) ) ) ;

    %% -----
    % create reliability MaskingImage for phase unwrapping based on the magnitude:
    % average across the 2 acquisitions and normalize
    Mag.img = ( Img{1,1}.img + Img{2,1}.img )/2 ; 

    echoMasks = ones( size( Mag.img ) ) ;

    for iEcho = 1 : nEchoes
        tmp = Mag.img(:,:,:,iEcho) ;
        Mag.img(:,:,:,iEcho) = Mag.img(:,:,:,iEcho) /max( tmp(:) ) ;

        echoMasks(:,:,:,iEcho) = TE(iEcho)*ones( Fields.getgridsize() ) < TE_max ;
        Phase.Hdr.MaskingImage(:,:,:,iEcho) = echoMasks(:,:,:,iEcho) .* ( Mag.img(:,:,:,iEcho) > Params.threshold ) ;
    end

    Phase.Hdr.MaskingImage = Mag.img > Params.threshold ;
    
    Phase.unwrapphase( Mag, Params ) ; 

    Fields = Phase.copy() ;
    Fields.scalephasetofrequency() ;
    %% -----
    % Field.img now contains nEchoes estimates of the field change per
    % unit-change in shim amplitude along its 4th dimension. 
    % We now need to combine them, accounting for problem regions where too
    % much intervoxel dephasing may have occurred
   
    % estimate the inter-voxel dephasing 

    % Fields now contains nEchoes estimates of the effective applied field introduced by the shim
    % before combining the estimates, exclude potential problem regions where excessive intervoxel
    % dephasing may have occured by using the the applied field estimate at the 1st echo (minimal dephasing):
    
    % use the absolute spatial gradient in the phase at the first echo (minimal dephasing)    
    [Dx, Dy, Dz] = createdifferenceoperators( Fields.getgridsize(), Fields.getvoxelspacing(), 1 ) ; 

    df1   = Fields.img(:,:,:,1 ) ;
    % effective gradient [units: Hz/mm]
    Gx_eff = reshape( abs(Dx*df1(:)), Fields.getgridsize() ) ;
    Gy_eff = reshape( abs(Dy*df1(:)), Fields.getgridsize() ) ;
    Gz_eff = reshape( abs(Dz*df1(:)), Fields.getgridsize() ) ;

    voxelSpacing = Fields.getvoxelspacing() ;

    % 2pi phase dispersion across a voxel will result in total signal loss;
    % however, the need for unwrapping imposes a stricter limit on phase
    % accumulation: wraps cannot occur more frequently than ~ 
    TE_max = 1000*pi./( voxelSpacing(1)*Gx_eff + voxelSpacing(2)*Gy_eff + voxelSpacing(3)*Gz_eff )/4 ;

    echoWeights = zeros( size( Fields.img ) ) ;
    
    for iEcho = 1 : nEchoes
        % exclude voxels wherever the estimated max TE is exceeded: 
        echoMask = TE(iEcho)*ones( Fields.getgridsize() ) < TE_max ;

        % average reference map to be weighted by the echo-time x avg. magnitude(TE)
        echoWeights(:,:,:,iEcho) = TE(iEcho)*echoMask.*Mag.img(:,:,:,iEcho) ;
    end

    %ESTIMATE + SAVE DBDI SNR

    % Weighted-average of across echoes, converted to Hz/unit-shim:
    dBdI(:,:,:,iChannel) = sum( echoWeights .* Fields.img, 4 )./sum( echoWeights, 4 )...
        /( Params.currents( iChannel, 2) - Params.currents( iChannel, 1) ) ;
    
    % % smoothing: [3x3x3] magnitude-weighted Gaussian filtering
    %
    % Fields.filter( Mag.img ) ; 
    
        
    % for iEcho = 1 : nEchoes
    %     
    %
    %
    %     
    %     Mag.img = ( Img{ iEcho, 1, 1 }.img + Img{ iEcho, 1, 2 }.img )/2 ;
    %     Mag.img = Mag.img./max(Mag.img(:)) ;
    %     
    %
    %
    %
    %
    %     PhaseTE     = Phase{1}.copy() ;
    %     % extract iEcho
    %     PhaseTE.img = PhaseTE.img(:,:,:,iEcho) ;
    %     
    %     MagTE     = Mag{1}.copy() ;
    %
    %
    % Phase.img = angle( Img{ iEcho, 1, 2 }.img .*exp(i*Img{ iEcho, 2, 2 }.img) ...
    %     .* conj( Img{ iEcho, 1, 1 }.img .*exp(i*Img{ iEcho, 2, 1 }.img) ) ) ;
    %
    %     
    %     % if Params.nEchoes == 1 % 2 echoes acquired but only 1 input (phase *difference* image)
    %     %     assert( myisfield( Img{ 1, 2, iCurrent }.Hdr.MrProt, 'alTE') && numel( Img{ 1, 2, iCurrent }.Hdr.MrProt.alTE ) >= 2, ...
    %     %        'Expected Siemens FM phase difference image with at least 2 TEs specified in the DICOM header.' )
    %     %     % replace echo time stored in siemens header with the echo time difference between the 1st 2 echoes:       
    %     %
    %     %     Img{ 1, 2, iCurrent}.Hdr.EchoTime = ... 
    %     %         ( Img{ 1, 2, iCurrent }.Hdr.MrProt.alTE(2) - Img{ 1, 2, iCurrent }.Hdr.MrProt.alTE(1) )/1000 ; % [units : ms]
    %     %
    %     % end 
    %
    %
    %
    % for iEcho = 1 : Params.nEchoes
    %     
    %     for iCurrent = 1  : Params.nCurrents
    %         % -------
    %         % PROCESS GRE FIELD MAPS
    %         Img{ iEcho, 1, iCurrent } = MaRdI( Params.dataLoadDirectories{ iEcho, 1, iCurrent, iChannel }  ) ; % mag
    %         Img{ iEcho, 2, iCurrent } = MaRdI( Params.dataLoadDirectories{ iEcho, 2, iCurrent, iChannel }  ) ; % phase
    %     end
    %     
    %     if Params.nEchoes == 1 % 2 echoes acquired but only 1 input (phase *difference* image)
    %         assert( myisfield( Img{ 1, 2, iCurrent }.Hdr.MrProt, 'alTE') && numel( Img{ 1, 2, iCurrent }.Hdr.MrProt.alTE ) >= 2, ...
    %            'Expected Siemens FM phase difference image with at least 2 TEs specified in the DICOM header.' )
    %         % replace echo time stored in siemens header with the echo time difference between the 1st 2 echoes:       
    %
    %         Img{ 1, 2, iCurrent}.Hdr.EchoTime = ... 
    %             ( Img{ 1, 2, iCurrent }.Hdr.MrProt.alTE(2) - Img{ 1, 2, iCurrent }.Hdr.MrProt.alTE(1) )/1000 ; % [units : ms]
    %     
    %     end 
    %
    %     TE( iEcho ) = Img{ iEcho, 2, 1 }.Hdr.EchoTime ;
    %
    %     % ASSUMES :
    %     %   NO TX FREQUENCY SHIFT BETWEEN ACQUISITIONS
    %     % TODO
    %     %   adjust Field in the event Tx Freq. changed between scans/shim currents
    %     Mag     = Img{ iEcho, 1, 2 }.copy() ;
    %     Phase   = Img{ iEcho, 2, 2 }.copy() ;
    %     
    %     Mag.img = ( Img{ iEcho, 1, 1 }.img + Img{ iEcho, 1, 2 }.img )/2 ;
    %     Mag.img = Mag.img./max(Mag.img(:)) ;
    %     
    %     Phase.img = angle( Img{ iEcho, 1, 2 }.img .*exp(i*Img{ iEcho, 2, 2 }.img) ...
    %         .* conj( Img{ iEcho, 1, 1 }.img .*exp(i*Img{ iEcho, 2, 1 }.img) ) ) ;
    %     
    %     Phase.Hdr.MaskingImage = Mag.img > Params.threshold ;
    %
    %     Phase.unwrapphase( Mag, Params ) ; 
    %     Field = Phase.copy();
    %     Field.scalephasetofrequency() ;
    %
    %     % convert to Hz/A
    %     Field.img = Field.img/( Params.currents( iChannel, 2) - Params.currents( iChannel, 1) ) ;
    %     
    %     % smoothing: [3x3x3] magnitude-weighted Gaussian filtering
    %     Field.filter( Mag.img ) ; 
    %     
    %     fieldMaps( :,:,:, iEcho ) = Field.img ;
    %
    %     % mask out later echoes from field estimate if later echoes deviate
    %     % more than 50 Hz from the field map based on the early echoes:
    %     masks( :,:,:, iEcho, iChannel ) = Field.Hdr.MaskingImage & ( fieldMaps( :,:,:, iEcho ) ~= 0 ) & ( abs( fieldMaps(:,:,:, iEcho ) - fieldMaps(:,:,:, 1) ) < 50 ) ;
    %     
    % end 
    %
    % fieldMaps( ~masks(:,:,:,:,iChannel) ) = NaN ; 
    % dBdI(:,:,:, iChannel)  = median( fieldMaps, 4, 'omitnan' ) ;

end 

snr = mean(abs(Fields.img),4)./std(Fields.img,[],4) ;
snr = mean(abs(Mag.img),4)./std(Mag.img,[],4) ;

dBdI( ( isnan(dBdI) | isinf(dBdI) ) ) = 0 ;

Hdr = Field.Hdr ;
Hdr.MaskingImage = sum( sum( masks, 4 )>0, 5 ) == Params.nChannels ; % spatial support

img = dBdI .* repmat( Hdr.MaskingImage, [1 1 1 Params.nChannels] ) ;




disp( 'Forming interpolant...' )
disp( '(Computation time depends on input image size. This may take a few minutes.)' ) ;

% The following avoids the error from scatteredInterpolant when one
% attempts to form a 3d interpolant from a 2d input: 
isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
r0          = [X(:) Y(:) Z(:)] ;

Interpolant = scatteredInterpolant() ;
Interpolant.Points = r0(:, isValidDim0) ;



end
% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Static=true, Hidden=true)
% =========================================================================
function [ Params ] = assigndefaultparameters( Params, Specs )
%ASSIGNDEFAULTPARAMETERS  
% 
% Params = ASSIGNDEFAULTPARAMETERS( Params )
% 
% Add default parameters fields to Params without replacing values (unless empty)
%
% DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
% DEFAULT_INSTITUTIONNAME         = 'IUGM' ;
% DEFAULT_STATIONNAME             = 'MRC35049' ;

DEFAULTS.pathToShimReferenceMaps = [] ;
DEFAULTS.InstitutionName         = 'IUGM' ;
DEFAULTS.StationName             = 'MRC35049' ;

if ( nargin == 0 ) || isempty( Params )
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

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
function [ Field, Params ] = parseinput( Inputs )
%PARSEINPUT
% 
% Simple parser returns the optional user inputs Field and Params irrespective
% of their input order (convenient).
%
% [ Field, Params ] = PARSEINPUT( Inputs )

Field  = [] ;
Params = [] ;

nArgin = length( Inputs ) ;

if (nArgin > 0)
    if (nArgin <= 2)
        for iArg = 1 : nArgin
            switch class( Inputs{iArg} ) 
                case 'struct'
                    Params = Inputs{iArg} ;
                case 'FieldEval'
                    Field = Inputs{iArg} ;
            end
        end
    else
        error('Too many input arguments. Should be <=2. See help ShimOpt') ;
    end
end

if isempty( Params ) 
    Params.dummy = [] ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true)
% =========================================================================
function [ img, Hdr, Interpolant ] = loadshimreferencemaps( pathToShimReferenceMaps )
%LOADSHIMREFERENCEMAPS
%
% [ img, Hdr, Interpolant ] = LOADSHIMREFERENCEMAPS( filename )

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

if myisfieldfilled( RefMaps, 'Interpolant' ) 
    Interpolant = RefMaps.Interpolant ;
else
    Interpolant = [] ;
end

end
% =========================================================================

% =========================================================================
% =========================================================================
end

end
