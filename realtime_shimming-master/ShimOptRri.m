classdef ShimOptRri < ShimOpt
%SHIMOPT - Shim Optimization
%
% .......
% 
% Usage
%
% Shim = ShimOpt( )
% Shim = ShimOpt( Params )
% 
% Defaults 
% 
% Params.pathToShimReferenceMaps =
%   '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007.mat'
%
% Params.ProbeSpecs = [] ;
%   .ProbeSpecs is a parameters struct for ProbeTracking(). 
%   See HELP ProbeTracking() for more information.
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
%               For realtime shimming, relates vector relating field to pressure
%               [units Hz/Pa]
%           .dcCurrentsOffsets
%               For realtime shimming, vector of "y-intercept" currents 
%               (i.e. currents for pressure = 0 Pa)
%               [units A]
%
%       .Probe
%           Object of type ProbeTracking
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    ProbeTracking
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
% Updated::20161122::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
%
% .....
% .Probe
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
% OPTIMIZESHIMCURRENTS()
%   Run (unconstrained) Cg solver;
%   If solution achievable given system constraints, return solution;
%   Else, run fmincon given constraints & return that solution instead;
%
% .....
% EXTENDHARMONICFIELD()
%   Write function to check if field map is exists + is reasonable over
%   the shim voi. if some portion is missing, a harmonic extension could
%   be performed to fill in some of the missing field.
% 
% .....
% CALIBRATEREALTIMEUPDATES()
%   Routine to extract relevant pressure measurements, process field maps,
%   set model parameters, etc.
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
    Probe ; % object of type ProbeTracking
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( Params )
%SHIMOPT - Shim Optimization

DEFAULT_PATHTOSHIMREFERENCEMAPS = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007.mat';
DEFAULT_PROBESPECS = [] ;

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'ProbeSpecs' ) || isempty(Params.ProbeSpecs)
   Params.ProbeSpecs = DEFAULT_PROBESPECS ;
end

ShimUse.display(['\n Preparing for shim ...  \n\n'...
        'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;

load( Params.pathToShimReferenceMaps ) ;

%%-----
% dB/dI linear 'Current-to-Field' operator
Shim.img              = SpineShim.img ;
Shim.Hdr              = SpineShim.Hdr ;

if myisfield( SpineShim, 'mask' )
    Shim.Hdr.MaskingImage = SpineShim.mask ;
end

Shim.Field = [ ] ;       
Shim.Model = [ ] ; 
Shim.Probe = ProbeTracking( Params.ProbeSpecs )  ; 

end
% =========================================================================
function [] = delete( Shim )
%DELETE  
% DELETE( Shim )
% 
% Destructor. Calls Probe.deletecomport( ) 

if ~isempty( Shim.Probe )
    Shim.Probe.delete();
end

clear Shim ;

end
% =========================================================================
function Shim = calibraterealtimeupdates( Shim, Params )
%CALIBRATEREALTIMEUPDATES
% 
% CALIBRATEREALTIMEUPDATES asks user to select the median pressure from the
% pressure logs corresponding to the inspired & expired field maps. From
% these, and the associated optimal currents for the 2 respiratory states,
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
%           .pressureLog
%
%       .Expired
%           .currents
%           .pressureLog
        
ShimUse.display( ['\n ------- \n ' ...
    'Determine median pressure over inspired apnea : \n \n '] ) ;

Params.Inspired.medianPressure = ...
    Shim.Probe.userselectmedianpressure( Params.Inspired.pressureLog ) ; 
ShimUse.display( ...
    ['Median pressure : ' num2str( Params.Inspired.medianPressure )] ) ;

ShimUse.display( ['\n ------- \n ' ...
    'Determine median pressure over expired apnea : \n \n '] ) ;

Params.Expired.medianPressure = ...
    Shim.Probe.userselectmedianpressure( Params.Expired.pressureLog ) ; 
ShimUse.display( ...
    ['Median pressure : ' num2str( Params.Expired.medianPressure )] ) ;

Shim.setcouplingcoefficients( ...
    Params.Inspired.currents, Params.Expired.currents, ...
    Params.Inspired.medianPressure, Params.Expired.medianPressure ) ;

Shim.setdccurrentoffsets( ...
    Params.Inspired.currents, Params.Expired.currents, ...
    Params.Inspired.medianPressure, Params.Expired.medianPressure ) ;

ShimUse.display( ['\n ------- \n ' ...
    'Optimal DC current offsets (in amperes): ' num2str(Shim.Model.dcCurrentOffsets') '\n \n '] ) ;

Shim.setupdateoperator() ;

end
% =========================================================================
function Shim = setcouplingcoefficients( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pressureInspired, pressureExpired )
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
    A*(currentsInspired - currentsExpired)/(pressureInspired - pressureExpired) ;

end
% =========================================================================
function Shim = setdccurrentoffsets( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pressureInspired, pressureExpired )
%SETDCCURRENTOFFSETS
% 
% Compute and set optimal shim DC current offsets (bias)
%
% Shim = SETDCCURRENTOFFSET( Shim, iIn, iEx, pIn, pEx  )
%
% Sets Shim.Model.dcCurrentsOffsets

Specs = ShimSpecs();

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

shimFieldOffset = ...
    A*(currentsExpired*pressureInspired - currentsInspired*pressureExpired)/(pressureInspired - pressureExpired);

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

[X,Y,Z] = Field.getvoxelpositions ;
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
function Shim = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Shim = OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]

Specs = ShimSpecs();

if nargin < 1
    error('Function requires at least 1 argument of type ShimOpt')
elseif nargin == 1
    Params.dummy = [];
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end

% Params for conjugate-gradient optimization
CgParams.tolerance     = 1E-6 ;
CgParams.maxIterations = 10000 ;    

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

b = M*(-Shim.Field.img(:)) ;

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*b, ... % effective solution vector
                            zeros( [Specs.Amp.nActiveChannels 1] ), ... % initial model (currents) guess
                            CgParams ) ;
% ------- 
% TODO: CHECK SOL. VECTOR FOR CONDITIONS 
% + allow more optional params for optimization 
isCurrentSolutionOk = false ;

if ~isCurrentSolutionOk
    
    [X0,X1,X2,X3] = ShimCom.getchanneltobankmatrices( ) ;

    Options = optimset(...
        'DerivativeCheck','off',...
        'GradObj','on',...
        'Display', 'off',... %'iter-detailed',...
        'MaxFunEvals',36000,...
        'TolX',1e-11,...
        'TolCon',1E-8);

    A=M*A;
    tic
    [Shim.Model.currents] = fmincon( ...
        @shimcost,...
        ones(Specs.Amp.nActiveChannels,1),...
        [],...
        [],...
        [],...
        [],...
        -Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
        Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
        @first_order_norm,...
        Options);
    toc
    
end

function [f, df] = shimcost( currents )
    
     y=A*currents - b;
     f=y'*y;
     df=2*A'*y;
    
end

function [C, Ceq] = first_order_norm( currents )
% C(x) <= 0
% (e.g. x = currents)
%
    C   = 0; 
    Ceq = [];
    waterLevel = 1E-8;
    % split currents up into banks
    % -------
    % bank 0
    i0 = X0*currents;

    % -------
    % bank 1
    i1 = X1*currents;

    % -------
    % bank 2
    i2 = X2*currents;

    % -------
    % bank 3
    i3 = X3*currents;

    % Overall abs current cannot exceed Specs.Amp.maxCurrentPerBank (e.g. 20 A)
    C(1) = sum( abs(i0) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(2) = sum( abs(i1) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(3) = sum( abs(i2) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(4) = sum( abs(i3) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 

    % pos. current cannot exceed Specs.Amp.maxCurrentPerRail (e.g. + 10 A) 
    C(5) = abs(sum( ((i0>0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(6) = abs(sum( ((i1>0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(7) = abs(sum( ((i2>0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(8) = abs(sum( ((i3>0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    
    % neg. current cannot be below Specs.Amp.maxCurrentPerRail (e.g. - 10 A) 
    C(9)  = abs(sum( ((i0<0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(10) = abs(sum( ((i1<0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(11) = abs(sum( ((i2<0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(12) = abs(sum( ((i3<0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 

end
    
Shim = Shim.setforwardmodelfield ;

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
%   where UO * pressure = currentsUpdate

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
%    .......................
%   
%   The following Params.fields are supported
%
%       .maxAbsField 
%           maximum absolute voxel value assumed to represent an accurate
%           field measurement. Voxels with abs-values greater than this
%           might stem from errors in the unwrapping.
%           default: 500 Hz
%
%       .maxFieldDifference
%           maximum absolute voxel-wise difference assumed to be valid between 
%           Field1 & Field2 (e.g. 'inspired field' vs. 'expired field')
%           default: 150 Hz (See Verma T, Magn Reson Med, 2014)
%
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
% Truncation Shim.Field.Hdr.MaskingImage) operator (e.g. M*b, 'picks out' the
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
function currents = computerealtimeupdate( Shim, Params )
% COMPUTEREALTIMEUPDATE
%
% currents = COMPUTEREALTIMEUPDATE( Shim, pressures, Params ) ;

if Params.extrapolationOrder == 0

    currents = Shim.Model.dcCurrentOffsets + ...
        Shim.Model.updateOperator * Shim.Probe.Data.pressure(end) ; 

% else
%     currents = Shim.Model.dcCurrentOffsets + Shim.Model.updateOperator * ...
%         predictpressure( Shim.Probe.Data.pressure, Params.extrapolationDelay,  )
%    % do not pass by value all of the pressure reading... (unnecesary copying increases with duration of experiment) 

end        

% function [predictedPressure] = predictpressure( pressures, delay, period )
% % PREDICTPRESSURE 
%
% if Params.extrapolationOrder == 1
%
%     predictedPressure = pressures(end) + (Params.extrapolationDelay/Shim.Opt.Probe.Specs.arduinoPeriod/1000)* ...
%             (pressure(end) - pressure(end-1)) );
%
% elseif Params.extrapolationOrder == 2
%
%     % predictedPressure = pressures(end) + (Params.extrapolationDelay/Shim.Opt.Probe.Specs.arduinoPeriod/1000)* ...
%     %         (pressure(end) - pressure(end-1)) );
% end

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
methods(Static)

% =========================================================================
function [Field, Extras] = mapfield( Params )
% MAPFIELD
%
% Field = MAPFIELD( Params ) 
%
% Params must contain the following fields
%
%   .Path.Mag.echo1
%   .Path.Mag.echo2
%   .Path.Phase.echo1
%   .Path.Phase.echo2
%   .threshold (as a fraction (<1) of max measured magnitude intensity)
%       determines the phase region to be unwrapped (i.e. areas of low signal
%       are ignored)
%   .isFilteringField 

DEFAULT_ISCORRECTINGPHASEOFFSET        = true ;
DEFAULT_ISUNWRAPPINGECHOESINDIVIDUALLY = false ; 
DEFAULT_ISFILTERINGFIELD               = false ;

if ~myisfield( Params, 'isCorrectingPhaseOffset' ) || isempty( Params.isCorrectingPhaseOffset ) 
    Params.isCorrectingPhaseOffset = DEFAULT_ISCORRECTINGPHASEOFFSET ;
end

if ~myisfield( Params, 'isUnwrappingEchoesIndividually' ) || isempty( Params.isUnwrappingEchoesIndividually ) 
    Params.isUnwrappingEchoesIndividually = DEFAULT_ISUNWRAPPINGECHOESINDIVIDUALLY ;
end

if ~myisfield( Params, 'isFilteringField' ) || isempty( Params.isFilteringField ) 
    Params.isFilteringField = DEFAULT_ISFILTERINGFIELD ;
end

Extras = [] ;

MagEcho1   = MaRdI(  Params.Path.Mag.echo1  ) ;
MagEcho2   = MaRdI(  Params.Path.Mag.echo2  ) ;

PhaseEcho1 = MaRdI(  Params.Path.Phase.echo1  ) ;
PhaseEcho2 = MaRdI(  Params.Path.Phase.echo2  ) ;

%-------
% Mask 
PhaseEcho1.Hdr.MaskingImage = ( MagEcho1.img ) > Params.threshold ;
PhaseEcho2.Hdr.MaskingImage = PhaseEcho1.Hdr.MaskingImage ;

if ~Params.isUnwrappingEchoesIndividually

    img            = MagEcho1.img .* exp(-i*PhaseEcho1.img) ;
    img(:,:,:,2)   = MagEcho2.img .* exp(-i*PhaseEcho2.img) ;
    PhaseEcho1.img = angle( img(:,:,:,2) ./ img(:,:,:,1) ) ;
    PhaseEcho1.unwrapphase( ) ;
    PhaseEcho1.Hdr.EchoTime = (PhaseEcho1.Hdr.EchoTime - PhaseEcho2.Hdr.EchoTime);
    Field = PhaseEcho1.scalephasetofrequency() ;


else
    % ------
    % Unwrap phase (Using Abdul-Rahman method)
    PhaseEcho1 = PhaseEcho1.unwrapphase( ) ;
    PhaseEcho2 = PhaseEcho2.unwrapphase( ) ;

    %-------
    % Compute field map
    mask = logical( PhaseEcho1.Hdr.MaskingImage .* PhaseEcho2.Hdr.MaskingImage ) ;

    Params.isCorrectingPhaseOffset =false;

    if Params.isCorrectingPhaseOffset    

        PhaseEcho2Plus2Pi              = PhaseEcho2.copy() ;
        PhaseEcho2Plus2Pi.img( mask )  = PhaseEcho2.img( mask ) + 2*pi ;

        PhaseEcho2Minus2Pi             = PhaseEcho2.copy() ;
        PhaseEcho2Minus2Pi.img( mask ) = PhaseEcho2.img( mask ) - 2*pi ;

        Field0        = PhaseEcho1.mapfrequencydifference( PhaseEcho2 ) ;
        FieldPlus2Pi  = PhaseEcho1.mapfrequencydifference( PhaseEcho2Plus2Pi ) ;
        FieldMinus2Pi = PhaseEcho1.mapfrequencydifference( PhaseEcho2Minus2Pi ) ;

        fieldNorms    = [0 0 0];
        fieldNorms(1) = norm( Field0.img( mask ) ) ;
        fieldNorms(2) = norm( FieldPlus2Pi.img( mask ) ) ;
        fieldNorms(3) = norm( FieldMinus2Pi.img( mask ) ) ;

        [~,iField] = min( fieldNorms ) ;

        if iField == 1

            Field = Field0;

        elseif iField == 2

            ShimUse.display('Correcting field offset: Adding 2pi to phase of 2nd echo.')
            PhaseEcho2 = PhaseEcho2Plus2Pi ;
            Field = FieldPlus2Pi;

        elseif iField == 3

            ShimUse.display('Correcting field offset: Subtracting 2pi from phase of 2nd echo.')
            PhaseEcho2 = PhaseEcho2Minus2Pi ;
            Field = FieldMinus2Pi;

        end
    end

    Field = Field0;
    Field.Hdr.MaskingImage = mask ;  
    Extras.PhaseEcho1 = PhaseEcho1 ;
    Extras.PhaseEcho2 = PhaseEcho2 ;

end

Params.isExtrapolatingField = true ;
Params.filterRadius = 3 ;

if Params.isFilteringField

    disp('Filtering field map...') 

    gridSizeImg = Field.getgridsize(); 
    Field = Field.zeropad( [3 3 3], 'both' ) ;
    
    mask = Field.Hdr.MaskingImage ;

    [Extras.LocalField, Field] = Field.extractharmonicfield( Params ) ;

    Field = Field.cropimg( gridSizeImg ) ;
    % dbstop in ShimOpt at 770   
    % mask = dilater( mask, [5 5 5] ) ;
    % mask(:,:,end) = [] ;    
    % Tmp.voxelSize            = Field.getvoxelsize() ;
    % Tmp.expansionOrder       = 1 ;
    % Tmp.radius               = 10 ;
    % Tmp.isDisplayingProgress = true ;
    % [xEp,A,M] = extendharmonicfield( Field.Hdr.MaskingImage .* Field.img, logical(mask), logical(Field.Hdr.MaskingImage), Tmp ) ;
    %
    % Field.img = (Field.Hdr.MaskingImage .* Field.img) + xEp ;
    % Field.Hdr.MaskingImage = (Field.img ~=0) ;
    %
    % NiiOptions.filename = './precrop'
    % nii(Field.img, NiiOptions) ;
    %
    % NiiOptions.filename = './postcrop'
    % nii(Field.img, NiiOptions) ;

end


end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

