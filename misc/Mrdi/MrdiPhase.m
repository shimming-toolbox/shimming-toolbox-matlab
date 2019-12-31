classdef MrdiPhase < MrdiInfo 
%MrdiPhase  MR DICOM Image Phase 
%
% Member methods for phase images 
% 
% e.g.
%   unwrapphase()
%   ...etc.
% 
% For documentation, type doc MrdiPhase
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Phase = MrdiPhase( Img )
   
    if ( nargin < 1 ) || isempty( Img )
        return ;
    elseif ( nargin ~= 1 ) || ~isa( Img, 'MaRdI' ) || ~Img.isphase()
        error('Constructor requires 1 input argument: a Mardi object initialized with phase data. See doc Mardi.') ;
    end

    Phase.copyproperties( Img ) ;   

end    
% =========================================================================
function [] = unwrapphase( Phase, varargin )
%UNWRAPPHASE
%
% Interface to SUNWRAP, to FSL Prelude, or to Abdul-Rahman's path-based phase unwrapper
%
% .......
%   
% Usage
%
%   [] = UNWRAPPHASE( Phase )
%   [] = UNWRAPPHASE( Phase, Mag )
%   [] = UNWRAPPHASE( Phase, Mag, Options )
%   
%   Phase and Mag are objects of type Mardi.
% 
%   Options is a struct that can contain the following fields:
%
%   .unwrapper 
%       == 'AbdulRahman_2007' [default if Phase.img is 3D] : 
%           calls UNWRAP3D. 
%           Not permitted if Phase.img is 2D (will default to SUNWRAP)
%           If Mag is supplied, Phase.Hdr.MaskingImage must
%           See HELP UNWRAP3D for description of permitted Options 
%
%       == 'Sunwrap' [default if Phase.img is 2D] : 
%           calls SUNWRAP. See HELP SUNWRAP for more details.
%
%       == 'FslPrelude' : calls PRELUDE. 
%           See HELP PRELUDE for description of permitted Options 
% 
%   .threshold  [default = 0.01] 
%       Relative threshold of Mag.img used to define the unwrapping region in
%       the SUNWRAP case and for the other 2 cases when Phase.Hdr.MaskingImage
%       is undefined (Mag.getreliabilitymask() is called).
%
%   NOTE: UNWRAP3D and PRELUDE support an Options.mask input to which
%   Phase.Hdr.MaskingImage will always be assigned. To assign the mask
%   manually, run Phase.setmaskingimage before calling Phase.unwrapphase. 


%% ------
% check inputs, assign defaults:
assert( strcmp( Phase.Hdr.PixelComponentPhysicalUnits, '0000H' ), 'SCALEPHASE2RAD() before unwrapping.' ) ;

if nargin > 1
    for iArg = 1 : length( varargin )
        switch class( varargin{iArg} )
            case 'Mardi'
                Mag = varargin{iArg} ;
                assert( Mag.ismagnitude(), 'See HELP Mardi.unwrapphase' ) ;
            case 'struct'
                Options = varargin{iArg} ;
            otherwise
                error('Unrecognized input. See HELP Mardi.unwrapphase') ;
        end
    end
end

DEFAULTS.threshold     = 0.01 ;
DEFAULTS.isLinearizing = false ;

if ( size( Phase.img, 3 ) > 1 )
    is3d = true ;
    DEFAULTS.unwrapper = 'AbdulRahman_2007' ;    
else
    is3d = false ;
    DEFAULTS.unwrapper = 'Sunwrap' ;
end

Options.dummy = [];
Options       = assignifempty( Options, DEFAULTS ) ;

if ~any( strcmp( Options.unwrapper, {'Sunwrap','sunwrap','FslPrelude','AbdulRahman_2007'} ) )
    error('Unrecognized "Options.unwrapper" input. See HELP Mardi.unwrapphase') ;
elseif strcmp( Options.unwrapper, 'AbdulRahman_2007' ) && ( size( Phase.img, 3 ) == 1 )
    warning('Options.unwrapper = AbdulRahman_2007 is incompatible with 2d images. Using Sunwrap method.') ;
    Options.unwrapper = 'Sunwrap' ;
end

%% ------
nVolumes = Phase.getnumberofmeasurements() ;
nEchoes  = size( Phase.img(), 4 ) ;

if ~exist( 'Mag' )
    assert( strcmp( Options.unwrapper, 'AbdulRahman_2007' ), ...
        ['Missed required input: corresponding Magnitude Mardi-object. ' ...
        ' (Applies to Sunwrap and FslPrelude cases.) See HELP Mardi.unwrapphase'] ) ;

    assert( myisfieldfilled( Phase.Hdr, 'MaskingImage' ), ...
        ['Either a logical masking array must be defined in Phase.Hdr.MaskingImage ' ...
         'OR the input must include the corresponding Magnitude image. See HELP Mardi.unwrapphase'] ) ;

else
    assert( Phase.iscoincident( Mag ), ['Inputs Mag.img and Phase.img must correspond '...
        '(respective voxel positions and number of measurements should be identical'] )

    if ~myisfieldfilled( Phase.Hdr, 'MaskingImage' )
        mask = Mag.getreliabilitymask( Options.threshold ) ; 
        % indexing for potential bug, e.g. in Siemens dual echo field mapping, where nEchoes magnitude input > nEchoes phase input...
        % TODO: find a more elegant way of handling this...
        Phase.setmaskingimage( mask(:,:,:,1:nEchoes,1:nVolumes) ) ; 
    end
end

if ( nEchoes > 1 ) && Options.isLinearizing
    % phase evolution between 1st 2 echoes, to be used to correct 2*pi offsets
    % of individual echoes after unwrapping
    PhaseDiff      = Phase.copy() ;

    img            = Mag.img(:,:,:,1,:) .* exp( i*Phase.img(:,:,:,1,:) ) ;
    img(:,:,:,2,:) = Mag.img(:,:,:,2,:) .* exp( i*Phase.img(:,:,:,2,:) ) ;

    PhaseDiff.img  = angle( img(:,:,:,2,:) .* conj(img(:,:,:,1,:) ) ) ;

    PhaseDiff.Hdr.EchoTime = Phase.getechotime(2) - Phase.getechotime(1) ; % [units : ms]

    mask = Phase.Hdr.MaskingImage(:,:,:,1,:) & Phase.Hdr.MaskingImage(:,:,:,2,:) ;
    PhaseDiff.setmaskingimage( logical( mask ) & ~isnan( PhaseDiff.img ) ) ;

    PhaseDiff.unwrapphase( Mag, Options ) ;

    % if PhaseDiff is not centered around 0, shift it:
    mask      = PhaseDiff.Hdr.MaskingImage(:,:,:,1,:) ;
    phaseDiff = PhaseDiff.img(:,:,:,1,:) ; 
    n         = round( median( phaseDiff(mask) )/pi ) ;
    PhaseDiff.img(:,:,:,1,:) = PhaseDiff.img(:,:,:,1,:) - n*pi ;
end    


%% ------
for iVolume = 1 : nVolumes

    display(['Unwrapping volume ' num2str(iVolume) ' of ' num2str(nVolumes) ] ) ;    

    for iEcho = 1 : nEchoes
        
        if nEchoes > 1
            display(['Unwrapping echo ' num2str(iEcho) ' of ' num2str(nEchoes) ] ) ;    
        end

        switch Options.unwrapper

            case 'AbdulRahman_2007'
                
                Phase.img(:,:,:,iEcho, iVolume) = unwrap3d( Phase.img(:,:,:,iEcho, iVolume), logical(Phase.Hdr.MaskingImage(:,:,:,iEcho, iVolume)), Options ) ;

            case 'FslPrelude'

                if myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage)
                    Options.mask = single( Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume) ) ;
                end
                
                Options.voxelSize = Phase.getvoxelspacing() ;   

                Phase.img(:,:,:,iEcho, iVolume) = prelude( Phase.img(:,:,:,iEcho, iVolume), Mag.img(:,:,:,iEcho, iVolume), Options ) ;

            case {'Sunwrap', 'sunwrap'}
                
                iMag      = Mag.img(:,:,:,iEcho,iVolume) ;
                iMag      = iMag./max(iMag(:)) ;
                Phase.img(:,:,:,iEcho, iVolume) = sunwrap( iMag .* exp( 1i* Phase.img(:,:,:,iEcho,iVolume) ), Options.threshold ) ;

        end

        if Options.isLinearizing && ( iEcho == 1 ) 
            % if phase(TE1) is not centered around 0, shift it:
            mask  = Phase.Hdr.MaskingImage(:,:,:,1,iVolume) ;
            phase = Phase.img(:,:,:,iEcho,iVolume) ; 
            n     = fix( median( phase(mask) )/pi ) ;
            if n ~= 0
                display( ['Recentering phase of 1st echo to sit between [-pi,pi] (global subtraction of ' num2str(n) 'pi)'] ) 
                Phase.img(:,:,:,iEcho,iVolume) = Phase.img(:,:,:,iEcho,iVolume) - n*pi*double(mask) ;
            end
        end

    end

    if ( nEchoes > 1 ) && Options.isLinearizing
    %% -----
    % Correct inter-echo wraps  
        for iEcho = 1 : nEchoes 
            % (assuming phase(TE=0) offset = 0):
            phaseEstimate = (Phase.getechotime(iEcho)/PhaseDiff.Hdr.EchoTime)*PhaseDiff.img ;

            err = phaseEstimate - Phase.img(:,:,:,iEcho,iVolume) ;
            err = median( err( Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume ) ) ) ; 

            % When the median deviation between estimate and unwrapped
            % measurement exceeds pi, correct the measurement by adding the
            % appropriate pi-multiple:
            n  = ( abs(err) > pi ) .* round( err/pi )  ;

            Phase.img(:,:,:,iEcho,iVolume) = Phase.img(:,:,:,iEcho,iVolume) + n*pi*Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume ) ;
        end
    end

end

Phase.img = double( Phase.img ) ;

%% -----
% if time series, correct for potential wraps between time points:
if nVolumes > 1
    display('Correcting for potential temporal phase wraps...')
    
    for iEcho = 1 : nEchoes
        % correct temporal wraps by comparing to phaseEstimate:
        phaseEstimate = median( Phase.img(:,:,:,iEcho,:),  5, 'omitnan' ) ;

        % normalize the estimate so its spatial median is within [-pi,pi]
        tmpMask = logical( prod( Phase.Hdr.MaskingImage(:,:,:,iEcho,:), 5 ) ) ;
        n       = fix( median( phaseEstimate(tmpMask) )/pi ) ;
        phaseEstimate = phaseEstimate - n*pi ;
        
        % Wherever the absolute deviation from the estimate exceeds pi,
        % correct the measurement by adding the appropriate pi-multiple:
        for iVolume = 1 : nVolumes
            dPhase = phaseEstimate - Phase.img( :,:,:,iEcho,iVolume )  ;
            n      = ( abs(dPhase) > pi ) .* fix( dPhase/pi ) ;
            Phase.img( :,:,:,iEcho,iVolume ) = Phase.img(:,:,:,iEcho,iVolume ) + n*pi ;
        end
    end
end

%% -----
% update header
Img.Hdr.ImageType         = 'DERIVED\SECONDARY\P\' ; 
Img.Hdr.SeriesDescription = ['phase_unwrapped_' Options.unwrapper ] ; 

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=protected)
% =========================================================================
function [] = scalephasetofrequency( Img, undoFlag )
%SCALEPHASETOFREQUENCY
%
% Converts unwrapped phase [units:rad] to field [units: Hz]
% 
%   Field = scalephasetofrequency( UnwrappedPhase )
%   Phase = scalephasetofrequency( Field, -1 )
%   
%   The 'undo' mode with -1 as the 2nd argument scales from Hz back to rad

assert( ~Img.ismagnitude(), 'Input cannot be a magnitude image' ) ;

te = Img.getechotime()*(1E-3) ;

if (nargin < 2) || (undoFlag ~= -1)
    assert( Img.isphase() && strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0000H' ), ...
        'Expected Phase Img input with voxel values in radians' ) ;
    
    scalingFactors = 1./(2*pi*te) ;
    Img.Hdr.PixelComponentPhysicalUnits = '0005H' ; % i.e. Hz

elseif (undoFlag == -1)
    assert( strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0005H' ), ...
        'Expected Field Img input with voxel values in Hz' )
    
    scalingFactors = 2*pi*te ;
    Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none
end

for iEcho = 1 : numel(te) 
    Img.img(:,:,:,iEcho,:) = scalingFactors(iEcho) * Img.img(:,:,:,iEcho,:) ;
end

end
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
