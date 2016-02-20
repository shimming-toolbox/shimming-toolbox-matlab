classdef ShimOpt
%SHIMOPT
% 
% Shim Optimization
%
% =========================================================================
% *** TODO 
% .....
% Generally:
%   Lots of overlap between ShimOpt and MaRdI methods. 
%   Can ShimOpt be a MaRdI subclass? 
%
%   Rename Shims.mask as Shims.support ?
%       Shims.Operator as Shims.operator ?
% .....
% INTERPTOCALIBRATIONGRID()
%   Rename as INTERPTOREFERENCEGRID()
%   (check conflicts)
%
%% .....
% OPTIMIZESHIMCURRENTS()
%   Run (unconstrained) CG solver;
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
% extractharmonicfield() & optimizeshimcurrents()
%   
%   documentation + optional params
% .....
% RESLICEIMG()
%   griddata takes too f-ing long.
%   write interp function in cpp
%
% =========================================================================

properties
    Hdr ;
    Field ;
    Model ;
    Parameters ;
    Operator ;
    mask ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( pathToCalibrationInfo )
%SHIMOPT - Shim Optimization
%
% Shim = ShimOpt( pathToCalibrationInfo )
%
%   Shim contains fields
%
%       .Hdr
%           Info Re: calibration data
%
%       .Field
%           Object of type MaRdI pertaining to field distribution to be shimmed
%
%       .Model
%           .currents  (i)
%           .field     (Ai)
%
%       .Operator    
%           'Current-to-Field' matrix operator (A)
%
%       .Parameters
%
%       .mask
% .......
    fprintf(['\n Preparing for shim \n\n #####  \n\n'...
            'Loading SpineShim calibration info from ' pathToCalibrationInfo '\n\n']) ;

    load( pathToCalibrationInfo ) ;

    Shim.Hdr  = SpineShim.Hdr ;

    nVoxels = numel( SpineShim.img(:,:,:, 1) ) ;

    Shim.Parameters.nActiveChannels = size( SpineShim.img, 4 ) ;
    % Parameters for conjugate-gradient optimization
    Shim.Parameters.CG.tolerance     = 1E-6 ;
    Shim.Parameters.CG.maxIterations = 10000 ;    

    %%-----
    % dBdI linear 'Current-to-Field' operator
    Shim.Operator = SpineShim.img ;
    % Shim.Operator = zeros( nVoxels, Shim.Parameters.nActiveChannels ) ; 
    
    % for channel = 1 : Shim.Parameters.nActiveChannels
    %     Shim.Operator(:, channel) = reshape( SpineShim.img(:,:,:, channel), [nVoxels 1] ) ;
    % end

    Shim.mask = SpineShim.mask ;
    Shim.Field = [ ] ;       
    Shim.Model = [ ] ; 

end
% =========================================================================
function Shim = extractharmonicfield( Shim )
    %EXTRACTHARMONICFIELD
    % Extract (smooth) harmonic field via RESHARP (Sun, H. Magn Res Med, 2014)
    %
    % 
    [localField, reducedMask ] = resharp( Shim.Field.img, ...
                                     Shim.Field.Hdr.MaskingImage, ... 
                                     Shim.Field.getvoxelsize(), ...
                                     2*Shim.Field.getvoxelsize(), ...
                                     0) ;

    reducedMask                = shaver(reducedMask, 1) ;

    Shim.Field.img             = reducedMask .* ( Shim.Field.img - localField ) ;

    Shim.Field.Hdr.MaskingImage = reducedMask ;

end
% =========================================================================
function Shim = computecurrentsoffset( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pressureInspired, pressureExpired )
%COMPUTECURRENTSOFFSET
%
%   Syntax
%
%   Shim = COMPUTECURRENTSOFFSET( Shim, iIn, iEx, pIn, pEx  )
%
%       creates field Shim.Model.currentsOffset
% ------
Specs = ShimSpecs();

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

shimFieldOffset = ...
    A*(currentsInspired*pressureExpired - currentsExpired*pressureInspired)/(pressureInspired - pressureExpired);

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currentsOffset = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*M*shimFieldOffset, ... % effective solution vector
                            zeros( [Shim.Parameters.nActiveChannels 1] ), ... % initial model (currents) guess
                            Shim.Parameters.CG ) ;

end
% =========================================================================
function Shim = definecouplingcoefficients( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pressureInspired, pressureExpired )
%DEFINECOUPLINGCOEFFICIENTS
%
%   iIn & iEx are vectors of optimal currents for inspired and expired fields.
%   pIn & pEx the associated pressures (scalars)
%
%   Syntax
%
%   Shim = DEFINECOUPLINGCOEFFICIENTS( Shim, iIn, iEx, pIn, pEx )
%
%       creates field Shim.Model.couplingCoefficients
% ------
Specs = ShimSpecs();


A = Shim.getshimoperator ;

Shim.Model.couplingCoefficients = ...
    A*(currentsExpired - currentsInspired)/(pressureInspired - pressureExpired) ;

end
% =========================================================================
function [Shim] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
%   OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]
%
% ------
Specs = ShimSpecs();

if nargin < 1
    error('Function requires at least 1 argument of type ShimOpt')
elseif nargin == 1
    Params.dummy = [];
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end


A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

b = M*(-Shim.Field.img(:)) ;

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*b, ... % effective solution vector
                            zeros( [Shim.Parameters.nActiveChannels 1] ), ... % initial model (currents) guess
                            Shim.Parameters.CG ) ;

% ------- 
% check sol vector for conditions 
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
    [currents] = fmincon( ...
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

    Shim.Model.currents = currents ;
    Shim = Shim.forwardmodelfield ;
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

    % Overall abs current cannot exceed 20 A / bank 
    C(1) = sum( abs(i0) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(2) = sum( abs(i1) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(3) = sum( abs(i2) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(4) = sum( abs(i3) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 

    % pos. current cannot exceed +10 A / bank 
    C(5) = abs(sum( ((i0>0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(6) = abs(sum( ((i1>0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(7) = abs(sum( ((i2>0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(8) = abs(sum( ((i3>0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    
    % neg. current cannot be below -10 A / bank 
    C(9)  = abs(sum( ((i0<0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(10) = abs(sum( ((i1<0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(11) = abs(sum( ((i2<0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(12) = abs(sum( ((i3<0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 

end

end
% =========================================================================
function nVoxels = getnumberofvoxels( Shim )
    nVoxels = prod( Shim.getgridsize( ) ) ;

end
% =========================================================================
function fieldOfView = getfieldofview( Shim )
    fieldOfView = [ Shim.Hdr.PixelSpacing(1) * double( Shim.Hdr.Rows - 1 ), ...
                    Shim.Hdr.PixelSpacing(2) * double( Shim.Hdr.Columns - 1 ), ...
                    Shim.Hdr.SpacingBetweenSlices * double( Shim.Hdr.NumberOfSlices - 1) ] ;

end
% =========================================================================
function voxelSize = getvoxelsize( Shim )
    voxelSize = [ Shim.Hdr.PixelSpacing(1) Shim.Hdr.PixelSpacing(2) Shim.Hdr.SpacingBetweenSlices ] ;

end
% =========================================================================
function gridSize = getgridsize( Shim )
    gridSize = [ Shim.Hdr.Rows, Shim.Hdr.Columns, Shim.Hdr.NumberOfSlices ] ;

end
% =========================================================================
function [rHat,cHat,sHat] = getdirectionalcosines( Shim ) 
% GETDIRECTIONALCOSINES
% 
% [r,c,s] = GETDIRECTIONALCOSINES( Shim ) 
%   r: directional cosine of image rows relative to dicom reference coordinate system axes
%   c: " " of image columns relative to " "
%   s: " " of image slices relatice to " "
      
    rHat = Shim.Hdr.ImageOrientationPatient(1:3) ;  % unit vector row dir
    cHat = Shim.Hdr.ImageOrientationPatient(4:6) ; % unit vector column dir
    sHat = cross( rHat, cHat ) ;  % unit vector slice direction
    
    % in general, the coordinate is increasing along slice dimension with slice index:
    [~,iS] = max( abs(sHat) ) ;
    if sHat(iS) < 0
        sHat = -sHat ;
    end

end 
% =========================================================================
function [X,Y,Z] = getvoxelpositions( Shim )
% GETVOXELPOSITIONS
% 
% [X,Y,Z] = GETVOXELPOSITIONS( Img ) 
%
%   Returns three 3D arrays of doubles, each element containing the
%   location [units: mm] of the corresponding voxel with respect to 
%   DICOM's 'Reference Coordinate System'.

    [rHat, cHat, sHat] = Shim.getdirectionalcosines() ; 
    
    % form DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
    %
    % If Anatomical Orientation Type (0010,2210) is absent or has a value of
    % BIPED, the x-axis is increasing to the left hand side of the patient. The
    % y-axis is increasing to the posterior side of the patient. The z-axis is
    % increasing toward the head of the patient.
    %
    % Arrays containing row, column, and slice indices of each voxel
    assert( ~myisfield( Shim.Hdr, 'AnatomicalOrientationType' ) || ...
            strcmp( Shim.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
            'Error: AnatomicalOrientationType not supported.' ) ;
            
    [R,C,S] = ndgrid( [0:1:Shim.Hdr.Rows-1], ...
                      [0:1:Shim.Hdr.Columns-1], ...
                      [0:1:(Shim.Hdr.NumberOfSlices-1)] ) ; 

    %-------
    % SCALE to physical by sample spacing 
    % (i.e. grid size, i.e. effective voxel size) 
    voxelSize = Shim.getvoxelsize() ;
    
    R = voxelSize(1)*double(R);
    C = voxelSize(2)*double(C);
    S = voxelSize(3)*double(S);

    %-------
    % ROTATE to align row direction with x-axis, 
    % column direction with y-axis, slice with z-axis
    X1 = cHat(1)*R + rHat(1)*C + sHat(1)*S;
    Y1 = cHat(2)*R + rHat(2)*C + sHat(2)*S;
    Z1 = cHat(3)*R + rHat(3)*C + sHat(3)*S;

    %-------
    % TRANSLATE w.r.t. origin 
    % (i.e. location of 1st element: .img(1,1,1))
    X = Shim.Hdr.ImagePositionPatient(1) + X1 ; 
    Y = Shim.Hdr.ImagePositionPatient(2) + Y1 ; 
    Z = Shim.Hdr.ImagePositionPatient(3) + Z1 ; 

end
% =========================================================================
function Shim = resliceimg( Shim, X_1, Y_1, Z_1, interpolationMethod ) 
%RESLICEIMG
%
%   Shim = RESLICEIMG( Shim, X, Y, Z )
%   Shim = RESLICEIMG( Shim, X, Y, Z, interpolationMethod ) 
%
%   X, Y, Z MUST refer to X, Y, Z patient coordinates (i.e. of the DICOM
%   reference coordinate system)
%   
%   Optional interpolationMethod is a string supported by griddata().
%   See: help griddata  
%

    %%------ 
    % Reslice to new resolution
    if nargin < 5
        interpolationMethod = 'linear' ;
    end
    
    [X_0, Y_0, Z_0] = Shim.getvoxelpositions( ) ;
  
    Tmp = zeros( [size(X_1) Shim.Parameters.nActiveChannels] ) ;

    for iChannel = 1 : Shim.Parameters.nActiveChannels    
        disp( iChannel ) ;

        Tmp(:,:,:,iChannel) = ...
            griddata( X_0, Y_0, Z_0, Shim.Operator(:,:,:,iChannel), X_1, Y_1, Z_1, interpolationMethod ) ;
    end

    Shim.Operator = Tmp ;
    % if new positions are outside the range of the original, 
    % interp3/griddata replaces array entries with NaN
    Shim.Operator( isnan( Shim.Operator ) ) = 0 ; 

    % ------------------------------------------------------------------------
    
    % ------------------------------------------------------------------------
    % Update header
    Shim.Hdr.ImageType = 'DERIVED\SECONDARY\REFORMATTED' ;

   
    Shim.Hdr.ImagePositionPatient( 1 ) = X_1(1) ; 
    Shim.Hdr.ImagePositionPatient( 2 ) = Y_1(1) ;
    Shim.Hdr.ImagePositionPatient( 3 ) = Z_1(1) ;

    %-------
    % Rows 
    Shim.Hdr.Rows = size(Shim.Operator, 1) ;
    
    dx = X_1(1,2,1) - X_1(1,1,1) ;
    dy = Y_1(1,2,1) - Y_1(1,1,1) ;
    dz = Z_1(1,2,1) - Z_1(1,1,1) ;  
    
    Shim.Hdr.PixelSpacing(1) = ( dx^2 + dy^2 + dz^2 )^0.5 ;
    
    Shim.Hdr.ImageOrientationPatient(1) = dx/Shim.Hdr.PixelSpacing(1) ;
    Shim.Hdr.ImageOrientationPatient(2) = dy/Shim.Hdr.PixelSpacing(1) ;
    Shim.Hdr.ImageOrientationPatient(3) = dz/Shim.Hdr.PixelSpacing(1) ;

    %-------
    % Columns 
    Shim.Hdr.Columns = size(Shim.Operator, 2) ;       
    
    dx = X_1(2,1,1) - X_1(1,1,1) ;
    dy = Y_1(2,1,1) - Y_1(1,1,1) ;
    dz = Z_1(2,1,1) - Z_1(1,1,1) ;  
    
    Shim.Hdr.PixelSpacing(2) = ( dx^2 + dy^2 + dz^2 )^0.5 ;
 
    Shim.Hdr.ImageOrientationPatient(4) = dx/Shim.Hdr.PixelSpacing(2) ;
    Shim.Hdr.ImageOrientationPatient(5) = dy/Shim.Hdr.PixelSpacing(2) ;
    Shim.Hdr.ImageOrientationPatient(6) = dz/Shim.Hdr.PixelSpacing(2) ;
   
    %-------
    % Slices
    Shim.Hdr.NumberOfSlices       = size(Shim.Operator, 3) ;
    Shim.Hdr.SpacingBetweenSlices = ( (X_1(1,1,2) - X_1(1,1,1))^2 + ...
                                     (Y_1(1,1,2) - Y_1(1,1,1))^2 + ...
                                     (Z_1(1,1,2) - Z_1(1,1,1))^2 ) ^(0.5) ;
    
    Shim.Hdr.SliceLocation = dot( [X_1(1) Y_1(1) Z_1(1)],  ... 
        cross( Shim.Hdr.ImageOrientationPatient(4:6), Shim.Hdr.ImageOrientationPatient(1:3) ) ) ;

end
% =========================================================================
function [] = assessshimvolume( Shim )
    fprintf( ['\n Comparing Defined Vs. Desired (Image) shim volumes\n'...
             '(In mm, relative to isocentre)\n...'] )

    [X,Y,Z] = Shim.getvoxelpositions( ) ;
    X = X( logical(Shim.mask) ) ;
    Y = Y( logical(Shim.mask) ) ;
    Z = Z( logical(Shim.mask) ) ;
    
    fprintf( ['\n Shim volume defined over ' ...
                 '\n in Z (read [rows]):    ' num2str( [ min( X(:) ) max( X(:) ) ] ) ...
                 '\n in X (p.e. [columns]): ' num2str( [ min( Y(:) ) max( Y(:) ) ] ) ...
                 '\n in Y (slice):          ' num2str( [ min( Z(:) ) max( Z(:) ) ] ) ...
                 '\n']) ;

    [X,Y,Z] = Shim.Field.getvoxelpositions( ) ;
    X = X( logical(Shim.Field.Hdr.MaskingImage) ) ;
    Y = Y( logical(Shim.Field.Hdr.MaskingImage) ) ;
    Z = Z( logical(Shim.Field.Hdr.MaskingImage) ) ;

    fprintf( ['\n Image volume defined over ' ...
                 '\n in Z (read [rows]):    ' num2str( [ min( X(:) ) max( X(:) ) ] ) ...
                 '\n in X (p.e. [columns]): ' num2str( [ min( Y(:) ) max( Y(:) ) ] ) ...
                 '\n in Y (slice):          ' num2str( [ min( Z(:) ) max( Z(:) ) ] ) ...
                 '\n']) ;

end
% =========================================================================
function Shim = forwardmodelfield( Shim )
%FORWARDMODELFIELD
    
    nVoxels = Shim.getnumberofvoxels() ;

    A = zeros( nVoxels, Shim.Parameters.nActiveChannels ) ; 
    
    for channel = 1 : Shim.Parameters.nActiveChannels
        A(:, channel) = reshape( Shim.Operator(:,:,:, channel), [nVoxels 1] ) ;
        % A(:, channel) = reshape( Shim.Field.Hdr.MaskingImage .* Shim.Operator(:,:,:, channel), [nVoxels 1] ) ;
    end

    Shim.Model.field = reshape( A*Shim.Model.currents, size( Shim.Field.img ) ) ;

end
% =========================================================================
function [Img] = interptocalibrationgrid( Shim, Img )
%RCS = (DICOM) Reference Coordinate System
    [RCS_X, RCS_Y, RCS_Z] = Shim.getvoxelpositions( ) ;

    Img = Img.resliceimg( RCS_X, RCS_Y, RCS_Z ) ;
end
% =========================================================================
function M = gettruncationoperator( Shim )
% GETTRUNCATIONOPERATOR
%
% M = GETTRUNCATIONOPERATOR( Shim ) ;
%
% Truncation (masking) operator (e.g. M*b, 'picks out' the voi voxels from 
% vector b)

nVoxelsImg = numel( Shim.Field.Hdr.MaskingImage ) ;
nVoxelsVoi = nnz( Shim.Field.Hdr.MaskingImage ) ;

indicesVoi = find( Shim.Field.Hdr.MaskingImage(:) ) ;

M = sparse( [1:nVoxelsVoi], indicesVoi, ones([nVoxelsVoi 1]), nVoxelsVoi, nVoxelsImg ) ;

% % M = sparse( find(sum(I, 2) ~= 0 ), unique(indexEp), ones([nEpRecoverable 1]), nRows, nVoxelsImg ) ;
%     
% M = NaN  
end
% =========================================================================
function A = getshimoperator( Shim )
% GETSHIMOPERATOR
%
% A = GETSHIMOPERATOR( Shim ) ;
%
%   where A * vectorOfShimCurrents = shimField

nVoxelsImg = Shim.getnumberofvoxels() ;

A = zeros( nVoxelsImg, Shim.Parameters.nActiveChannels ) ; 

for channel = 1 : Shim.Parameters.nActiveChannels
    A(:, channel) = reshape( Shim.Operator(:,:,:, channel), [nVoxelsImg 1] ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)

function [Field, Extras] = mapfield( Params )
% MAPFIELD
%
% Field = MAPFIELD( Params ) 
%
% Params must contain the following fields
%   .Path.Mag.echo1
%   .Path.Mag.echo2
%   .Path.Phase.echo1
%   .Path.Phase.echo2
%   .threshold     
%
%       .threshold (as a fraction (<1) of max measured magnitude intensity)
%           determines the phase region to be unwrapped (i.e. areas of 
%           low signal are ignored)
%       
% =========================================================================

% =========================================================================

PhaseEcho1 = MaRdI(  Params.Path.Phase.echo1  ) ;
MagEcho1   = MaRdI(  Params.Path.Mag.echo1  ) ;

PhaseEcho2 = MaRdI(  Params.Path.Phase.echo2  ) ;
MagEcho2   = MaRdI(  Params.Path.Mag.echo2  ) ;

%-------
% Mask 
PhaseEcho1.Hdr.MaskingImage = ( MagEcho1.img ) > Params.threshold ;
PhaseEcho2.Hdr.MaskingImage = ( MagEcho2.img ) > Params.threshold ;

% -------------------------------------------------------------------------
% Unwrap phase 
% (Using Abdul-Rahman metho)
PhaseEcho1 = PhaseEcho1.unwrapphase( ) ;
PhaseEcho2 = PhaseEcho2.unwrapphase( ) ;

%-------
% Compute field map
Field = PhaseEcho1.mapfrequencydifference( PhaseEcho2 ) ;
Field.Hdr.MaskingImage = PhaseEcho1.Hdr.MaskingImage .* PhaseEcho2.Hdr.MaskingImage ; 

Extras.PhaseEcho1 = PhaseEcho1 ;
Extras.PhaseEcho2 = PhaseEcho2 ;



end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
