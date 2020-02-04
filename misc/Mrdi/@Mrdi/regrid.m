function [F] = regrid( Img, X_Ep, Y_Ep, Z_Ep, varargin ) 
%REGRID  Interpolate an image (Mrdi-type object)
% 
% In general, REGRID() uses MATLAB's scatteredInterpolant class. 
% The exception is when the image input (Img.img) is 2d and the target
% output (prescribed by inputs X,Y,Z) is a volume. This scenario is
% incompatible with scatteredInterpolant, and nearest-neighbor substitution is
% used instead.
%
% -----   
% Basic Usage
%
% [] = REGRID( Img, X, Y, Z )
% [] = REGRID( Img, X, Y, Z, mask )
% 
% Inputs:
%
% X, Y, Z:  
%       2d or 3d arrays (size=output image grid) describing the X, Y, Z patient
%       coordinates (i.e. of the DICOM reference coordinate system) of the
%       target (output) voxels. In general, if one is interpolating from one
%       image grid (Img) to another (Mrdi-type object Img2), these arrays are
%       obtained by the call: [X,Y,Z] = Img2.Grid.gridpositions()
% 
% mask: [Optional, default = true(size output image grid)]
%       A logical array (size=output image grid) specifying the subset of the
%       output voxels that are of interest. (i.e. voxels in the output image
%       with a corresponding mask entry == FALSE will simply be assigned zero).
%       Note: Specifying the region of interest for extrapolation with this
%       variable can greatly accelerate the interpolation!
%
% -----   
%
% Advanced Usage TODO
%
%   [F] = REGRID( Img, X, Y, Z, mask, F ) 
% 
%   case:
%       interpolationMethod [default='linear']
%       is a string supported by the scatteredInterpolant constructor.
%   F is the object of type 'scatteredInterpolant' used for interpolation.


%% -----
%NOTE: Terminology: 
% 'Ip' = interpolant/initial point
% 'Ep' = extrapolant/end point

F = [] ;

DEFAULT_INTERPOLATIONMETHOD  = 'linear' ;
DEFAULT_ISFORMINGINTERPOLANT = true ;

[X_Ip, Y_Ip, Z_Ip] = Img.Grid.gridpositions( ) ;

%% -----
% Parse and check inputs
if Img.Grid.comparegrids( X_Ip, Y_Ip, Z_Ip, X_Ep, Y_Ep, Z_Ep ) 
    warning('Voxel positions are already identical. Not interpolating.');
    return ;
end
        
isFormingInterpolant = true ; 

if nargin < 5

    interpolationMethod  = DEFAULT_INTERPOLATIONMETHOD ;
    isFormingInterpolant = DEFAULT_ISFORMINGINTERPOLANT ;

elseif nargin >= 5 
    if islogical( varargin{1} ) ;
        maskEp = varargin{1} ;
    end

    if nargin == 6
        if ischar( varargin{2} )
            interpolationMethod = varargin{2} ;
        elseif isa( varargin{2}, 'scatteredInterpolant' ) ;
            F = varargin{2} ;
        else
            error( 'Unknown input. See HELP MrdiProc.resliceimg' ) ;
        end
    end
end

isUsingScatteredInterpolant = [] ;

if (ndims(Img.img) > 1) && (ndims(Img.img) <= 5)

    gridSizeIp = Img.gridSize ;
    
    if gridSizeIp(3) > 1
        isUsingScatteredInterpolant = true ;
    else
        isUsingScatteredInterpolant = false ;
    end
else
    error('Dimensions of input Img.img must be >= 2, and <= 5') ;
end

if ndims( X_Ep ) == 2 % interpolating down to 2d single-slice
    gridSizeEp = [ size(X_Ep) 1 ] ;
elseif ndims( X_Ep ) == 3
    gridSizeEp = size(X_Ep) ;
else
    error('Expected 2d or 3d target interpolation grid')
end

%% -----
% Define variables
gridSizeIp = size( X_Ip ) ;

if myisfieldfilled( Img.Hdr, 'MaskingImage' ) 
    maskIp = logical( sum( sum( Img.Hdr.MaskingImage, 5 ), 4 ) ) ;
else
    maskIp = true( gridSizeIp ) ;
end

gridSizeEp = size( X_Ep ) ;
if ndims(gridSizeEp) == 2
    gridSizeEp = [gridSizeEp 1] ;
end

if ~exist('maskEp')
    if ~isUsingScatteredInterpolant 
        warning('No logical mask provided: For faster results, restrict the target/output voxels to those of interest by providing this mask!') ;
    end
    maskEp = true( gridSizeEp ) ;
end

iEp = find( maskEp(:) ) ; % indices of target voxels

iNearest = zeros( size(iEp) ) ;
nEp      = length( iEp ) ;

nImgVolumesDim4 = size(Img.img, 4 ) ; % nEchoes
nImgVolumesDim5 = size(Img.img, 5 ) ; % nMeasurements
nImgVolumes     = nImgVolumesDim4 * nImgVolumesDim5 ;

imgOut = zeros( [gridSizeEp nImgVolumesDim4 nImgVolumesDim5] ) ;

%% -------
if isFormingInterpolant 
    tic
    disp( 'Forming interpolant...' )
    disp( '(Computation time depends on input image size. This may take a few minutes.)' ) ;

    if isUsingScatteredInterpolant

        % The following avoids the error from scatteredInterpolant when one
        % attempts to form a 3d interpolant from a 2d input: 
        isValidDim0 = [ numel(unique(X_Ip(:))) numel(unique(Y_Ip(:))) numel(unique(Z_Ip(:))) ] > 1 ;
        r0          = [X_Ip(:) Y_Ip(:) Z_Ip(:)] ;
        
        isValidDim1 = [ numel(unique(X_Ep(:))) numel(unique(Y_Ep(:))) numel(unique(Z_Ep(:))) ] > 1 ;
        r1          = [X_Ep(:) Y_Ep(:) Z_Ep(:)] ;

        if nnz( isValidDim0 ) == 2
            
            assert( all( isValidDim1 == isValidDim0 ), ... 
                'Query points should sit within the same plane as the interpolant points' ) ;
                
            % coordinate of interpolant plane along normal dim:
            qn0 = unique( r0(:, ~isValidDim0) ) ;
            % coordinate of query plane along same dim:
            qn1 = unique( r1(:, ~isValidDim1) ) ;

            % This could instead be a warning? (e.g. if the 2d planes are indeed very close, should interp still be performed?)
            assert( qn0 == qn1, 'Query points should sit within the same plane as the interpolant points' ) ;

            % exclude the coordinate along the normal dim from the interpolation
            r1 = r1(:, isValidDim1) ; 

        end
        
        F                     = scatteredInterpolant() ;
        F.Method              = interpolationMethod ;
        F.ExtrapolationMethod = 'none' ;
        F.Points              = r0(:, isValidDim0) ;
    
    else % Map nearest neighbours
        
        % truncate IP voxels
        X_Ip = X_Ip(maskIp) ;
        Y_Ip = Y_Ip(maskIp) ;
        Z_Ip = Z_Ip(maskIp) ;

        for iR = 1 : nEp
            [~,iNearest(iR)] = min( sqrt( ( X_Ip - X_Ep( iEp(iR) ) ).^2 + ...
                                   ( Y_Ip - Y_Ep( iEp(iR) ) ).^2 + ...
                                   ( Z_Ip - Z_Ep( iEp(iR) ) ).^2 ) ) ;
        end
    end
    toc
end

%% -----
disp('Reslicing...')
tic
for iImgDim4 = 1 : nImgVolumesDim4
    for iImgDim5 = 1 : nImgVolumesDim5
        disp( ['Reslicing image volume...' num2str(iImgDim4*iImgDim5) ' of ' num2str(nImgVolumes) ]) ;
       
        imgIp = Img.img(:,:,:, iImgDim4, iImgDim5 ) ;
        imgEp = zeros( gridSizeEp ) ;
      
        if isUsingScatteredInterpolant  
            F.Values = imgIp(:) ;
            imgEp    = reshape( F( r1 ), gridSizeEp ) ;
        
        else % Nearest-neighbor substitution
            imgIp = imgIp(maskIp) ;
            for iR = 1 : nEp 
                imgEp( iEp(iR) ) = imgIp( iNearest(iR) ) ;
            end
        end

        imgOut(:,:,:, iImgDim4, iImgDim5 ) = imgEp ;
    end
end
toc

imgOut( isnan( imgOut ) ) = 0 ; 

Img.img = imgOut ; 

Img.Hdr.MaskingImage = Img.img ~= 0 ;
%% -----------------------------------------------------------------------

% ------------------------------------------------------------------------
% Update image + header

% Img.Hdr.ImageType = 'DERIVED\SECONDARY\REFORMATTED' ;

Img.Info.updateimage( img, X_Ep, Y_Ep, Z_Ep ) ;

% Img.Hdr.ImagePositionPatient( 1 ) = X_Ep(1) ; 
% Img.Hdr.ImagePositionPatient( 2 ) = Y_Ep(1) ;
% Img.Hdr.ImagePositionPatient( 3 ) = Z_Ep(1) ;
%
% %-------
% % Rows 
% Img.Hdr.Rows = size(Img.img, 1) ;
%
% dx = X_Ep(2,1,1) - X_Ep(1,1,1) ;
% dy = Y_Ep(2,1,1) - Y_Ep(1,1,1) ;
% dz = Z_Ep(2,1,1) - Z_Ep(1,1,1) ;  
%
% % vertical (row) spacing
% Img.Hdr.PixelSpacing(1) = ( dx^2 + dy^2 + dz^2 )^0.5 ; 
%
% % column direction cosine (expressing angle btw column direction and X,Y,Z axes)
% Img.Hdr.ImageOrientationPatient(4) = dx/Img.Hdr.PixelSpacing(1) ;
% Img.Hdr.ImageOrientationPatient(5) = dy/Img.Hdr.PixelSpacing(1) ;
% Img.Hdr.ImageOrientationPatient(6) = dz/Img.Hdr.PixelSpacing(1) ;
%
% %-------
% % Columns 
% Img.Hdr.Columns = size(Img.img, 2) ;       
%
% dx = X_Ep(1,2,1) - X_Ep(1,1,1) ;
% dy = Y_Ep(1,2,1) - Y_Ep(1,1,1) ;
% dz = Z_Ep(1,2,1) - Z_Ep(1,1,1) ;  
%
% % horizontal (column) spacing
% Img.Hdr.PixelSpacing(2) = ( dx^2 + dy^2 + dz^2 )^0.5 ;
%
% % row direction cosine (expressing angle btw column direction and X,Y,Z axes)
% Img.Hdr.ImageOrientationPatient(1) = dx/Img.Hdr.PixelSpacing(2) ;
% Img.Hdr.ImageOrientationPatient(2) = dy/Img.Hdr.PixelSpacing(2) ;
% Img.Hdr.ImageOrientationPatient(3) = dz/Img.Hdr.PixelSpacing(2) ;
%
% %-------
% % Slices
%
% if size( Img.img, 3 ) > 1
%     Img.Hdr.SpacingBetweenSlices = ( (X_Ep(1,1,2) - X_Ep(1,1,1))^2 + ...
%                                      (Y_Ep(1,1,2) - Y_Ep(1,1,1))^2 + ...
%                                      (Z_Ep(1,1,2) - Z_Ep(1,1,1))^2 ) ^(0.5) ;
% else
%     Img.Hdr.SpacingBetweenSlices = 0 ;
% end
%
% Img.setslicenormalvector() ; % redefines sHat
%
% R = Img.rotationMatrix( ) ;  
% Img.Hdr.SliceLocation = dot( Img.Hdr.ImagePositionPatient, R(:,3) ) ;

end

