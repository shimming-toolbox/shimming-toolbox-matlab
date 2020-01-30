function Interpolant = getinterpolant( Img, method, extrapolationMethod )
%GETINTERPOLANT     Return scatteredInterpolant object
% 
% GETINTERPOLANT returns an instance of Matlab's scatteredInterpolant class,
% useful for interpolating between different image grids (voxel positions).
% 
% Interpolant = GETINTERPOLANT( Img ) 
% Interpolant = GETINTERPOLANT( Img, method ) 
% Interpolant = GETINTERPOLANT( Img, method, extrapolationMethod ) 
%
% Default Interpolant property assignments:
%
%   .Method = 'linear' [i.e. the *interpolation* method]
%
%   .ExtrapolationMethod = 'none' 
%
%   .Values = [vectorized voxel values of 1st echo/measurement, i.e.: Img.img(:,:,:,1,1)]
%
% Note that all 3 properties can be reassigned at any point upon return.
%
% For info on the 2 optional arguments, 
% See also scatteredInterpolant

DEFAULT_METHOD              = 'linear' ;
DEFAULT_EXTRAPOLATIONMETHOD = 'none' ;

if ( nargin < 3 ) || isempty( extrapolationMethod )
    
    extrapolationMethod = DEFAULT_EXTRAPOLATIONMETHOD ;
    
    if ( nargin < 2 ) || isempty( method )
        method = DEFAULT_METHOD ;
    end
end

disp( 'Forming interpolant...(Computation time is proportional to image size. This may take a few minutes.)' )

Interpolant = scatteredInterpolant() ;

Interpolant.Method              = method ;
Interpolant.ExtrapolationMethod = extrapolationMethod ;

[X,Y,Z] = Img.getvoxelpositions() ;

% The following avoids the error from scatteredInterpolant when one
% attempts to form a 3d interpolant from a 2d input: 
isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
r0          = [X(:) Y(:) Z(:)] ;

Interpolant.Points = r0(:, isValidDim0) ;

v = Img.img(:,:,:,1,1) ;
Interpolant.Values = v(:) ;

end

