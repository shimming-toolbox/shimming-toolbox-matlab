classdef MrdiInterp < matlab.mixin.SetGet 
%MrdiInterp  MR DICOM Image Interpolation 
%
% DRAFT - class not yet implemented
   
% properties( Access = protected, Hidden = true )
properties( SetAccess = immutable, Hidden = true )
    
    img {mustBeNumeric} ;
    
    Grid {MrdiGrid} ;

end

properties( AbortSet = true )
    
    interpolationMask {mustBeNumericOrLogical} = true ;

end

properties( AbortSet = true, Transient = true )

    extrapolationMask {mustBeNumericOrLogical} = true ;

end

properties( AbortSet = true )

    interpolationMethod char {mustBeMember( interpolationMethod, {'linear','nearest','natural'})} = 'linear' ;
    
    extrapolationMethod char {mustBeMember( extrapolationMethod, {'linear','nearest','none'} )} = 'none' ;

end

properties( Access = protected, Hidden = true )
    Engine = scatteredInterpolant() ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Interp = MrdiInterp( Img, Params )
    
    if nargin == 0
        return ;
    elseif isa( Img, 'Mrdi' )
        Interp.img    = Img.img ;
        Interp.Grid   = Img.Grid ;
        Interp.maskIp = true( Interp.Grid.size ) ;
    else
        error('Invalid inputs') ;
    end

    if nargin == 2
        names = fieldnames( Params ) ;
        
        for iName = 1 : length( names ) 
            Interp.set( names(iName), Params(iName ) ) ;
        end
    end

    % Interp.Engine = scatteredInterpolant() ;

end
% =========================================================================    
function interpolationMethod = get.interpolationMethod( Interp )
    
    interpolationMethod = Interp.Engine.Method ;

end
% =========================================================================    
function extrapolationMethod = get.extrapolationMethod( Interp )
    
    extrapolationMethod = Interp.Engine.ExtrapolationMethod ;

end
% =========================================================================    
function [] = set.interpolationMethod( Interp, interpolationMethod )
   
    Interp.interpolationMethod = interpolationMethod ;
    Engine.Method = interpolationMethod

end
% =========================================================================    
function img1 = interpolate( Interp, img0, varargin )
    
    [X_Ip, Y_Ip, Z_Ip] = Interp.Grid.gridpositions( ) ;

    %% -----
    % Parse and check inputs
    if Img.Grid.comparegrids( X_Ip, Y_Ip, Z_Ip, X_Ep, Y_Ep, Z_Ep ) 
        warning('Voxel positions are already identical. Not interpolating.');
        return ;
    end

    img0 = img0(:) ;
    
    if numel(img0)
    Interp.Engine.Points = img0(:) ;

end
% function Interpolant = getinterpolant( Img, method, extrapolationMethod )
% %GETINTERPOLANT     Return scatteredInterpolant object
% % 
% % GETINTERPOLANT returns an instance of Matlab's scatteredInterpolant class,
% % useful for interpolating between different image grids (voxel positions).
% % 
% % Interpolant = GETINTERPOLANT( Img ) 
% % Interpolant = GETINTERPOLANT( Img, method ) 
% % Interpolant = GETINTERPOLANT( Img, method, extrapolationMethod ) 
% %
% % Default Interpolant property assignments:
% %
% %   .Method = 'linear' [i.e. the *interpolation* method]
% %
% %   .ExtrapolationMethod = 'none' 
% %
% %   .Values = [vectorized voxel values of 1st echo/measurement, i.e.: Img.img(:,:,:,1,1)]
% %
% % Note that all 3 properties can be reassigned at any point upon return.
% %
% % For info on the 2 optional arguments, 
% % See also scatteredInterpolant
%
% DEFAULT_METHOD              = 'linear' ;
% DEFAULT_EXTRAPOLATIONMETHOD = 'none' ;
%
% if ( nargin < 3 ) || isempty( extrapolationMethod )
%     
%     extrapolationMethod = DEFAULT_EXTRAPOLATIONMETHOD ;
%     
%     if ( nargin < 2 ) || isempty( method )
%         method = DEFAULT_METHOD ;
%     end
% end
%
% disp( 'Forming interpolant...(Computation time is proportional to image size. This may take a few minutes.)' )
%
% Interpolant = scatteredInterpolant() ;
%
% Interpolant.Method              = method ;
% Interpolant.ExtrapolationMethod = extrapolationMethod ;
%
% [X,Y,Z] = Img.Grid.gridpositions() ;
%
% % The following avoids the error from scatteredInterpolant when one
% % attempts to form a 3d interpolant from a 2d input: 
% isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
% r0          = [X(:) Y(:) Z(:)] ;
%
% Interpolant.Points = r0(:, isValidDim0) ;
%
% v = Img.img(:,:,:,1,1) ;
% Interpolant.Values = v(:) ;
%
% end


