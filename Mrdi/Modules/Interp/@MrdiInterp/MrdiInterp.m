classdef MrdiInterp 
%MrdiInterp  MR DICOM Image Interpolation 
%
% DRAFT - class not yet implemented
   
% properties( Access = protected, Hidden = true )

properties( SetAccess=immutable )
    
    Grid MrdiGrid ;
    
    img {mustBeNumeric} ;

end

properties( AbortSet )
    
    interpolationMask {mustBeNumericOrLogical} = true ;

end

properties( AbortSet, Transient )

    extrapolationMask {mustBeNumericOrLogical} = true ;

end

properties( AbortSet )

    interpolationMethod char {mustBeMember( interpolationMethod, {'linear','nearest','natural'})} = 'linear' ;
    
    extrapolationMethod char {mustBeMember( extrapolationMethod, {'linear','nearest','none'} )} = 'none' ;

end

properties( Access=protected, Hidden )

    Engine   = scatteredInterpolant() ;
    Engine2d = [] ;

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Interp = MrdiInterp( Img, Params )
    
    if nargin == 0
        return ;
    elseif nargin == 1
        Params = struct([]) ;
    end

    assert( isa( Img, 'Mrdi' ), 'First input must be a Mrdi-type image object' ) ;
    assert( isstruct( Params ), 'Second input must be a parameters struct' ) ;
    
    % Grid property SetAccess=immutable (must be defined here in the constructor)
    Interp.Grid = Img.Grid ;
    Interp.img  = Img.img ;

    Interp = Interp.configure( Params ) ;

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
    
    % if numel(img0)
    % Interp.Engine.Points = img0(:) ;

end
% =========================================================================

% =========================================================================
% =========================================================================
end

methods( Access=protected, Hidden )
% =========================================================================    
function Interp = configure( Interp, Params )
%CONFIGURE  

DEFAULTS.interpolationMethod = 'linear' ;

if Interp.Grid(3) == 1
    DEFAULTS.extrapolationMethod = 'nearest' ;
else
    DEFAULTS.extrapolationMethod = 'linear' ;

    Params = assignifempty( Params, DEFAULTS ) ;
% DEFAULTS3D.interpolationMask   = 
% Interp.interpolationMethod char {mustBeMember( interpolationMethod, {'linear','nearest','natural'})} = 'linear' ;
%     
%     extrapolationMethod char {mustBeMember( extrapolationMethod, {'linear','nearest','none'} )} = 'none' ;
%
% if Interp.Grid.size(3) > 1
%     
%     Interp.Engine = scatteredInterpolant() ;


% if nargin == 2
    %     Interp.assigndefaults( Params ) ;
%         assert( isstruct( Params )
%     if nargin == 2
%         Interp.maskIp = true( Interp.Grid.size ) ;
%         names = fieldnames( Params ) ;
%         
%         for iName = 1 : length( names ) 
%             Interp.set( names(iName), Params(iName ) ) ;
%         end
%     end

    % Interp.Engine = scatteredInterpolant() ;

end
% =========================================================================    

end
% =========================================================================    
% =========================================================================    
end

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

[X,Y,Z] = Img.Grid.gridpositions() ;

% The following avoids the error from scatteredInterpolant when one
% attempts to form a 3d interpolant from a 2d input: 
isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
r0          = [X(:) Y(:) Z(:)] ;

Interpolant.Points = r0(:, isValidDim0) ;

v = Img.img(:,:,:,1,1) ;
Interpolant.Values = v(:) ;

end
