function [] = fslview( img, Parameters, varargin )
%FSLVIEW view 3d array in fslview
%
%   Syntax
%
%   fslview( img )
%   fslview( img1, Parameters )
%   fslview( img1, Parameters, img2, img3, ... )
%
%   Description
%
%   creates temporary nifti conversion of img and passes it to fslview
%
%   
%       Parameters (object) to fslview input 
%       (see manual 'man fslview')
%       
%   The following Parameters fields are supported
%   
%   .voxelSize (3-component vector)
%       default:
%           [1 1 1] ;
%
%   .l (color map: Full list avaible in FSLView GUI)  
%       (e.g. one of: Greyscale, Red-Yellow, Blue-Lightblue,
%       Red, Green, Blue, Yellow, Pink, Hot, Cool, Copper, ...)
%       default: 
%           Parameters.l = 'Greyscale'
%
%   .b (brightness)
%       if scalar: 
%           brightness = [-Abs|.b|,+Abs|.b|] 
%       if 2 component vector: 
%           brightness = [.b(1), .b(2)] 
%       default: 
%           brightness = [min(img), max(img)] 
%
%
% =========================================================================
% Ryan Topfer 2012
% topfer@ualberta.ca

% TODO
% variable brightness and colorschme for each overlay

% =========================================================================
% Check inputs // Assign defaults
% =========================================================================
DEFAULT_COLORMAP  = 'Greyscale' ;
DEFAULT_VOXELSIZE = [1 1 1] ;

if nargin == 1 || isempty(Parameters)
    Parameters.dummy = [] ; %using defaults
end

if  ~myisfield( Parameters, 'voxelSize' ) || isempty(Parameters.voxelSize)
    Parameters.voxelSize = DEFAULT_VOXELSIZE ;
end

%------ 
% Brightness levels

if ~myisfield( Parameters, 'b' ) || isempty(Parameters.b)
    Parameters.b = [min(img(:)) max(img(:))] ;
    
elseif max( size( Parameters.b ) ) == 1 ;
    Parameters.b(2) = abs(Parameters.b) ;
    Parameters.b(1) = -abs(Parameters.b(1)) ;
end

viewParameters = [' -b "' num2str(Parameters.b(1)) ',' num2str(Parameters.b(2)) '"'];

%------ 
% Color scheme

if myisfield( Parameters, 'l' )
    if isempty(Parameters.l)
        Parameters.l = DEFAULT_COLORMAP ; 
    end
    viewParameters = [ viewParameters ' -l ' Parameters.l ] ;
end

filenames = {['view3dtmp1.nii']} ; 
call2Fslview = strcat( filenames(1), viewParameters ) ;

% =========================================================================

% =========================================================================
Parameters.filename  = filenames{1} ;
nii( img, Parameters ) ;

%------ 
% Additional input images?
nImgs = nargin - 2 ;

for iImg = 1 : nImgs
   
    filenames( iImg + 1 ) = {['view3dtmp' int2str( iImg + 1 ) '.nii']} ;
    Parameters.filename = filenames( iImg + 1 ) ;
    nii( varargin{iImg}, Parameters ) ;

    call2Fslview = strcat( call2Fslview, {' '}, Parameters.filename, viewParameters) ;
end
    
call2Fslview = ['fslview ', char(call2Fslview), ' & '] ;

system( call2Fslview ) ;
pause(5)

% =========================================================================
% Clean up
% =========================================================================
delete('view3dtmp*') ;

end
