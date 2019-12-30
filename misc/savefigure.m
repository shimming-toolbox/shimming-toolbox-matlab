function [Params] = savefigure( img, Params )
%SAVEFIGURE Write .png image to file using the 'export_fig' tool
% 
% Usage
%
%   [] = SAVEFIGURE( img, Parameters )
%   Parameters = SAVEFIGURE( img, Parameters )
%   
%   Returns employed Parameters struct.
%
% Inputs
%   
%   The following Parameter.fields are supported
%
%   .filename
%       default = './tmp'
%
%   .colormap
%       default = 'gray'
%
%   .scaling
%       default = [min(img) max(img)]
%
%   .magnification
%       default = 1
%
%   .isColorbar
%       default = false
%
%   .isBackgroundTransparent
%       default = false

DEFAULTS.filename                = './tmp' ;
DEFAULTS.colormap                = 'gray' ;
DEFAULTS.magnification           = 1 ;
DEFAULTS.isColorbar              = false ;
DEFAULTS.isBackgroundTransparent = false ;

%% ========================================================================
% Check inputs
% =========================================================================
if nargin < 1 || isempty(img)
    disp('Error: function requires at least 1 argument (2D image-matrix)')
    help(mfilename);
    return;  
end

if nargin < 2 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'scaling' ) || isempty(Params.scaling)
    Params.scaling = [min(img(:)) max(img(:))] ;
    if ( Params.scaling(2) == Params.scaling(1) )
        Params.scaling(2) = Inf ;
    end
end

Params = assignifempty( Params, DEFAULTS ) ;

[~,~,extension] = fileparts( Params.filename ) ;

if ~strcmp( extension, '.png' )
    Params.filename = [Params.filename '.png' ] ;
end

%% ========================================================================
% Create figure
% =========================================================================

figure('units','normalized','outerposition',[0 0 1 1])

imagesc( img, Params.scaling ) ; 
colormap(Params.colormap); 
axis image ;

if Params.isColorbar
    colorbar;
    hcb=colorbar;
    set(hcb,'YTick',[])
end

set(gca,'XTick',[]) % Remove the ticks in the x axis
set(gca,'YTick',[]) % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure

if Params.isBackgroundTransparent
    export_fig(Params.filename, '-png', '-transparent', ['-m' num2str(Params.magnification)]) 
else
    export_fig(Params.filename, '-png', ['-m' num2str(Params.magnification)]) 
end

% close gcf

% crop out border voxels ?
img = imread( Params.filename, 'png' ) ;
imwrite( img( 5:end-3, 5:end-3, : ), Params.filename, 'png' ) ;

end
