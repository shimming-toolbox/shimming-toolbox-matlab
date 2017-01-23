classdef ShimEval < ShimOpt 
%SHIMEVAL - Shim Evaluation 
%
% .......
% 
% Usage
%
% Shim = ShimOpt( )
%
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
% ShimEval is a ShimOpt subclass [ShimEval < ShimOpt < MaRdI]
%     
% =========================================================================
% Updated::20160925::ryan.topfer@polymtl.ca
% =========================================================================

properties
    ShimmedField;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimEval( Params )
%SHIMEVAL - Shim Evaluation 

if nargin < 1
    Params = [];
end

Shim = Shim@ShimOpt( Params ) ;

Shim.ShimmedField = [] ;

end
% =========================================================================
function Shim = histogramfield( Shim )
%histogram field
end
% =========================================================================
function [] = saveimages( Shim, Params )
%SAVEIMAGES
%
% Params can have the following fields
%   
%   .imgSlice 
%       index of the slice to image 
%       default: 1
%
%   .mediaSaveDir 
%       directory name in which images are saved 
%       default: ['./img_ ' datestr(now,30)] '/']
%
%   .scaling
%       brightness window limits in Hz 
%       default: [min(img(:)) max(img(:))]
%   
%   .colorMap
%       default: 'gray'
%
%   .isColorBar
%       default: false 

DEFAULT_MEDIASAVEDIR = ['./img_' datestr(now,30) '/'] ;
DEFAULT_IMGSLICE = 1 ;

if  ~myisfield( Params, 'mediaSaveDir' ) || isempty(Params.mediaSaveDir)
    Params.mediaSaveDir = DEFAULT_MEDIASAVEDIR ;
end

if  ~myisfield( Params, 'imgSlice' ) || isempty(Params.imgSlice)
    Params.imgSlice = DEFAULT_IMGSLICE ;
end

mkdir( Params.mediaSaveDir ) ;

% -------
img = double( Shim.Field.img ) ;
Params.filename = [Params.mediaSaveDir 'original_field'] ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* Shim.Field.img ) ;
Params.filename = [Params.mediaSaveDir 'original_field_masked'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice), Params ) ;

% -------
img = double( Shim.ShimmedField.img ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_field'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* Shim.ShimmedField.img ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_field_masked'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.ShimmedField.img - Shim.Field.img ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_diff'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* (Shim.ShimmedField.img - Shim.Field.img ) ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_diff_masked'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Model.field ) ;
Params.filename = [Params.mediaSaveDir 'predicted_shim_field'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* Shim.Model.field ) ;
Params.filename = [Params.mediaSaveDir 'predicted_shim_field_masked'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% ------
img = double(Shim.Field.img + Shim.Model.field) ;
Params.filename = [Params.mediaSaveDir 'predicted_shimmed_field'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ) , Params ) ;

% -------
img = double(Shim.Field.Hdr.MaskingImage .* ( Shim.Field.img + Shim.Model.field ) )  ;
Params.filename = [Params.mediaSaveDir 'predicted_shimmed_field_masked'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;

% ------
Params.scaling  = [0 max(abs(Params.scaling))] ;
Params.colormap = 'gray' ;
img = double(abs((Shim.Field.img + Shim.Model.field) - Shim.ShimmedField.img)) ;
Params.filename = [Params.mediaSaveDir 'discrepancy'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ) , Params ) ;

% ------
img = double(abs(Shim.Field.Hdr.MaskingImage .* ((Shim.Field.img + Shim.Model.field) - Shim.ShimmedField.img))) ;
Params.filename = [Params.mediaSaveDir 'discrepancy_masked'];
% nii( img, Params ) ;
writeimg( img( :,:, Params.imgSlice ) , Params ) ;



% -------
Params.scaling  = [0 1] ;
Params.colormap = 'gray' ;
img =  double(Shim.Field.Hdr.MaskingImage);
Params.filename = [Params.mediaSaveDir 'shim_mask'];
% nii( Shim.Field.Hdr.MaskingImage, Params ) ;
writeimg( img( :,:, Params.imgSlice ), Params ) ;




end
% =========================================================================
function Results = assessfielddistribution( Shim, mask )
%ASSESSSHIM
%
% Results = ASSESSFIELDDISTRIBUTION( Shim, mask )
% 
% mask is a binary array the same size as Shim.img, Shim.Field.img, etc.
% indicating the region of interest over which shimming was performed. 
%
% Results contains fields
%
%   .shimVolume
%       shimmed volume of interest [units: cm^3]
%
    
if nargin < 2 || isempty(mask)
    mask = Shim.Field.Hdr.MaskingImage ;
end

mask = logical( mask ) ;

Results.shimVolume = nnz( mask ) .* prod( 0.1*Shim.getvoxelsize() )  ; % [units: cm^3]

originalField     = Shim.Field.img( mask ) ;
shimmedModelField = Shim.Field.img( mask ) + Shim.Model.field( mask ) ;
shimmedField      = Shim.ShimmedField.img( mask ) ; 

Results.Std.Original.data  = std( originalField ) ;
Results.Std.Shimmed.model  = std( shimmedModelField ) ;
Results.Std.Shimmed.data   = std( shimmedField ) ;

Results.PercentImprovement.model  = ...
    100*( 1 - Results.Std.Shimmed.model/Results.Std.Original.data ) ;

Results.PercentImprovement.data  = ...
    100*( 1 - Results.Std.Shimmed.data/Results.Std.Original.data ) ;

Results.percentDiscrepancy = ...
    median( 100*( abs(shimmedModelField - shimmedField)./abs(shimmedModelField) ) ) ;

end
% =========================================================================
function Shim = setshimmedfield( Shim, ShimmedField )
%SETSHIMMEDFIELD
%
% Shim = SETSHIMMEDFIELD( Shim, ShimmedField )
%
% Sets Shim.ShimmedField
%
% ShimmedField is a MaRdI type object with .img in Hz

Shim.ShimmedField = ShimmedField ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
function [] = writeimg( img, Params )
%WRITEIMG
%
%   Description
%   
%   Write .png image to file using the 'export_fig' tool
%
%    .......................
%
%   Syntax
%
%   WRITEIMG( img, Parameters )
%
%    .......................
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
%
% =========================================================================
% Created : Sat 23 Jan 2016 20:18:54 EST 
% topfer@ualberta.ca
% =========================================================================

DEFAULT_FILENAME      = './tmp.bin' ;
DEFAULT_COLORMAP      = 'gray' ;
DEFAULT_MAGNIFICATION = 1 ;
DEFAULT_ISCOLORBAR    = false ;
DEFAULT_ISBACKGROUNDTRANSPARENT = false ;
% =========================================================================
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

if  ~myisfield( Params, 'filename' ) || isempty(Params.filename)
    Params.filename = DEFAULT_FILENAME ;
end

if  ~myisfield( Params, 'colormap' ) || isempty(Params.colormap)
    Params.colormap = DEFAULT_COLORMAP ;
end

if  ~myisfield( Params, 'magnification' ) || isempty(Params.magnification)
    Params.magnification = DEFAULT_MAGNIFICATION ;
end

if  ~myisfield( Params, 'isColorbar' ) || isempty(Params.isColorbar)
    Params.isColorbar = DEFAULT_ISCOLORBAR ;
end

if  ~myisfield( Params, 'isBackgroundTransparent' ) || isempty(Params.isBackgroundTransparent)
    Params.isBackgroundTransparent = DEFAULT_ISBACKGROUNDTRANSPARENT ;
end

if  ~myisfield( Params, 'scaling' ) || isempty(Params.scaling)
    Params.scaling = [min(img(:)) max(img(:))] ;
end

[~,~,extension] = fileparts( Params.filename ) ;

if ~strcmp( extension, '.png' )
    Params.filename = [Params.filename '.png' ] ;
end

% =========================================================================
% Create figure
% =========================================================================


figure
Params

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

end

end
