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
% Updated::20170204::ryan.topfer@polymtl.ca
% =========================================================================

properties
    ShimmedField;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimEval( Shims )
%SHIMEVAL - Shim Evaluation 

if nargin < 1
    Params = [];
end

Shim = Shim@ShimOpt( Params ) ;

Shim.ShimmedField = [] ;

end
% =========================================================================
function Shim = histogramfield( Shim )
%HISTOGRAMFIELD



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
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* Shim.Field.img ) ;
Params.filename = [Params.mediaSaveDir 'original_field_masked'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice), Params ) ;

% -------
img = double( Shim.ShimmedField.img ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_field'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* Shim.ShimmedField.img ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_field_masked'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.ShimmedField.img - Shim.Field.img ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_diff'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* (Shim.ShimmedField.img - Shim.Field.img ) ) ;
Params.filename = [Params.mediaSaveDir 'shimmed_diff_masked'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Model.field ) ;
Params.filename = [Params.mediaSaveDir 'predicted_shim_field'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% -------
img = double( Shim.Field.Hdr.MaskingImage .* Shim.Model.field ) ;
Params.filename = [Params.mediaSaveDir 'predicted_shim_field_masked'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% ------
img = double(Shim.Field.img + Shim.Model.field) ;
Params.filename = [Params.mediaSaveDir 'predicted_shimmed_field'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ) , Params ) ;

% -------
img = double(Shim.Field.Hdr.MaskingImage .* ( Shim.Field.img + Shim.Model.field ) )  ;
Params.filename = [Params.mediaSaveDir 'predicted_shimmed_field_masked'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;

% ------
Params.scaling  = [0 max(abs(Params.scaling))] ;
Params.colormap = 'gray' ;
img = double(abs((Shim.Field.img + Shim.Model.field) - Shim.ShimmedField.img)) ;
Params.filename = [Params.mediaSaveDir 'discrepancy'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ) , Params ) ;

% ------
img = double(abs(Shim.Field.Hdr.MaskingImage .* ((Shim.Field.img + Shim.Model.field) - Shim.ShimmedField.img))) ;
Params.filename = [Params.mediaSaveDir 'discrepancy_masked'];
% nii( img, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ) , Params ) ;



% -------
Params.scaling  = [0 1] ;
Params.colormap = 'gray' ;
img =  double(Shim.Field.Hdr.MaskingImage);
Params.filename = [Params.mediaSaveDir 'shim_mask'];
% nii( Shim.Field.Hdr.MaskingImage, Params ) ;
MaRdI.writeimg( img( :,:, Params.imgSlice ), Params ) ;




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


Results.Std.Original.data  = std( originalField ) ;
Results.Std.Shimmed.model  = std( shimmedModelField ) ;

Results.PercentImprovement.model  = ...
    100*( 1 - Results.Std.Shimmed.model/Results.Std.Original.data ) ;

if myisfield( Shim, 'ShimmedField' ) && ~isempty( Shim.ShimmedField ) 
    shimmedField      = Shim.ShimmedField.img( mask ) ; 
    
    Results.Std.Shimmed.data   = std( shimmedField ) ;
    Results.PercentImprovement.data  = ...
        100*( 1 - Results.Std.Shimmed.data/Results.Std.Original.data ) ;
    Results.percentDiscrepancy = ...
        median( 100*( abs(shimmedModelField - shimmedField)./abs(shimmedModelField) ) ) ;
end


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
% =========================================================================
function [] = plotfieldhistogram( field, Params )
%PLOTFIELDHISTOGRAM
% 
% Params
%
%   .xLimits
%        horizontal limits in Hz 
%        [default : max(abs(field(:))) .* [ -1 1 ] ]
%
%   .yLimits
%        vertical limits as count #
%
%   .fontSize
%       [default : 16]
%   .yAxisLocation
%       [default : 'left'] 
%
%   .textContent
%        string appears in textbox (e.g. 'Standard Deviation = 5 Hz')
%        [default : '']
%
%   .textPosition
%        [X Y]
%        position of textbox relative to graph's origin 


DEFAULT_ISTEXTBOX     = false ;
DEFAULT_YAXISLOCATION = 'left' ;
DEFAULT_FONTSIZE     = 16 ;

% -------
% Check inputs
% -------
if nargin < 1 || isempty(field) 
    disp('Error: function require 1 argument (field vector, i.e. delta B0)')
    help(mfilename);
    return;  
end

if nargin < 2 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'xLimits' ) || isempty( Params.xLimits )
    Params.xLimits = max(abs(field(:))) .* [ -1 1 ] ;
end

assert( length(Params.xLimits) == 2 ) ;


if  ~myisfield( Params, 'xTick' ) 
    midPoint = (max(Params.xLimits) + min(Params.xLimits))/2 ;
    Params.xTick = [round(Params.xLimits(1)/2) midPoint round(Params.xLimits(2)/2)] ;
end


if  ~myisfield( Params, 'yAxisLocation' ) || isempty(Params.yAxisLocation)
    Params.yAxisLocation = DEFAULT_YAXISLOCATION ;
end

if  ~myisfield( Params, 'fontSize' ) || isempty(Params.fontSize)
    Params.fontSize = DEFAULT_FONTSIZE ;
end

if  ~myisfield( Params, 'textPosition' ) || isempty(Params.textPosition) || ...
    ~myisfield( Params, 'textContent' ) || isempty(Params.textContent)
    Params.isTextBox = DEFAULT_ISTEXTBOX ;
else
    Params.isTextBox = true;
end

scnsize = get(0,'ScreenSize');
figure('Color',[1 1 1],...
    'Position',scnsize-[-30 -40 60 120]);

if ~myisfield( Params, 'binWidth') || isempty( Params.binWidth )
    Histo = histogram( field, 'BinMethod', 'sqrt' ) ;
    Params.nBins = ceil( sqrt(numel(field(:)))) ;
else
    Histo = histogram( field, Params.nBins, 'BinWidth', Params.binWidth );
end

if ~myisfield( Params, 'yLimits' )
    Params.yLimits = [ min( Histo.Values ) round(1.05*max( Histo.Values ) ) ] ;
end

if  ~myisfield( Params, 'yTick' ) 
    midPoint = round((max(Params.yLimits) + min(Params.yLimits))/2) ;
    Params.yTick = [Params.yLimits(1) midPoint Params.yLimits(2)] ;
end


disp('Plot Params :')
disp(' ')
disp(' .xLimits')
disp([' ' num2str(Params.xLimits)])
disp(' .yLimits')
disp([' ' num2str(Params.yLimits)])
disp(' .nBins')
disp([' ' num2str(Params.nBins)])

Histo


ylabel('Counts');
xlabel('Field (Hz)')

axis([Params.xLimits Params.yLimits]) ;

ax = gca ;

ax.XTick = Params.xTick ; 
ax.YTick = Params.yTick ;

ax.FontSize = Params.fontSize ;
ax.AmbientLightColor =[1 1 1];

ax.YAxisLocation = Params.yAxisLocation ; 

if Params.isTextBox
    text(Params.textPosition(1),Params.textPosition(2),...
        Params.textContent,...
        'FontSize', Params.fontSize ) ;
end

end

end
% =========================================================================
% =========================================================================

end
