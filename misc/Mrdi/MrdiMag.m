classdef MrdiMag < MrdiInfo 
%MrdiMag  MR DICOM Image Magnitude 
%
% Member methods for magnitude images 
% 
% e.g.
%   getreliabilitymask()
%   segmentspinalcanal()
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
function Mag = MrdiMag( Img )
   
    if ( nargin < 1 ) || isempty( Img )
        return ;
    elseif ( nargin ~= 1 ) || ~isa( Img, 'MaRdI' ) || ~Img.ismagnitude()
        error('Constructor requires 1 input argument: a Mardi object initialized with phase data. See doc Mardi.') ;
    end

    Mag.copyproperties( Img ) ;   

end    
% =========================================================================
function mask = getreliabilitymask( Mag, threshold )
%GETRELIABILITYMASK
% 
%  mask = getreliabilitymask( Mag )
%  mask = getreliabilitymask( Mag, threshold )
% 
%  For each echo and each measurement, GETRELIABILITYMASK normalizes
%  Mag.img(:,:,:,iEcho,iMeasurement) and returns a logical mask wherein
%  elements are assigned TRUE whenever the corresponding normalized magnitude
%  voxel is > threshold
%
%  By default, threshold = 0.01

assert( Mag.ismagnitude(), 'Function expected magnitude image input' ) ;

DEFAULTS.threshold = 0.01 ;

if ( nargin == 1 ) || isempty( threshold )
    threshold = DEFAULTS.threshold ;
end

mask = false( size( Mag.img ) ) ;

for iVolume = 1 : size( mask, 5 ) 
    for iEcho = 1 : size( mask, 4 ) 
        mag = Mag.img(:,:,:,iEcho, iVolume) ;
        mask(:,:,:,iEcho,iVolume) = ( mag./max(mag(:)) ) > threshold ;
    end
end

end
% =========================================================================
function [mask, weights] = segmentspinalcanal( Img, Params )
%SEGMENTSPINALCANAL
% 
% segment T2* multiecho data using the Spinal Cord Toolbox (must be installed + in path)
%
% [ mask, weights ] = SEGMENTSPINALCANAL( Img, Params )
%
% Params
%
%   .dataLoadDir 
%       DICOM folder
%   
%   .dataSaveDir 
%
%   .isUsingPropsegCsf
%       [default = false]
%
% NOTE
%   The protocol is basically that of Topfer R, et al. Magn Reson Med, 2018. 
%   It hasn't been tested extensively for different acquisition prtocols/systems

mask = false ;

if nargin < 2 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'dataLoadDir' ) || isempty(Params.dataLoadDir)
    [Params.dataLoadDir,~,~] = fileparts( Img.Hdr.Filename ) ;
end

[mask, weights] = Mardi.segmentspinalcanal_s( Params ) ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [mask, weights] = segmentspinalcanal_s( Params )
%SEGMENTSPINALCANAL_S
% 
% segment T2* multiecho data using the Spinal Cord Toolbox
%
% [ mask, weights ] = SEGMENTSPINALCANAL_S( Params )
%
% Params
%
%   .dataLoadDir 
%       DICOM folder
%   
%   .dataSaveDir
%       [default = './gre_seg/'] 
%
%   .isUsingPropsegCsf
%       [default = false]
%
%   .centerlineMethod
%       Method used to obtain the centerline: 
%           'midfov': mask is centered in the middle of the axial FOV
%           'spinalcord': mask follows the spinal cord centerline
%       [default = 'midfov']
%
% NOTE
%   The protocol is basically that of Topfer R, et al. Magn Reson Med, 2018. 
%   It hasn't been tested extensively for different acquisition protocols or systems 
%
% TODO
% SEGMENTSPINALCANAL_S
%   is the static form of Mardi.segmentspinalcanal( Img, Params )
%
%   I (RT) was hoping Matlab would allow the 2 identically named methods (as in C)
%   given that the static form takes only 1 arg, and the other form requires 2...
%
%   --> either rectify this if pos. or change the method names or another alt.

%   .isForcingOverwrite
%       if .nii or .gz files exist already in dataSaveDir they will be deleted
%       [default = false]
%

mask = false ;

DEFAULT_DATALOADDIR        = [] ; % path to dicom folder
DEFAULT_DATASAVEDIR        = './gre_seg/'
DEFAULT_ISFORCINGOVERWRITE = false ;
DEFAULT_ISUSINGPROPSEGCSF  = true ; %use the propseg -CSF option
DEFAULT_CENTERLINEMETHOD   = 'midfov' ;
DEFAULT_CYLINDERSIZE       = 40 ;
DEFAULT_GAUSSIANSIZE       = 20 ;

if nargin < 1 || isempty(Params) || ~myisfield( Params, 'dataLoadDir' ) || isempty(Params.dataLoadDir)
    error('Function requires struct Params. with Params.dataLoadDir defined. See documentation.')
end

if  ~myisfield( Params, 'dataSaveDir' ) || isempty(Params.dataSaveDir)
    Params.dataSaveDir = DEFAULT_DATASAVEDIR ;
elseif ( Params.dataSaveDir(end) ~= '/' )
    Params.dataSaveDir(end+1) = '/';
end

if  ~myisfield( Params, 'isForcingOverwrite' ) || isempty(Params.isForcingOverwrite)
    Params.isForcingOverwrite = DEFAULT_ISFORCINGOVERWRITE ;
end

if  ~myisfield( Params, 'cylinderSize' ) || isempty(Params.cylinderSize)
    Params.cylinderSize = DEFAULT_CYLINDERSIZE ;
end

if  ~myisfield( Params, 'gaussianSize' ) || isempty(Params.gaussianSize)
    Params.gaussianSize = DEFAULT_GAUSSIANSIZE ;
end

if  ~myisfield( Params, 'centerlineMethod' ) || isempty(Params.centerlineMethod)
    Params.centerlineMethod = DEFAULT_CENTERLINEMETHOD ;
end

Params.tmpSaveDir = [ Params.dataSaveDir 'tmp_sct_' datestr(now, 30) '/'] ;
mkdir( Params.tmpSaveDir )

% if ~Params.isForcingOverwrite & exist( Params.dataSaveDir )
%     error('Params.dataSaveDir should not exist, or use input option Params.isForcingOverwrite == true')
% end

if  ~myisfield( Params, 'isUsingPropsegCsf' ) || isempty(Params.isUsingPropsegCsf)
    Params.isUsingPropsegCsf = DEFAULT_ISUSINGPROPSEGCSF ;
end

if ~exist( Params.dataSaveDir )
    mkdir( Params.dataSaveDir ) ;
end

dicm2nii( Params.dataLoadDir, Params.tmpSaveDir )

% rename
system( ['mv ' Params.tmpSaveDir '*.nii.gz ' Params.tmpSaveDir 't2s_allEchoes.nii.gz'] ) ;

% average across echoes
system( ['sct_maths -i ' Params.tmpSaveDir 't2s_allEchoes.nii.gz -mean t -o ' Params.tmpSaveDir 't2s.nii.gz'] ) ;

% switch between methods for obtaining a pixel location per slice
if Params.centerlineMethod == 'midfov'
    % create a vertical line centered in the axial FOV
    system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p center -size 1 -f box -o ' Params.tmpSaveDir 't2s_centerline.nii.gz' ] ) ;
elseif Params.centerlineMethod == 'spinalcord'
    % get cord centerline
    % TODO: make the param -c an input Params.
    system( ['sct_get_centerline -i ' Params.tmpSaveDir 't2s.nii.gz -c t2 -o ' Params.tmpSaveDir 't2s_centerline'] ) ;
end

% create a binary cylindrical mask around the centerline
system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p centerline,' ...
    Params.tmpSaveDir 't2s_centerline.nii.gz -size ' ...
    num2str(Params.cylinderSize) 'mm -f cylinder -o ' ...
    Params.tmpSaveDir 't2s_seg.nii.gz' ] ) ;

% create a soft gaussian mask around the centerline
system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p centerline,' ...
    Params.tmpSaveDir 't2s_centerline.nii.gz -size ' ...
    num2str(Params.gaussianSize) 'mm -f gaussian -o ' ...
    Params.tmpSaveDir 't2s_weights.nii.gz' ] ) ;

% unzip the image volumes we wish to keep
system( ['gunzip ' Params.tmpSaveDir 't2s_seg.nii.gz -df'] ) ;
system( ['gunzip ' Params.tmpSaveDir 't2s_weights.nii.gz -df'] ) ;

% delete the other images
system( ['rm ' Params.tmpSaveDir '*.nii.gz'] ) ;

% move segmentation + weights
system( ['mv ' Params.tmpSaveDir 't2s_seg.nii ' Params.dataSaveDir 'gre_seg.nii'] ) ;  
system( ['mv ' Params.tmpSaveDir 't2s_weights.nii ' Params.dataSaveDir 'gre_weights.nii'] ) ;  
system( ['rm -r ' Params.tmpSaveDir] ) ;

Mask = load_untouch_nii( [ Params.dataSaveDir 'gre_seg.nii' ] );
mask = Mask.img ;
mask = logical(permute( mask, [2 1 3] )) ;
mask = flipdim( mask, 1 ) ;

Weights = load_untouch_nii( [ Params.dataSaveDir 'gre_weights.nii' ] );
weights = Weights.img ;
weights = double(permute( weights, [2 1 3] )) ;
weights = flipdim( weights, 1 ) ;
% normalize
weights = weights - min(weights(:)) ;
weights = weights/max(weights(:)) ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
end
