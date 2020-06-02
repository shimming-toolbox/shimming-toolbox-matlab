function [img,info,json] = read_nii( niiFile, options )
%READ_NII Load NIfTI image, header, and (if present) the accompanying .json sidecar
%     
%     [img, info]       = read_nii( niiFile )
%     [img, info, json] = read_nii( niiFile, options )
% 
% The input `niiFile` is the path to the NIfTI image as a string scalar or
% character vector. When called with 2 output arguments, the function is
% equivalent short-hand for
%     
%     info = niftiinfo( niiFile ); img = niftiread( info ) 
%
% The function checks the `niiFile` parent-folder for the presence of a
% sidecar (an identically named file, but with a .json file extension).
% When such a file exists, the 3rd output is returned as a struct via 
% `json = jsondecode( fileread( jsonFile ) );` otherwise, `json = []`.
% 
% __OPTIONS__
% 
% The function accepts an `options` struct of parameters (for now, only one) 
% as a 2nd argument for which the `.rescale` field can be assigned:
%
% | `options.rescale` |  Effect                                                |
% | ------------------| -------------------------------------------------------|
% | 'off'             | Rescaling disabled                                     |
% | 'basic'           | Rescale according to nii header info                   | 
% | 'auto' [default]  | Rescale and convert to physical units when possible    |
%
% __NOTE__  
% For now, the sole effect of `'auto'` is to convert Siemens raw phase images
% to physical units (radians), which requires converting from their original
% integer type (between [0,4095]) to a 32-bit "single" float (between [-pi,pi)). 
% The json sidecar must be available to verify the Manufacturer and ImageType
% entries. Otherwise, and for all other image inputs `'auto'` reverts to
% `'basic'`.

DEFAULTS.rescale = 'auto';

narginchk(1, 2);

if nargin == 1 || ~isfield( options, 'rescale' ) || isempty( options.rescale )
    options.rescale = DEFAULTS.rescale;
else
    assert( isstruct(options), 'Optional 2nd input must be a struct' );
    assert( ismember( options.rescale, ["off" "basic" "auto"] ), ...
        'Invalid assignment to `options.rescale`: Value must be "off", "basic", or "auto"]' );
end

info = niftiinfo( niiFile ); 
img  = niftiread( info ); 

% `extractBefore` should get the correct filename in both .nii and .nii.gz cases
jsonFile = strcat( extractBefore( niiFile, '.nii' ), '.json' ); 

if isfile( jsonFile )
    if isOctave()
      json = loadjson( jsonFile );
    else % matlab
      json = jsondecode( fileread( jsonFile ) );  
    end
else
    json = [];
end

% -------------------
%% Optional rescaling
if strcmp(options.rescale, 'off')
    return
end

if strcmp(options.rescale, 'basic')
    [img, info] = rescale( img, info );
    return
end

% NOTE: Other approaches to rescaling and/or converting from raw file values
% could be added (including cases where the json sidecar is unavailable)
if isfield( json, 'Manufacturer') && strcmp(json.Manufacturer, 'Siemens') ...
  && strcmp( check_json_image_type(json), 'phase' )

    [img, info] = convert_siemens_phase( img, info );
else 
    [img, info] = rescale( img, info );
end

end

% -----------------------------------------------------------------------------
%% Local functions
% -----------------------------------------------------------------------------
function [img, info] = convert_siemens_phase( img0, info0 )
%CONVERT_SIEMENS_PHASE Return phase img in rad 
%    
%    [img, info] = convert_siemens_phase( img0, info0 )
%
% Converts from integer-type (between [0,4095]) to a 32-bit "single" float
% (between [-pi,pi)). 

PHASE_SCALING_SIEMENS = 4096;

if (info0.AdditiveOffset == -PHASE_SCALING_SIEMENS) && (info0.MultiplicativeScaling == 2) 

    [img, info] = rescale( img0, info0 );

    img = single(img)*(pi/PHASE_SCALING_SIEMENS);

    %% Update header: both Matlab simplified and original `raw` portions
    info.Datatype     = 'single';
    info.BitsPerPixel = 32;
    % NB: 16 is the NIfTI code for float; bitpix = number of bits
    info.raw.datatype = 16;
    info.raw.bitpix   = 32;
else
    warning( 'The nii header differs from that expected of Siemens phase data.\n%s', ...
              'Output values (units) are effectively unknown' );
    [img, info] = rescale( img0, info0 );
end

end

% -----------------------------------------------------------------------------
function [imgType] = check_json_image_type( json )
%CHECK_JSON_IMAGE_TYPE Check json; return 'magnitude', 'phase' or 'unknown'
%    
%    imgType = check_json_image_type( json )
%
% Checks the `ImageType` field of struct `json` (derived from the json sidecar) 
% and returns a char vector `imgType` which may be:  
%
% - `magnitude` for magnitude images
% - `phase` for phase images
% - `unknown` otherwise 
%
% __NOTE__
% The check is hardly exhaustive and has only been tested for Siemens MRI data.
assert( nargin<2 )

if isempty(json)
    imgType = []; 
    return
end
    
assert( isstruct(json) );

imgType = 'unknown'; % default

if myisfield( json, 'ImageType' )

    isPhase = any( ismember(json.ImageType, "P") );
    isMag   = any( ismember(json.ImageType, "M") ); 

    if isPhase && isMag 
        % Both true: json file and/or DICOM issue (hopefully this doesn't occur?) 
        warning('Ambiguous ImageType entry in json file: Indicates magnitude AND phase?');
    elseif isPhase
        imgType = 'phase';
    elseif isMag
        imgType = 'magnitude';
    end

    if strcmp( imgType, 'unknown' ) && ~isfield( json, 'Manufacturer' ) || ~strcmp( json.Manufacturer, 'Siemens' )
        warning('Unknown image type. Possibly due to images sourced from non-Siemens MRI')
    end
end

end
% -----------------------------------------------------------------------------
function [img, info] = rescale( img0, info0 )
%RESCALE Rescale image according to NIfTI header info
%      
%     [img, info] = rescale( img, info )
    
    img1 = info0.AdditiveOffset + info0.MultiplicativeScaling .* img0;
   
    %% Check for possible integer overflow/type issues
    if isequaln( img1, double(info0.AdditiveOffset) + double(info0.MultiplicativeScaling) * double(img0) )

        img  = img1;
        info = info0;
        
        %% Update header: both Matlab simplified and original `raw` portions
        info.MultiplicativeScaling = 1;
        info.AdditiveOffset        = 0;
        info.raw.scl_slope         = 1;
        info.raw.scl_inter         = 0;

    else
        warning( ['Aborting default image rescaling to avoid integer overflow.\n'
                  'Original NIfTI header may contain errors.\n '] ) 
    
        img  = img0;
        info = info0;
    end
end
