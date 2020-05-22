function [img,info,json] = read_nii( niiFile )
%READ_NII Load NIfTI image and header; and, when present, the accompanying .json sidecar
%     
%     [img, info]       = read_nii( niiFile )
%     [img, info, json] = read_nii( niiFile )
% 
% The input `niiFile` is the path to the NIfTI image as a string scalar or
% character vector. When called with 2 output arguments, the function is
% equivalent short-hand for
%     
%     info = niftiinfo( niiFile ); img = niftiread( info ) 
%
% When called with the 3rd output argument, the function checks the parent
% folder of `niiFile` for an identically named file but with a .json file extension.
% When such a file is present, the 3rd output is returned as a struct via 
% `json = jsondecode( fileread( jsonFile ) );` otherwise, `json = []`.
%arguments
%    niiFile(1,:) string {mustBeStringScalarOrCharVector, mustBeFile} ;
%end

info = niftiinfo( niiFile ); 
img  = niftiread( info ); 

if nargout < 3
    return;
end

[folder, name] = fileparts( niiFile );
jsonFile       = fullfile( folder, name + ".json" ); 

if isfile( jsonFile )
    json = jsondecode( fileread( jsonFile ) );     
else
    json = [];
end

end
