classdef MrdiIo
%MrdiIo  MR DICOM Image I/O 
% 
% Class for writing/reading images, and saving/instantiating Mrdi
% objects to/from file.



%% Properties
% 
% * FileExt : Struct of image file extensions 
%       .dicom: known dicom file extensions
%       .nifti: known nifti file extensions
%       .supportedInputs: file types supported for Mrdi object instantiation (restricted now to dicoms)
%
% * searchDir : Search directory for image files [ default: "./" ]
%       Determines whether searchDir subdirectories are included in image search [default: TRUE]
    
    % % File extensions of images to search for [string array, default = [".dcm" ".IMA" ]
    % % (searchFileExt must be a member of MrdiIo.FileExt.supportedInputs, i.e. the default values)
    % searchFileExt {mustBeMember( searchFileExt, [".dcm" ".IMA" ] )} = MrdiIo.FileExt.supportedInputs ;
    %
    % % Byte limit on image loading (loading aborts when exceeded) [default: 2E9]
    % %
    % % NOTE: Matlab memory() function exists to query available memory but
    % % it currently only works for Windows...
    % nBytesMax(1,1) {mustBePositive} = 2*( 1E9 ) ;
    %
    % % Determines whether MrdiIo.make typecasts Mrdi objects to known subclasses [default: TRUE]
    % isConverting(1,1) {mustBeNumericOrLogical} = true ;

properties( Constant = true ) 

    % Struct of image file extensions 
    % dicom: known dicom file extensions
    % nifti: known nifti file extensions
    % supportedInputs: file types supported for Mrdi object instantiation (restricted now to dicoms)
    FileExt = struct( 'dicom', [ ".dcm" ".IMA" ], ...
                      'nifti', [".nii" ], ...
            'supportedInputs', [".dcm" ".IMA" ] ) ;

end

properties

    % Search directory for image files [ default: "./" ]
    searchDir(1,:) {mustBeStringOrChar} = "./" ;
    
    % Determines whether searchDir subdirectories are included in image search [default: TRUE]
    isSearchRecursive(1,1) {mustBeNumericOrLogical} = true ;
    
    % File extensions of images to search for [string array, default = [".dcm" ".IMA" ]
    % (searchFileExt must be a member of MrdiIo.FileExt.supportedInputs, i.e. the default values)
    searchFileExt {mustBeMember( searchFileExt, [".dcm" ".IMA" ] )} = MrdiIo.FileExt.supportedInputs ;
    
    % Byte limit on image loading (loading aborts when exceeded) [default: 2E9]
    %
    % NOTE: Matlab memory() function exists to query available memory but
    % it currently only works for Windows...
    nBytesMax(1,1) {mustBePositive} = 2*( 1E9 ) ;
   
    % Determines whether MrdiIo.make typecasts Mrdi objects to known subclasses [default: TRUE]
    isConverting(1,1) {mustBeNumericOrLogical} = true ;

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Io = MrdiIo( )
    return;
end
% =========================================================================

end
% =========================================================================
% =========================================================================    
methods(Static) 
    %.....
    [Hdr]        = dicominfosiemens( varargin )
    %.....
    [List]       = findimagefiles( searchDir, fileExt, isSearchRecursive )
    %.....
    [Imgs, Hdrs] = loadandsortimages( List, nBytesMax )
    %.....
    [Imgs]       = make( varargin )
    %.....
    [Defaults]   = makedefaults( varargin )
    %.....
    []           = write( Img, saveDirectory, imgFormat, isSavingSingleNiis )
end
% =========================================================================
% =========================================================================
methods(Static, Hidden=true) 
    %.....
    [Imgs, Hdrs] = loadandsortdicoms( List ) 
    %.....
    [isValid]    = assertheadervalidity( Hdr ) 
end
% =========================================================================
% =========================================================================

end
