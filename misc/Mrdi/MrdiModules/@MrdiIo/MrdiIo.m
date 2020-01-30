classdef MrdiIo
%MrdiIo  MR DICOM Image I/O 
% 
% Helper class for writing/reading images, and saving/instantiating Mrdi
% objects to/from file.

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    

properties( Constant = true ) % constant until there is reason to change

    % Image file formats supported for object instantiation
    supportedInputs = [ MrdiIo.dicomExt ] ; % just Dicoms for now. 

    % % Supported output formats for writing objects to file
    % supportedOutputs = { '.dcm' ; '.mrdi' } ; 

end

properties( Constant = true, Hidden = true ) % constant until there is reason to change

    % Known Dicom file extensions 
    dicomExt = { '.dcm' ; '.IMA' } ; 

    % Known Nifti file extensions
    niftiExt = { '.nii' } ;

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
    [List] = findimagefiles( imgPath, fileExt )
    %.....
    [Imgs, Hdrs] = loadandsortimages( imgPath, fileExt, nBytesMax )
    %.....
    [Imgs] = make( varargin )
    %.....
    [] = write( Img, saveDirectory, imgFormat, isSavingSingleNiis )
end
% =========================================================================
% =========================================================================
methods(Static, Hidden=true) 
    %.....
    [Imgs, Hdrs] = loadandsortdicoms( List ) 
end
% =========================================================================
% =========================================================================


end
