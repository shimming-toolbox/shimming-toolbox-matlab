classdef (Sealed, Hidden) Maker < handle
%img.Maker Loads and sorts image files for the construction of Img objects 
%
% img.Maker facilitates the construction of image objects by
% finding, organizing, and loading image files.
%
% Further, (TODO) it will handle typecasting of generic images into specific subtypes when appropriate
% (e.g. if files are phase images, `img.Maker` should convert them).
% 
% __DEV NOTE__
% This class should only ever need to be accessed via the `img.make()`
% function: the public interface to create image objects). 
%
% As such, class methods could technically all be local functions within
% `img.make`; however, doing so would *make* that single file extremely
% complicated. Delegating and compartmentalizing the various tasks involved
% with loading image files to the methods of this class seems likely to be more
% maintainable. 
%
% Also, it constitutes a step toward keeping the interface and
% implementation [separate.](https://oop.tech-academy.co.uk/interface-vs-implementation/)

    
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
    FileExt = struct( 'dicom', [ ".dcm" ".IMA" ], ...
                      'nifti', ".nii" ) ;

    % File types (extensions) supported for Img object instantiation (restricted now to dicoms)
    supportedInputs = [".dcm" ".IMA"] ;

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Io = Maker( )
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
    [data, Hdrs] = loadandsortdicoms( List ) 
    %.....
    [data, Hdrs] = loadandsortimages( List, nBytesMax )
end
% =========================================================================
% =========================================================================

end
