function [imgs, Hdrs] = load_and_sort_images( List )
%%LOAD_AND_SORT_IMAGES Load and sort image files
%
%     [imgs, Hdrs] = load_and_sort_images( List ) 
%
% Loads a `List` of image files (returned from `imutils.io.find_image_files`)
% and sorts images + headers into 2 cell arrays: 
%
% —`imgs` contains numeric arrays of image data (e.g. voxel values); 
%
% —`Hdrs` contains struct arrays of headers.
%
% When `List` contains images from multiple acquisition series, each
% series becomes a separate element of the returned cells.
%
% Each series array (i.e. cell-element `i`) is sorted according to 
% slice, echo time, measurement number, channel number:
%
% For `imgs{i}`, these dimensions are respectively 3, 4, 5, and 6;
%
% For `Hdrs{i}`, they are respectively 1, 2, 3, and 4.
%
% __NOTE__
% Returns are issued as cell arrays presuming it may be useful to create
% several `Img` objects at once, i.e.: User provides a path to `img.make()`
% containing several (e.g. scan series) image subdirectories, in which case,
% successive series are likely to possess different numbers of images (so,
% arrays `imgs{i+1}` and `Hdrs{i+1}` would be sized differently from their
% preceding elements). 
% For now, loading will be restricted to 1 series (directory) at a time.
% 
% __ETC__ 
%
% See also 
% imutils.io.find_image_files 
    arguments
        List(:,1) struct ;
    end

imgs = [] ;
Hdrs = [] ;

if isempty(List)
    warning('List of image files was empty.') ; 
    return ;
elseif ~isfield( List, 'bytes' ) || ~isfield( List, 'name' )
    error(['Invalid input: List should be a struct returned from ' ...
           'imutils.io.find_image_files()'])
end

%% Extract specific file types (e.g. Dicom) as they will need to be handled differently
filenames = string( [ { List(:).name } ] ) ;

%% Dicoms: 
if any( endsWith( filenames, [".dcm" ".IMA"] ) )
    [imgs, Hdrs] = imutils.io.load_and_sort_dicoms( ...
        List( endsWith( filenames, [".dcm" ".IMA"] ) ) ) ;
end

%% Etc: (if other supported file types are added)

end

% DEV NOTE: 
% Consider +optional argument limitting the data size loaded to memory.
% Not adding now for simplicity + it may not be needed.
% Also note that the Matlab `memory` function can query available memory but
% it currently only works on Windows.
%
% e.g.
%     [imgs, Hdrs] = loadandsortimages( List, nBytesMax ) 
%
% nBytesMax(1,1) {mustBePositive} = 2E9 ;
%    nBytesMax=[2E9]
%      Byte limit on load. The function aborts if the limit is exceeded.
%
%% Check total file size 
% assert( sum( [ List(:).bytes ] ) < nBytesMax, [ 'Byte limit exceeded [%s GB].\n ' ...
%     'Select a new image path or byte limit and run again.' ], num2str( nBytesMax/1E9 ) ) ;
