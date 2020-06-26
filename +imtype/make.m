function [img] = make( searchDir ) 
%MAKE Constructs an `img` object from file 
%     
%     img = make( searchDir ) 
%
% Initializes an `img` image object with the image data found in `searchDir`. 
%
% __INPUTS__
%   
%   `searchDir` 
%     directory to search for images (string or char array)
    arguments
        searchDir(1,:) string {mustBeFileOrFolder};
    end

%% Find image files 
List = imutils.io.find_image_files( searchDir );

%% Read in files
[imgs, Hdrs] = imutils.io.load_and_sort_images( List );

% %% ----- 
% % Initialize Img objects 
% for iSeries = 1 : numel( imgs )
%     Imgs{iImg} = Img( imgs{iSeries}, Hdrs{iSeries} ) ;
% end

% When the image type is known to correspond to a specific Img
% child class, the object should be typecast accordingly.
%
% isConverting(1,1) {mustBeNumericOrLogical} = true ;
% Determines whether MrdiIo.make typecasts Mrdi objects to known subclasses [default: TRUE]


end
