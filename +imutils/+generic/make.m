function [Imgs] = make( searchDir ) 
%MAKE Constructs an `Img` object from file 
%     
%     [Imgs] = make( searchDir ) 
%
% Initializes an `Img` image object with the image data found in `searchDir`. 
%
% __INPUTS__
%   
%   `searchDir` 
%     directory to search for images (string or char array)
    arguments
        searchDir(1,:) string {mustBeFileOrFolder};
    end
   
Maker = img.Maker;

%% Find image files 
List = Maker.findimagefiles( searchDir );

%% Read in files
[imgs, Hdrs] = Maker.loadandsortimages( List );

% %% ----- 
% % Initialize Img objects 
% for iSeries = 1 : numel( imgs )
%     Imgs{iImg} = MaRdI( imgs{iSeries}, Hdrs{iSeries} ) ;
% end

% When the image type is known to correspond to a specific Img
% child class, the object should be typecast accordingly.
%
% isConverting(1,1) {mustBeNumericOrLogical} = true ;
% Determines whether MrdiIo.make typecasts Mrdi objects to known subclasses [default: TRUE]


end
