function [imgs, Hdrs] = load_and_sort_dicoms( List ) 
%LOAD_AND_SORT_DICOMS  Load and sort dicoms from file list
%     
%     [imgs, Hdrs] = load_and_sort_dicoms( List )
%
% Accepts a `List` of dicom file names to load, and sorts the images and headers
% into 2 cell arrays which respectively contain numeric arrays of images, and
% struct arrays of headers. Each acquisition series is assigned to a separate
% element of the returned cells.
%
% `List` is a 1-D struct array of file paths returned from `imutils.io.find_image_files`.
%
% Sorting of the returned arrays occurs according to slice, echo time,
% measurement number, and channel. 
%
% For the image arrays, these dimensions are respectively 3, 4, 5, and 6;
%
% For the headers, they are respectively 1,2,3, and 4.
    arguments
        List(:,1) struct ;
    end

imgs = [] ;
Hdrs = [] ;

%% Load and sort headers
Hdrs    = loaddicomheaders( List ) ;
nSeries = numel( Hdrs ) ; 

Progress = CmdLineProgressBar( [ num2str(nSeries) ' acquisition series found. Sorting headers:' ]);

for iSeries = 1 : nSeries 
    Progress.print( iSeries, nSeries ) ;
    Hdrs{iSeries} = sortdicomheaders( Hdrs{iSeries} ) ;   
end

%% Augment data struct with Siemens private header courtesy of parse_siemens_shadow()
if strcmp( Hdrs{1}(1).Manufacturer, 'SIEMENS' )
    Progress = CmdLineProgressBar( [ 'Parsing Siemens private header for 1st image in series: ' ] );

    for iSeries = 1 : nSeries 
        Progress.print( iSeries, nSeries ) ;
        % NOTE: Ideally, dicominfosiemens() would be used for every DICOM file,
        % however, it is rather slow. Compromise: parse the private metadata
        % for the first image in the series only, as much of it will be
        % constant across the series (i.e. across DICOM files).
        ExtendedHdr             = dicominfosiemens( Hdrs{iSeries}(1) ) ;
        Hdrs{iSeries}(1).Img    = ExtendedHdr.Img ; 
        Hdrs{iSeries}(1).Ser    = ExtendedHdr.Ser ; 
        Hdrs{iSeries}(1).MrProt = ExtendedHdr.MrProt ; 
    end
end

imgs = loaddicomimages( Hdrs ) ;

end

% -----
%% Internal functions

function [ Hdrs ] = loaddicomheaders( List )
%LOADDICOMHEADERS  Return cell array of DICOM Hdr structs
%     
%     Hdrs = loaddicomheaders( List )
% 
% `loaddicomheaders` wraps to `dicominfo` for each element in `List` (a struct
% (e.g. returned from `dir()`) whereby fields .name and .folder describe the
% DICOM files to be loaded. 
%
% Headers from DICOMs belonging to the same series (indicated by the field
% .SeriesInstanceUID) are grouped into the same struct array, which becomes an
% element of the returned cell array Hdrs. 
% 
% Returned Hdrs is then a 1-D cell array, with as many elements as there are
% distinct acquisition series. 

    nImg = length( List ) ;
    
    Hdrs               = cell( 1, 1 ) ; % cell containing dicom header struct arrays
    nSeries            = 0 ; % number of series found
    seriesInstanceUIDs = "" ; % string array holding SeriesInstanceUID's
    Progress           = CmdLineProgressBar([ 'Loading image headers: ' ]);

    for iImg = 1 : nImg 
        
        Progress.print( iImg, nImg ) ;
        
        Hdr     = dicominfo( [List(iImg).folder filesep List(iImg).name] ) ;
        % Check if new Hdr corresponds to previously loaded series; 
        % if not, append new entry.
        iSeries = find( Hdr.SeriesInstanceUID == seriesInstanceUIDs ) ;
        
        if isempty( iSeries )
            nSeries = nSeries + 1 ;
            iSeries = nSeries ;
            seriesInstanceUIDs(nSeries) = string( Hdr.SeriesInstanceUID ) ;
            Hdrs{iSeries} = Hdr ;
        else
            Hdrs{iSeries}(end+1) = Hdr ;
        end
    end

end % loaddicomheaders

function [ HdrsSorted ] = sortdicomheaders( Hdrs )
%SORTDICOMHEADERS Reorders and reshapes 1-d DICOM header struct array 
%     
%     HdrsSorted = sortdicomheaders( Hdrs )
% 
% Input `Hdrs` should be a 1-d struct-array of DICOM Hdrs, in no particular
% order, where all elements (i.e. every included DICOM header) belong to the
% same acquisition series. 
%
% Output `HdrsSorted` has same number of elements as the input, but is
% reordered and arrayed such that its dimensions refer, respectively, to 
% ( Slice, Echo, Acquisition, Channel )
    arguments
        Hdrs {mustBeValidHeader} ;
    end

    nImgInSeries      = length( Hdrs ) ;

    nRows             = unique( [ Hdrs(:).Rows ] ) ;
    nColumns          = unique( [ Hdrs(:).Columns ] ) ;

    sliceLocations    = unique( [ Hdrs(:).SliceLocation ] ) ;
    nSlices           = length( sliceLocations ) ;

    echoTimes         = unique( [ Hdrs(:).EchoTime ] ) ;
    nEchoes           = length( echoTimes ) ;

    acquisitionNumber = unique( [ Hdrs(:).AcquisitionNumber ] ) ;
    nAcquisitions     = length( acquisitionNumber ) ;

    channels          = unique( string( { Hdrs(:).Private_0051_100f } ) ) ;
    nChannels         = length( channels ) ;

    % Initialize
    HdrsSorted = Hdrs(1) ; 

    for iImg = 1 : nImgInSeries 
        
       iSlice        = find( Hdrs(iImg).SliceLocation == sliceLocations ) ;
       iEcho         = find( Hdrs(iImg).EchoTime == echoTimes ) ;
       iAcquisition  = find( Hdrs(iImg).AcquisitionNumber == acquisitionNumber ) ;
       iChannel      = find( Hdrs(iImg).Private_0051_100f == channels ) ;

       HdrsSorted(iSlice, iEcho, iAcquisition, iChannel) = Hdrs(iImg) ;
    
    end

end % sortdicomheaders

%% -----
function [ imgs ] = loaddicomimages( Hdrs )
%LOADDICOMIMAGES  Return cell array of dicom images 
%     
%     [ imgs ] = loaddicomimages( Hdrs )
%
% `Hdrs` is a cell array, wherein each element is a struct array of dicom headers.
% Return cell array `imgs` is the same size as `Hdrs`, and each element is an image
% array, its size a function of the corresponding `Hdrs{i}` entry

    nSeries  = numel( Hdrs ) ;
    imgs     = cell( 1, nSeries ) ;

    Progress = CmdLineProgressBar( [ 'Loading images: ' ] );

    for iSeries = 1 : nSeries 

        Progress.print( iSeries, nSeries ) ;

        [ nSlices, nEchoes, nMeasurements, nChannels ] = size( Hdrs{iSeries} ) ;

        for iSlice = 1 : nSlices
            for iEcho = 1 : nEchoes
                for iMeasurement = 1 : nMeasurements
                    for iChannel = 1 : nChannels
                        img(:,:,iSlice,iEcho,iMeasurement,iChannel) = ...
                            dicomread( Hdrs{iSeries}(iSlice, iEcho, iMeasurement, iChannel) ) ;
                    end
                end
            end
        end
    
        imgs{iSeries} = img ;
        clear img ;

    end

end % loaddicomimages

function [] = mustBeValidHeader( Hdrs )
%mustBeValidHeader Issues an error if an image header (struct) is invalid 
%     
%     [] = mustBeValidHeader( Hdrs )
%
% `Hdrs` can be a DICOM header struct as returned by `dicominfo`. If it does
% not contain the necessary data fields to permit the construction of a
% img-object, then an error is thrown.
    if nargin ~= 1 || ~isstruct( Hdrs )
        error('Function requires single input argument: Hdrs struct array.' ) ;
    end
    
    %% Check fields exist
    REQUIRED_FIELDS = { 'SeriesInstanceUID'; 
                        'Rows'; 
                        'Columns'; 
                        'PixelSpacing' ; 
                        'ImageOrientationPatient' ; 
                        'ImagePositionPatient' ; 
                      } ;
    
    
    if any( ~isfield( Hdrs, REQUIRED_FIELDS ) )
        missingFields = { REQUIRED_FIELDS{ ~isfield( Hdrs, REQUIRED_FIELDS ) } } ;
        error( [ 'Required Hdr fields not defined:\n .' strjoin( missingFields, '\n .' ) ], '%s' ) ; 
    end
    
    assert( isfield( Hdrs, 'SliceThickness' ) | isfield( Hdrs, 'SpacingBetweenSlices' ), ...
        'Expected Hdr entry for at least one of the fields .SliceThickness or .SpacingBetweenSlices' ) ;

    %% Check field entries that should possess a unique value
    assert( length( unique( string( { Hdrs(:).SeriesInstanceUID } ) ) ) == 1, ...
        'All elements of the DICOM Headers struct array must belong to the same acquisition series' ) ;

    assert( length(unique( [ Hdrs(:).Rows ] )) == 1, ...
        'Expected unique Hdr entry for .Rows' ) ;
    
    assert( length( unique( [ Hdrs(:).Columns ] ) ) == 1, ...
        'Expected unique Hdr entry for .Columns' ) ;

    % NOTE/TODO
    % The following conditions need to be met for instantiating an img.Grid
    % object, however, they may not pass for acquisitions with multiple slice
    % groups (e.g. localizer)...
    
    % pixelSpacing = [ Hdrs(:).PixelSpacing ] ;
    %
    % assert( ( length( unique( pixelSpacing(1,:) ) ) == 1 ) || ...
    %         ( length( unique( pixelSpacing(2,:) ) ) == 1 ), ...
    %     'Expected unique Hdr entry for both .PixelSpacing(1) and .PixelSpacing(2)' ) ;
    %
    % if myisfield( Hdrs, 'SliceThickness' )
    %     assert( length( unique( [ Hdrs(:).SliceThickness ] ) ) == 1, ...
    %         'Expected unique Hdr entry for .SliceThickness' ) ;
    % end
    %
    % if myisfield( Hdrs, 'SpacingBetweenSlices' )
    %     assert( length( unique( [ Hdrs(:).SpacingBetweenSlices ] ) ) == 1, ...
    %         'Expected unique Hdr entry for .SpacingBetweenSlices' ) ;
    % end

end % mustBeValidHeader

%% -----
function [img, Hdr] = reshapemosaic( img, Hdr )
%RESHAPEMOSAIC Reshape Siemens mosaic into volume array and remove padded zeros
%
% Adapted from [dicm2nii](http://www.mathworks.com/matlabcentral/fileexchange/42997/) by xiangrui.li@gmail.com 

error('DEPRECATED')
assert( ~isempty(strfind( Img.Hdr.ImageType, 'MOSAIC' ) ), 'Corrupt image header?' ) ;       

nImgPerLine = ceil( sqrt( Img.getnumberofslices() ) ); % always nImgPerLine x nImgPerLine tiles

nRows    = size(Img.img, 1) / nImgPerLine; 
nColumns = size(Img.img, 2) / nImgPerLine; 
nEchoes  = size(Img.img, 4) ;  
nVolumes = size(Img.img, 5) ;

img = zeros([nRows nColumns Img.getnumberofslices() nEchoes nVolumes], class(Img.img));

for iImg = 1 : Img.getnumberofslices()

    % 2nd slice is tile(1,2)
    r = floor((iImg-1)/nImgPerLine) * nRows + (1:nRows);     
    c = mod(iImg-1, nImgPerLine) * nColumns + (1:nColumns);

    img(:, :, iImg, :) = Img.img(r, c, :);
end

Img.img = img ;


% -----
% Update header

% -----
% Correct Hdr.ImagePositionPatient 
%   see: http://nipy.org/nibabel/dicom/dicom_mosaic.html

%-------
% Rotation and Scaling matrix: RS  
RS = Img.rotationMatrix * diag(Img.voxelSpacing) ;

Img.Hdr.ImagePositionPatient = Img.Hdr.ImagePositionPatient + ...
    RS(:,1)*(double(Img.Hdr.Rows) - nRows)/2 + ...
    RS(:,2)*(double(Img.Hdr.Columns) - nColumns)/2 ;


Img.Hdr.Rows    = uint16(nRows) ;
Img.Hdr.Columns = uint16(nColumns) ;

tmp = strfind( Img.Hdr.ImageType, '\MOSAIC' ) ;

if ~isempty(tmp) && (tmp > 1)
    Img.Hdr.ImageType = Img.Hdr.ImageType(1:tmp-1) ;
else
    Img.Hdr.ImageType = '' ;    
end

end %reshapemosaic
