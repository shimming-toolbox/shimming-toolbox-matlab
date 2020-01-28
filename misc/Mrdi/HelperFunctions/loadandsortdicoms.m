function [Imgs, Hdrs] = loadandsortdicoms( imgPath, fileExt, nBytesMax ) 
%LOADANDSORTDICOMS  Load and sort dicoms
% 
% LOADANDSORTDICOMS searches a directory for dicom images and sorts the
% images+headers into 2 cell arrays which respective contain numeric arrays of
% images, and struct arrays of headers. Each acquisition series is assigned to 
% a separate element of the returned cells.
%
% Sorting of the returned arrays occurs according to slice, echo time,
% measurement number, and channel. 
% For the image arrays, these dimensions are respectively 3, 4, 5, and 6;
% For the headers, they are respectively 1,2,3, and 4.
%
% [Imgs, Hdrs] = loadandsortdicoms( imgPath, fileExt, nBytesMax )
%
% Inputs
%
%   imgPath [default: './' ]
%       directory to search for images (char)
%
%   fileExt [default: { '.dcm' ; '.IMA' }]
%       file extensions of images (char or char-containing cell array)
%       (fileExt must be a member of the given default values)
%
%   nBytesMax [default: 2 * 1E9]
%       byte limit on load. Function aborts if limit is exceeded.
%       (Note: matlab 'memory' function exists to query available memory but
%       currently only works on Windows)
    arguments
        imgPath(1,:) char = pwd ;
        fileExt {mustBeMember( fileExt, {'.dcm' ; '.IMA'} )} = {'.dcm' ; '.IMA'} ;
        nBytesMax(1,1) {mustBePositive} = 2*( 1E9 ) ;
    end

Imgs = [] ;
Hdrs = [] ;

%% ----- 
% Find files
display( [ 'Searching for images in ' imgPath ] ) ; 

switch class( fileExt )
    case 'char'
        List = dir( [imgPath filesep '**' filesep '*' fileExt ] ) ;
  
    case 'cell'
        List = dir( [imgPath filesep '**' filesep '*' fileExt{1} ] ) ;

        for iExt = 2 : length( fileExt )
            List_iExt = dir( [imgPath filesep '**' filesep '*' fileExt{iExt} ] ) ;
            List(end+1:end+length(List_iExt)) = List_iExt ;
        end
    
    otherwise
        error('Invalid input: fileExt') ;
end

nImg   = length( List ) ;
nBytes = sum( [ List(:).bytes ] ) ;

display( [ num2str(nImg) ' images found [totaling ' num2str( nBytes ) ' bytes].' ] ) ; 

if isempty(List)
   return ;
end

assert( nBytes < nBytesMax, [ 'Byte limit exceeded [%s GB]. \n' ...
    'Select a different image path or byte limit and run again.' ], num2str( nBytesMax/1E9 ) ) ;

%% -----
% Load and sort headers
Hdrs    = loaddicomheaders( List ) ;
nSeries = numel(Hdrs) ; 

Progress = CmdLineProgressBar( [ num2str(nSeries) ' acquisition series found. Sorting headers:' ]);

for iSeries = 1 : nSeries 
    Progress.print( iSeries, nSeries ) ;
    Hdrs{iSeries} =sortdicomheaders( Hdrs{iSeries} ) ;   
end

%% -----
% Augment data struct with Siemens private header courtesy of parse_siemens_shadow()
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

%% -----
Imgs = loaddicomimages( Hdrs ) ;
return ;


%% =========================================================================
% Internal functions
%% =========================================================================
function [ Hdrs ] = loaddicomheaders( List )
%LOADDICOMHEADERS  Return cell array of DICOM Hdr structs
% 
% Hdrs = LOADDICOMHEADERS( List )
%
% List is a struct (e.g. returned from DIR()) whereby fields .name and .folder 
% describe the DICOM files to be loaded.
% 
% LOADDICOMHEADERS wraps to DICOMINFO() for each element in List. Headers from
% DICOMs belonging to the same series (as indicated by the field
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
        % Check if new Hdr corresponds to previously loaded series, and create new
        % entry in Hdrs if not.
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
%% =========================================================================
function [ HdrsSorted ] = sortdicomheaders( Hdrs )
%SORTDICOMHEADERS  Reorders+reshapes 1-d DICOM header struct array 
%
% HdrsSorted = SORTDICOMHEADERS( Hdrs )
% 
% Input Hdrs should be a 1-d struct-array of DICOM Hdrs, in no particular
% order, where all elements (i.e. every included DICOM header) belongs to the
% same acquisition series. 
%
% HdrsSorted has same number of elements as the input, but is reordered and
% arrayed such that its dimensions refer, respectively, to 
% ( Slice, Echo, Acquisition, Channel )

    if length( unique( string( { Hdrs(:).SeriesInstanceUID } ) ) ) ~= 1
        error('All elements of the DICOM Headers struct array must belong to the same acquisition series' ) ;
    end

    nImgInSeries  = length( Hdrs ) ;

    nRows = unique( [ Hdrs(:).Rows ] ) ; 
    assert( length(nRows) == 1, 'Expected single Hdr entry for .Rows' ) ;
    
    nColumns = unique( [ Hdrs(:).Columns ] ) ; 
    assert( length(nColumns) == 1, 'Expected single Hdr entry for .Columns' ) ;
    
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
%% =========================================================================
function [ Imgs ] = loaddicomimages( Hdrs )
%LOADDICOMIMAGES  Return cell array of dicom images 
% 
% [ Imgs ] = loaddicomimages( Hdrs )
%
% Hdrs is a cell array, wherein each element is a struct array of dicom headers.
% Return cell array Imgs is the same size as Hdrs, and each element is an image
% array, its size a function of the corresponding Hdrs{i} entry

    nSeries  = numel( Hdrs ) ;
    Imgs     = cell( 1, nSeries ) ;

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
    
        Imgs{iSeries} = img ;
        clear img ;

    end

end % loaddicomimages
%% =========================================================================

end
