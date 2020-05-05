function [] = sortdicoms( unsortedDir, sortedDir, isCopying )
%SORTDICOMS Arrange unsorted DICOMs into subdirectories according to acquisition series
%     
%     [] = sortdicoms( unsortedDir ) 
%     [] = sortdicoms( unsortedDir, sortedDir ) 
%     [] = sortdicoms( unsortedDir, sortedDir, isCopying ) 
%
% Copies and renames .dcm or .IMA images within `unsortedDir` into
% subdirectories based on the image headers. 
%
% `sortedDir` is the destination parent folder where series-/scan-specific
% subdirectories will be created if they do not already exist.  When called
% without a second argument, `sortedDir = [unsortedDir '/sorted']` by
% default.
%
% To move rather than copy the files, call the function with the third argument
% `isCopying = 0` (or `false`).
%
% __NOTE__
% This function is used to organize images transfered via Siemens 'real-time'
% socket protocol, by which DICOMs with abstruse filenames are transfered
% immediately upon reconstruction into a single directory (irrespective of
% acquisition series/scan). 

DEFAULT_SORTEDDIR = [ unsortedDir '/sorted/' ] ;
DEFAULT_ISCOPYING = true ;

%Dicom files inside the directory -----------------------------------------
listOfImages = dir( [ unsortedDir '/*.dcm'] );

if length(listOfImages) == 0
    % try .IMA
    listOfImages = dir( [unsortedDir '/*.IMA'] ) ;
end

nImages = length( listOfImages ) ;

assert( nImages~=0, 'No .dcm or .IMA files found in given directory' ) ;

if nargin < 2 || isempty( sortedDir )
    sortedDir = DEFAULT_SORTEDDIR ;
end

if ~exist( sortedDir, 'dir') 
    mkdir( sortedDir ); 
end

if nargin < 3 || isempty( isCopying )
    isCopying = DEFAULT_ISCOPYING ;
end

%Read Dicom files header and extract series names and numbers -------------
for iImage = 1 : nImages

    display( ['Sorting... ' num2str(100*iImage/nImages) ' %)' ])
    
    Hdr       = dicominfo( [unsortedDir '/' listOfImages(iImage).name] );
    Hdr.Img   = parse_siemens_shadow( Hdr ) ;
    seriesDir = [ num2str(Hdr.SeriesNumber) '_' Hdr.SeriesDescription '/' ];

    if ( Hdr.SeriesNumber < 10)
        seriesDir = [ '0' seriesDir ] ;
    end
    
    % Create directories for each series
    if ~exist( [ sortedDir seriesDir ], 'dir' )
        mkdir( [ sortedDir seriesDir ] );
    end
    
    echoDir = [ sortedDir seriesDir 'echo_' num2str( Hdr.EchoTime, 3 ) ] ;
         
    if ~exist( echoDir, 'dir' )
        mkdir( echoDir ) ;
    end
    
    iSlice   = Hdr.Img.ProtocolSliceNumber + 1 ;
    sliceStr = num2str( iSlice ) ;

    acqStr   = num2str( Hdr.AcquisitionNumber ) ;

    for ord = 3 : -1 : 1
        if iSlice < 10^ord
            sliceStr = [ '0' sliceStr ] ;
             acqStr  = [ '0' acqStr ] ;
        end
    end

    [~, ~, ext]    = fileparts( Hdr.Filename ) ;
    sortedFilename = fullfile( echoDir, strcat( Hdr.PatientName.FamilyName, '-', ...
                Hdr.Img.ImaCoilString, '-', sliceStr, '-', acqStr, ext ) ) ;
    
    if isCopying
        copyfile( Hdr.Filename, sortedFilename{1} ) ;
    else
        movefile( Hdr.Filename, sortedFilename{1} ) ;
    end

end
    
end
