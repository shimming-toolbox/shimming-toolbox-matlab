%% SORTIMAGES 
%%
function [] = sortimages( unsortedDicomDir, sortedDicomDir, isCopying )
%SORTIMAGES Arrange+rename unsorted DICOMs into organized subdirectories
%
%% Description 
%
% Siemens 'real-time' image transfer (i.e. via LAN socket) from 
% console to local PC simply dumps all image files into the shared folder
% without organizing acquisitions into subfolders or giving intelligible file
% names. sortimages() creates acquisition and echo subdirectories and renames
% image files according to the DICOM header.
% 
%% Usage
%
%    [] = sortimages( unsortedDicomDir ) 
%    [] = sortimages( unsortedDicomDir, sortedDicomDir ) 
%    [] = sortimages( unsortedDicomDir, sortedDicomDir, isCopying ) 
% 
% If sortedDicomDir is unspecified, a subdirectory ('sorted') is created
% within unsortedDicomDir to contain the sorted images.
%
% isCopying (Boolean) is TRUE by default. Set to 0 to move the files instead
% of copying. Note that the files will still be renamed.
%
%% TODO
%
% * Sort without call to parse_siemens_shadow() (slow!)

DEFAULT_SORTEDDICOMDIR = [ unsortedDicomDir '/sorted/' ] ;
DEFAULT_ISCOPYING      = true ;

%Dicom files inside the directory -----------------------------------------
listOfImages = dir( [ unsortedDicomDir '/*.dcm'] );

if length(listOfImages) == 0
    % try .IMA
    listOfImages = dir( [unsortedDicomDir '/*.IMA'] ) ;
end

nImages = length( listOfImages ) ;

assert( nImages~=0, 'No .dcm or .IMA files found in given directory' ) ;

if nargin < 2 || isempty( sortedDicomDir )
    sortedDicomDir = DEFAULT_SORTEDDICOMDIR ;
end

if ~exist( sortedDicomDir, 'dir') 
    mkdir( sortedDicomDir ); 
end

if nargin < 3 || isempty( isCopying )
    isCopying = DEFAULT_ISCOPYING ;
end

%Read Dicom files header and extract series names and numbers -------------
for iImage = 1 : nImages

    display( ['Sorting... ' num2str(100*iImage/nImages) ' %)' ])
    
    Hdr       = dicominfo( [unsortedDicomDir '/' listOfImages(iImage).name] );
    Hdr.Img   = parse_siemens_shadow( Hdr ) ;  
    seriesDir = [ num2str(Hdr.SeriesNumber) '_' Hdr.SeriesDescription '/' ];

    if ( Hdr.SeriesNumber < 10)
        seriesDir = [ '0' seriesDir ] ;
    end
    
    % Create directories for each series
    if ~exist( [ sortedDicomDir seriesDir ], 'dir' )
        mkdir( [ sortedDicomDir seriesDir ] );
    end
    
    echoDir = [ sortedDicomDir seriesDir 'echo_' num2str( Hdr.EchoTime, 3 ) ] ;
         
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
