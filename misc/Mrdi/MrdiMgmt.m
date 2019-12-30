classdef (Abstract) MrdiMgmt < handle
%MrdiMgmt MR DICOM Image File Management
% 
% Member methods pertain to managing (finding, sorting) .dcm/.IMA files
% 
% e.g.
%   sortimages()
%   findimages()
% 
% For documentation, type doc MrdiMgmt
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = MrdiMgmt( imgPath )

end
% =========================================================================    

end
% =========================================================================
% =========================================================================    
methods(Static)
% =========================================================================    
function list = findimages( imgDir )
%FINDIMAGES Lists .dcm OR .IMA files in a directory & 'echo_*' subdirectories
%
% list = FINDIMAGES( imageDirectory ) 
%
% Returns a cell array.

imgDir       = [ imgDir '/' ] ;

ListSubdirs  = dir( [ imgDir 'echo*'] );
nEchoSubdirs = length( ListSubdirs ) ;

if nEchoSubdirs > 0 
    
    List = finddcmorima( [imgDir ListSubdirs(1).name] ) ;
    nImgPerEcho = length(List) ;

    list = cell( nImgPerEcho*nEchoSubdirs, 1 ) ;

    for iImg = 1 : nImgPerEcho 
        list{iImg} = [ imgDir ListSubdirs(1).name '/' List(iImg).name ] ;
    end

    for iEcho = 2 : nEchoSubdirs
        ListiSubdir = finddcmorima( [imgDir ListSubdirs(iEcho).name] ) ;
        assert( length(ListiSubdir) == nImgPerEcho, 'Each echo subdirectory should contain the same number of images' ) ; 

        for iImg = 1 : nImgPerEcho
            list{ (iEcho-1)*nImgPerEcho + iImg} = [ imgDir ListSubdirs(iEcho).name '/' ListiSubdir(iImg).name ] ;
        end
    end

else
    List = finddcmorima( imgDir ) ;
    nImg = length(List) ;
    list = cell( nImg, 1 ) ;

    for iImg = 1 : nImg
        list{iImg} = [ imgDir List(iImg).name ] ;
    end
end


function List = finddcmorima( imgDir ) 
%Find .dcm OR .IMA files in imgDir
List   = dir( [ imgDir '/*.dcm'] );

if length(List) == 0 % try .IMA
    List = dir( [imgDir '/*.IMA'] ) ;
end

nImages = length( List ) ; 
assert( nImages ~= 0, 'No .dcm or .IMA files found in given directory' ) ;

end

end
% =========================================================================
function [] = sortimages( unsortedDicomDir, sortedDicomDir, isCopying )
%SORTIMAGES Arrange+rename unsorted DICOMs into organized subdirectories
%
% Description 
%
%   Siemens 'real-time' image transfer (i.e. via LAN socket) from 
%   console to local PC simply dumps all image files into the shared folder
%   without organizing acquisitions into subfolders or giving intelligible file
%   names. SORTIMAGES() creates acquisition and echo subdirectories and renames
%   image files according to the DICOM header.
% 
% Usage
%
% [] = SORTDATA( unsortedDicomDir ) 
% [] = SORTDATA( unsortedDicomDir, sortedDicomDir ) 
% [] = SORTDATA( unsortedDicomDir, sortedDicomDir, isCopying ) 
% 
% If sortedDicomDir is unspecified, a subdirectory ('sorted') is created
% within unsortedDicomDir to contain the sorted images.
%
% isCopying (Boolean) is TRUE by default. Set to 0 to move the files instead
% of copying. Note that the files will still be renamed.

%TODO : Sort without call to parse_siemens_shadow() (slow!)

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
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static, Hidden=true) 
% Possibly deprecated functions which may or may not be useful
% =========================================================================    
function fullDir = getfulldir( dataLoadDir, iDir )
%GETFULLDIR
% 
% fullDir = GETFULLDIR( parentDir, index ) 
%
%   Searches parentDir[ectory] for subdirectory beginning with index- (e.g.
%   .../dataLoadDir/[index]-* ) to return its full path.
%
%   (Useful for rapidly initializing MrdiImg with a dicom folder.)

if nargin < 2
    error('Function requires 2 input arguments: parent directory and index.')
end

if isnumeric(iDir)
    iDir = num2str( iDir ) ;
end

if length( iDir ) == 1
    iDir = ['0' iDir] ;
end

if ~strcmp( dataLoadDir(end), '/' )
    dataLoadDir(end+1) = '/' ;
end

Tmp = dir( [ dataLoadDir iDir '-*'] ) ;
fldrName = Tmp.name ;

[fullDir,~,~] = fileparts( [dataLoadDir fldrName '/'] ) ;

fullDir(end+1) = '/' ;

end
% =========================================================================
function [ studyDirs ] = tablestudy( sortedDicomDir )
%TABLESTUDY 
%
% Returns a cell array ( studyDirs ) pertaining to the study directory input
% ( sortedDicomDir ) where each element in the second column is a MaRdI-loadable 
% images series. (The 1st column is merely the row index.)
%
% e.g. Protocol to load MaRdI-object :
%
%   % omit semi-colon to display the full cell array (i.e. with the row indices)
%   [ studyDirs ] = MaRdI.tablestudy( sortedDicomDir ) 
%
%   % determine the row index of the series you want to load (e.g. 10):
%   Img = MaRdI( studyDirs{ 10, 2 } ) ;

assert( nargin == 1, 'Function requires sortedDicomDirectory as input argument.' ) ;

if ~strcmp( sortedDicomDir(end), '/' ) 
    sortedDicomDir(end+1) = '/' ;
end

studyDirs = cell( 0 ) ;

Tmp      = dir( [ sortedDicomDir ] );
Tmp      = Tmp( 3:end ) ; % ignore self ('.') and parent ('..') dirs
nEntries = length( Tmp ) ;

for iEntry = 1 : nEntries 

   if Tmp( iEntry ).isdir
   
       tmpSeriesSubdir = [ Tmp( iEntry ).name '/'] ;
    
        TmpEchoSubdirs = dir( [ sortedDicomDir tmpSeriesSubdir 'echo*' ] ) ;
        nEchoSubdirs   = length( TmpEchoSubdirs ) ;

        if nEchoSubdirs ~= 0

            for iEchoSubdir = 1 : nEchoSubdirs

                studyDirs{end+1, 2} = strcat( sortedDicomDir, tmpSeriesSubdir, TmpEchoSubdirs(iEchoSubdir).name )  ;
                iSeries = size( studyDirs, 1 ) ;
                studyDirs{ iSeries, 1 } = iSeries ;
            end

        % check if tmpSeriesSubdir itself contains images
        elseif length( dir( [ sortedDicomDir tmpSeriesSubdir '/*.dcm'] ) ) ~= 0 || ...
                length( dir( [ sortedDicomDir tmpSeriesSubdir '/*.IMA'] ) ) ~= 0 

           studyDirs{end+1, 2} = strcat( sortedDicomDir,tmpSeriesSubdir )  ;
            iSeries = size( studyDirs, 1 ) ;
            studyDirs{ iSeries, 1 } = iSeries ;

        end

   end
   
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
