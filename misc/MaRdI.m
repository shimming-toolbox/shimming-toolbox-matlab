classdef MaRdI < MrdiUtil & MrdiIo & MrdiMgmt  
%MaRdI Ma(t)-R-dI(com)
%
%   Dicom into Matlab for Siemens MRI data
%
% .......
%   
% Usage
%
%   Img = MaRdI( imgPath )
% 
%   where imgPath is the path to a single .dcm image OR a directory containing
%   the .dcm or .IMA images.
%
% Img contains public properties:
%
%   .img
%       Array of images: 
%       If multi-echo, the echo index is along the 4th dimension;
%       If multiple volume repetitions, the measurement index is along the 5th dimension.
%
%   .Aux
%       Aux-objects: auxiliary measurements (e.g. respiratory ProbeTracking)
%       By default .Aux is empty. To fill it, call Img.associateaux( Aux ) with a valid
%       Aux object. For more info, type: help MaRdI.associateaux( )
% 
% In addition to the read-only properties
%
%    .Hdr 
%       The full Siemens DICOM header corresponding to Img.img(:,:,1,1,1) 
%
%    .Hdrs 
%       Cell array of (truncated) DICOM headers courtesy of dicominfo().
%       (One entry for every image)
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

%% ========================================================================
%
% ..... 
% WRITE()
%   Saving as DICOM (and, by extension, NifTI) not rigorously/fully tested
% 
%%

% =========================================================================
% =========================================================================
properties
    img ; 
    Aux ;
end

properties(SetAccess={?MaRdI, ?MrdiIo, ?MrdiUtil, ?MrdiProc})
    Hdr ; % full Siemens DICOM header of 1st img (i.e. Img.img(:,:,1) )
    Hdrs ; % cell array of (truncated) DICOM headers courtesy of dicominfo()
    Ref ; % Reference properties - prior to manipulation
end

% properties(SetAccess=protected, Hidden = true)
%     Ref ; % Reference properties - prior to manipulation
% end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = MaRdI( imgPath )

Img.img  = [] ;
Img.Hdr  = [] ;
Img.Hdrs = [] ;
Img.Aux  = [] ;

if nargin == 1 
    
    if ~exist( imgPath )
        error( 'DICOM images not found. Check input path is valid.' ) ;
        return;
    
    elseif ( exist( imgPath ) == 2 ) % input is a single file (e.g. DICOM) 
        
        [~, ~, ext] = fileparts( imgPath ) ;
        
        assert( strcmp( ext, '.dcm' ) || strcmp( ext, '.IMA' ), ...
            'Input must be a path string leading to a single dicom image OR to a dicom-containing folder' )

        Img.img = double( dicomread( imgPath ) ) ;
        Img.Hdr = MaRdI.dicominfosiemens( imgPath ) ;

        %  Add/replace a few Hdr entries:
        
        if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
            Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
        end
       
        Img.setslicenormalvector() ;
    
    elseif ( exist( imgPath  ) == 7 ) % input is a directory (e.g. containing DICOMs)
        
        imgDir = [ imgPath '/' ] ; 

        imgList = MaRdI.findimages( imgDir ) ;
        nImages = length( imgList ) ;

        if ~nImages
            return;
        end
        
        % Read protocol info from 1st loaded image (arbitrary):
        SpecialHdr = MaRdI.dicominfosiemens( imgList{1} ) ;

        nRows      = SpecialHdr.Height ;
        nColumns   = SpecialHdr.Width ;
        
        if ~myisfield( SpecialHdr.MrProt, 'lRepetitions' ) 
            nVolumes = 1 ;
        else
            nVolumes = SpecialHdr.MrProt.lRepetitions + 1 ;
        end
        
        % containers for a few terms varying between dicoms:
        sliceLocations = zeros( nImages, 1 ) ; % [units: mm]
        echoTimes      = zeros( nImages, 1 ) ; % [units: ms]

        % Image headers, not yet organized:
        RawHdrs        = cell( nImages, 1 ) ; 
        
        % Read all the image headers and structure the MaRdI object accordingly:
        %
        % NOTE: It may well be possible to determine all the necessary info from
        % the complete Siemens header (e.g. SpecialHdr) of a single image;
        % practically, however, the following is easier since the header
        % information is abstrusely defined for some sequences.
        for iImg = 1 : nImages
            % using dicominfo() here as parse-siemens-shadow() takes much longer
            RawHdrs{ iImg }      = dicominfo( imgList{iImg} ) ;

            sliceLocations(iImg) = RawHdrs{ iImg }.SliceLocation ;
            echoTimes(iImg)      = RawHdrs{ iImg }.EchoTime ;
        end

        sliceLocations = sort( unique( sliceLocations ), 1, 'ascend' ) ; 
        nSlices        = length( sliceLocations ) ;
        
        echoTimes      = sort( unique( echoTimes ), 1, 'ascend' ) ; 
        nEchoes        = length( echoTimes ) ; 
        
        Img.img = zeros( nRows, nColumns, nSlices, nEchoes, nVolumes ) ;
        
        % copy of RawHdrs to be reorganized:
        Img.Hdrs = cell( nSlices, nEchoes, nVolumes ) ;
        
        for iImg = 1 : nImages

            iHdr    = RawHdrs{ iImg } ;

            iVolume = iHdr.AcquisitionNumber ;
            iSlice  = find( iHdr.SliceLocation == sliceLocations ) ;
            iEcho   = find( iHdr.EchoTime == echoTimes ) ;

            Img.Hdrs{ iSlice, iEcho, iVolume } = iHdr ;

            Img.img( :, :, iSlice, iEcho, iVolume ) = dicomread( imgList{iImg} ) ;

        end
        
        % Save the complete 1st Hdr
        Img.Hdr = MaRdI.dicominfosiemens( Img.Hdrs{1}.Filename ) ;
            
        Img.rescaleimg() ;

        if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
            Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
        end
        
        if ~isempty( strfind( Img.Hdr.ImageType, 'MOSAIC' ) ) 
            Img = reshapemosaic( Img ) ;
        end

        Img.setslicenormalvector() ;

    end
end

end
% =========================================================================
function [imgType] = getimagetype( Img ) 
%GETIMAGETYPE   Returns image type as string 
% 
% imgType = GETIMAGETYPE( Img )
%
% Returns either 'PHASE', 'MAGNITUDE', or 'UNKNOWN'

if strfind( Img.Hdr.ImageType, '\P\' )
    imgType = 'PHASE' ;
elseif strfind( Img.Hdr.ImageType, '\M\' )
    imgType = 'MAGNITUDE' ;
else
    imgType = 'UNKNOWN' ;
end

end
% =========================================================================
function [isMag] = ismagnitude( Img )
%ISMAGNITUDE    Returns TRUE if Img is a magnitude image, FALSE otherwise.
% 
% isMag = ISMAGNITUDE( Img )

switch Img.getimagetype()
    case 'MAGNITUDE'
        isMag = true ;
    otherwise
        isMag = false ;
end

end
% =========================================================================
function [isPhase] = isphase( Img )
%ISPHASE    Returns TRUE if Img is a phase image, FALSE otherwise.
%
% isPhase = ISPHASE( Img )

switch Img.getimagetype()
    case 'PHASE'
        isPhase = true ;
    otherwise
        isPhase = false ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=private)
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

   
end
% =========================================================================
% =========================================================================


end
