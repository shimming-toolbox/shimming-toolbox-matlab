function [ Hdr ] = dicominfosiemens( varargin )
%DICOMINFOSIEMENS  Return DICOM Hdr, augmented with fields from Siemens private header
%     
%      Hdr = dicominfosiemens( filename ) 
%      Hdr = dicominfosiemens( Hdr ) 
% 
% Wraps to `dicominfo` to return the standard header struct, and augments it
% with fields .Img, .Ser, .MrProt courtesy of `parse_siemens_shadow`.
%
% The input can be either the DICOM filename, or the corresponding Hdr struct as
% returned by `dicominfo( filename )`
% 
% __ETC__
% See also 
% dicominfo, parse_siemens_shadow
   
    if ( nargin ~= 1 ) && isempty( varargin{1} )
        error('Invalid input. See HELP dicominfosiemens') ;
    end
    
    switch class( varargin{1} )
        case {'char', 'string'}
            Hdr = dicominfo( varargin{1} ) ;
        case 'struct'
            Hdr = varargin{1} ;
        otherwise
            error('Invalid input. See HELP dicominfosiemens') ;
    end

    if strcmp( Hdr.Manufacturer, 'SIEMENS' )
        [ Hdr.Img, Hdr.Ser, Hdr.MrProt ] = parse_siemens_shadow( Hdr ) ;

        % parse_mrprot produces many warnings (at least for DICOMs from the Prisma).
        % This should suppress all but the 1st instance:
        [ lastMsg, lastWarnId ] = lastwarn( ) ;
        warning( 'off', lastWarnId ) ;
    else
        warning( ['Expected DICOM header field .Manufacturer == "SIEMENS".' ...
           ' Aborting call to parse_siemens_shadow()'] ) ;
        return ;
    end

end

