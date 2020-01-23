function [ Hdr ] = dicominfosiemens( varargin )
%DICOMINFOSIEMENS  Return standard DICOM Hdr struct with added fields from PARSE_SIEMENS_SHADOW(): .Img, .Ser, .MrProt 
%
% Hdr = DICOMINFOSIEMENS( filename ) 
% Hdr = DICOMINFOSIEMENS( Hdr ) 
% 
% Input can be either the DICOM filename, or the corresponding Hdr struct as
% returned by DICOMINFO( filename )
%
% See also DICOMINFO, PARSE_SIEMENS_SHADOW
   
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
        fprintf('\n') ;

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

