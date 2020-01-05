function [ Hdr ] = dicominfosiemens( filename )
%DICOMINFOSIEMENS  Return standard DICOM Hdr struct with added fields from PARSE_SIEMENS_SHADOW(): .Img, .Ser, .MrProt 
%
% Hdr = DICOMINFOSIEMENS( filename ) 
%
% See also DICOMINFO, PARSE_SIEMENS_SHADOW

    Hdr = dicominfo( filename ) ;

    % parses additional info (e.g. "MrProt" (re: protocol))
    [ Hdr.Img, Hdr.Ser, Hdr.MrProt ] = parse_siemens_shadow( Hdr ) ;
    fprintf('\n') ;

    % parse_mrprot produces many, many warning messages (at least for DICOMs from the Prisma).
    % This should suppress all but the 1st instance:
    [ lastMsg, lastWarnId ] = lastwarn( ) ;
    warning( 'off', lastWarnId ) ;

end

