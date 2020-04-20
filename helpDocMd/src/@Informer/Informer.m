classdef Informer 
%INFORMER Details the functionality of a .m file
% 
% The `Informer` class informs the `Documentor` class about the contents of .m
% files: providing the functional details needed to document them. 
%
% In general, the class works behind the scenes and would not be directly
% called upon by a user; however, an Informer object can be constructed
% independently as indicated below.
%
% __CONSTRUCTOR SYNTAX__
%     
%     Info = Informer( mFile ) ;
% 
% Creates an `Informer` object pertaining to .m file (a script, function, class
% method, or classdef file) pointed to by the path string `mFile`.
% If `mFile` contains multiple files, then `Info` is returned as an
% object-array.
% 
% `Informer` has but two properties: `mFile` and `Attributes`.
%
% `Attributes` is a struct containing all available functional details
% regarding `mFile`. 
%
% `Attributes` cannot be set directly, but is updated whenever `mFile` is set.
%
% The read-only fields of `Attributes` depend on the given .m-file type and
% should be fairly self-explanatory given the field names. More detail is
% available in the method documentation for Informer.getmattributes.
%
 

properties( AbortSet = true )

    % .m file path: Attributes property will update whenever mFile is set
    mFile {mustBeFile} = string( [ mfilename('fullpath') '.m'] ) ;

end

properties( SetAccess=private )
    
    % Functional description of the .m file
    %
    % Attributes is a struct of .m file attributes, returned from a call to
    % Informer.getmattributes( mFile ). It contains the following basic
    % fields:
    % 
    % - mType: Type of .m file: string scalar returned from Informer.mfiletype( mFile ).
    % Possibilities are: ["script","function","classdef","method","NA"]
    %   
    % - .Name: Name of the script, function, class or class method 
    % 
    % If the .m file is a valid MATLAB file (i.e. mType ~= "NA"), then
    % Attributes also contains fields:
    %
    % - .Description: Header line of help-text (string vector returned from
    %   Informer.extracthelpheader) 
    %
    % - .DetailedDescription: Body of help-text (string vector returned from
    %   Informer.extracthelpbody) 
    % 
    % ### References ###
    %
    % Remaining fields of Attributes vary depending on the type of file 
    % (i.e. Attributes.mType). For more info, see the method documentation:
    % Informer.getmattributes
    Attributes struct ;
    
end

 
% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function [ Info ] = Informer( mFile )
    
    if nargin == 0
        return ;
    else
        mustBeStringOrChar( mFile ) ;
    end

    if ~isfile( mFile )
        mFile = which( mFile ) ;
        if isempty( mFile )
            error('File not found') ;
        end
    end
   
    mFile = string( mFile ) ;    
    
    Info.mFile      = mFile(1) ;
    Info.Attributes = Informer.getmattributes(mFile(1)) ;

    if numel( mFile ) > 1
        
        fprintf( strcat("Retrieving file info (x/", num2str(numel(mFile)), "): 1..." ) ) ;
        
        for iM = 2 : numel( mFile )
            fprintf( [ num2str(iM) '...' ] ) ;
            Info(iM).mFile      = mFile(iM) ;
            Info(iM).Attributes = Informer.getmattributes(mFile(iM)) ;
        end
        
        fprintf( '\n' ) ;
    
        Info = Info' ;
    end

end
% =========================================================================    
function [Info] = set.mFile( Info, mFile )
   
   Info.mFile      = mFile ; 
   Info.Attributes = Informer.getmattributes( mFile ) ;

end
% =========================================================================    

% =========================================================================    
% =========================================================================    
end

% =========================================================================    
% =========================================================================    
methods( Static )
    %.....
    [mHelpBody] = extracthelpbody( mHelp )
    %.....
    [mHelpHeader] = extracthelpheader( mHelp, name )
    %.....
    [mHelp] = gethelptext( name )
    %.....
    [Att] = getmattributes( mFile )
    %.....
    [mType, mPath, mExist] = mfiletype( mFile )
end
% =========================================================================    
% =========================================================================    
methods( Static, Hidden )
    %.....
    [Att] = metainfo( Mc )
end
% =========================================================================    
% =========================================================================    

end
