classdef Informer
%INFORMER Details the functionality of a .m file
% 
% An INFORMER object describes a given .m file (i.e. a script, function, class
% method, or classdef file).
%
% The INFORMER class _informs_ the *Documentor* class about the content of .m
% files: providing the funtional details needed to document them.
%
% ### Syntax
%    
%    Info = Informer( mFile ) ;
% 
% ### Inputs
% 
% - mFile: full file path to the .m file of interest.
%
% ### Usage
% 
% In the general/anticipated use case, a user merely interacts with an Informer
% instance indirectly, as a member property of a Documentor object; however,
% an Informer object can be constructed independently as indicated above.
%
% Informer only has two properties public properties: mFile and Attributes,
% a struct storing all available information on the .m file. Info.Attributes
% cannot be set directly, but is updated each time Info.mFile is set. The
% read-only fields of Informer.Attributes depend on the given .m-file type and
% should be fairly self-explanatory given the field names; nevertheless, more
% detail is available in the method documentation for Informer.getmattributes.

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

   Info.mFile      = string( mFile ) ;    
   Info.Attributes = Informer.getmattributes( mFile ) ;

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
    [mHelp] = gethelptext( name )
    %.....
    [Att] = getmattributes( mFile )
    %.....
    [mHelpBody] = extracthelpbody( mHelp )
    %.....
    [mHelpHeader] = extracthelpheader( mHelp, name )
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
