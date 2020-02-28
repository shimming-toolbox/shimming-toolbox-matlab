classdef Informer
%INFORMER Details the functionality of a .m file
% 
% An INFORMER object instance describes a given .m file (i.e. a script,
% function, class method, or classdef file).
%
% The INFORMER class serves to _inform_ the *Documentor* class regarding the
% content of .m files: providing the funtional details needed to
% document them.
%
% ### Syntax ###
%
% Info = Informer( mPath ) ;
% 
% ### Inputs ###
% 
% - mPath: full file path to the .m file of interest.
%
% ### Usage ###
% 
% Informer only has two properties public properties: mPath and Attributes,
% a struct storing all available information on the .m file. Info.Attributes
% cannot be set directly, but is updated each time Info.mPath is set.

properties( AbortSet = true )

    % .m file path : Attributes property will update whenever mPath is set
    mPath {mustBeFile} = string( [ mfilename('fullpath') '.m'] ) ;
end

properties( SetAccess=private )
    
    % Functional description of the .m file
    %
    % Attributes is a struct of .m file attributes, returned from a call to
    % Informer.getmattributes( mPath ). It contains the following basic
    % fields:
    % 
    % - mType: Type of .m file: string scalar returned from Informer.mfiletype( mPath ).
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
function [ Info ] = Informer( mPath )
    
    if nargin == 0
        return ;
    else
        mustBeStringOrChar( mPath ) ;
    end

    if ~isfile( mPath )
        mPath = which( mPath ) ;
        if isempty( mPath )
            error('File not found') ;
        end
    end

   Info.mPath      = string( mPath ) ;    
   Info.Attributes = Informer.getmattributes( mPath ) ;

end
% =========================================================================    
function [Info] = set.mPath( Info, mPath )
   
   Info.mPath = mPath ; 
        
   Info.Attributes = Informer.getmattributes( mPath ) ;

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
    [mCode] = getcode( name )
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
