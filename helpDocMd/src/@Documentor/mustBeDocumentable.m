function [] = mustBeDocumentable( mFile )
%MUSTBEDOCUMENTABLE Issues an error if input path is not a valid script, function, or classdef .m file. 
%         
%     [] = mustBeDocumentable( mFile ) 
% 
% Throws an error if `mFile` does not meet the following requirements: 
%
% - `mFile` is (or, is convertible to) a string scalar
% - `mFile` is a valid file path, i.e. `isfile( mFile )` evaluates to `true`
% - `mFile` points to a valid script, function, or classdef file, i.e.
% `ismember( Informer.mfiletype( mFile ), ["script", "function","classdef"] ) == 1`
    arguments
        mFile(1,1) string {mustBeFile} ;
    end

    if ~ismember( Informer.mfiletype( mFile ), Documentor.mTypesSupported ) ;
        error('Input must be a path string to a valid script, function, or classdef .m file')
    end

end
