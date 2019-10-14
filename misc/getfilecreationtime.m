function [ creationTime ] = getfilecreationtime( filename ) 
%GETFILECREATIONTIME   Query system to return file creation time-stamp 
%
% creationTime = GETFILECREATIONTIME( filename ) 
% 
% Returned creationTime is a double formated as : yyyymmddHHMMSS.FFF
%
% NOTE//TODO: Implementation only supports MacOS// Add support for other operating systems!
%
% Based on forum comment from Tim Leuth:
% https://www.mathworks.com/matlabcentral/answers/288339-how-to-get-creation-date-of-files

if ismac
    [~, msg] = system( ['GetFileInfo ' filename ] ) ; 
    i = strfind( msg, 'created: ') + 9 ; 
    creationTime = msg(i:i+18) ;
    creationTime =  str2num( datestr( creationTime, 'yyyymmddHHMMSS.FFF') ) ;
else
    error('Unimplemented system. MacOS only! Todo.') ;
end

end
