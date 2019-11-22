function [ creationTime ] = getfilecreationtime( filename ) 
%GETFILECREATIONTIME   Query system to return file creation time-stamp 
%
% creationTime = GETFILECREATIONTIME( filename ) 
% 
% Returned creationTime is a string formated as : yyyymmddHHMMSS
%
% NOTE: Implementation only supports UNIX/LINUX
%
% TODO: Add support for other operating systems!

if isunix
    [~,creationTime] = system( ['date -r ' filename ' +"%Y%m%d%H%M%S"'] ) ; 
else
    error('Unimplemented system. Unix only! Todo.') ;
end

end
