function [Imgs] = make( varargin )
%MAKE
% 
%   Function to make generic Mrdi object and convert (when possible) to known child-class.    
% TODO
    [Imgs, Hdrs] = MrdiIo.loadandsortimages( varargin{:} )

end

