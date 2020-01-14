classdef (Abstract) MrdiUtil < handle
%MrdiUtil   MR DICOM Image Utility  
%
% Member methods for handling Mrdi objects 
% 
% e.g. public methods:
%
%   copy()
%   getproperties()
%
% e.g. protected methods:
%
%   copyproperties()
%   exist()
%
%   ...etc.
% 
% =========================================================================

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties
    img = [] ; 
    Aux = [] ;
% end
%
% properties(SetAccess={?MaRdI, ?MrdiIo, ?MrdiUtil, ?MrdiProc})
    Hdr = []; % full Siemens DICOM header of 1st img (i.e. Img.img(:,:,1) )
    Hdrs = []; % cell array of (truncated) DICOM headers courtesy of dicominfo()
    Ref =[]; % Reference properties - prior to manipulation
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = MrdiUtil( imgPath )

    % if ~isstr( imgPath )
    %     eval( ['help ' mfilename ] ) ;
    %     error('no no no') ;
    % end

end
% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================    
methods 
% =========================================================================
function ImgCopy = copy( Img )
%COPY   Return copy of a Mrdi-object.
% 
% ImgCopy = COPY( Img ) ;

constructorHandle = str2func( class( Img ) ) ;
ImgCopy = constructorHandle() ;

ImgCopy.copyproperties( Img ) ;

end
% =========================================================================
function Props = getproperties( Img )
%GETPROPERTIES  Return struct containing copies of a Mrdi-object's properties
% 
% Props = GETPROPERTIES( Img ) ;

C   = metaclass( Img ) ;
Tmp = C.Properties ;

for iProp = 1: length( Tmp )
    if ~Tmp{iProp}.Dependent
       Props.( Tmp{iProp}.Name ) = Img.( Tmp{iProp}.Name ) ;
    end
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=protected)
% =========================================================================
function [] = copyproperties( ImgCopy, ImgOriginal )
%COPYPROPERTIES  Copies properties from one MaRdI-object to another
%
% COPYPROPERTIES( Copy, Original ) 
%
% Copies the properties of Mrdi-object Original to Mrdi-object Copy 
% (overwriting its initial values)

assert( isa( ImgOriginal, 'MrdiUtil' ) ) ;
Props = ImgOriginal.getproperties() ;
names = fieldnames( Props ) ;

for iProp = 1: length( names )
    ImgCopy.( names{iProp} ) = ImgOriginal.( names{iProp} ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=protected, Hidden=true)
% =========================================================================
function [is] = exist( Img ) 
%EXIST  Return true (object exists)

is = true ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================


end
% =========================================================================
% =========================================================================

end
