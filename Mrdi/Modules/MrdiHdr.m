classdef MrdiHdr < dynamicprops
%MrdiHdr  MR DICOM Image Header

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% properties( SetAccess=protected, Hidden=true ) 
%     
%     % NOTE: 
%     % Best attributes for the .Hdrs property are open to debate:  
%     % Set access should probably be kept protected so nothing dubious happens
%     % to the metadata. However, re: GetAccess, though it is currently 'Hidden'
%     % (readable from anywhere --> simplifies debugging), ultimately, it should
%     % probably either be made visible to the user (if useful), or
%     % concealed entirely (simpler interface) with GetAccess=private.
%
%     % Struct array of metadata (e.g. DICOM headers)
%     Hdrs struct = struct( [] ) ;
%     
% end


% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Hdr = MrdiHdr( varargin )

    if nargin == 0
        return ;
    
    elseif isstruct( varargin{1} )
        
        names = fieldnames( varargin{1} ) ;

        for iProp = 1: length( names )
            Hdr.addprop( names{iProp} ) ;
            Hdr.( names{iProp} ) = varargin{1}.( names{iProp} ) ;
        end
    
    else
        error('Invalid input') ;
    end

end
% =========================================================================



end
% =========================================================================
% =========================================================================    


end
