function [] = draftdoccontent( Dr, detailLevel )  
%DRAFTDOCCONTENT Generates and assigns formatted documentation to `docContent`
%     
%     Dr.draftdoccontent( )
%     Dr.draftdoccontent( detailLevel )
%
% Calls `Informer` to retrieve info on the file assigned to `Dr.mFile` and
% assigns a formatted string vector to `Dr.docContent` which can then be
% printed to file via `Dr.printdoc()`.

% `draftdoc` uses the with whatever additional information can
% be gathered about the file assigned to `Dr.mFile` and 

    % Configures the degree of detail to include in the documentation
    % 0: The documentation is the same as the header comment
    % 1: Add additional info when possible (e.g. class property attributes, hidden members, etc.)
    % NOTE: only partially implemented! 
% `isDetailed` is a Boolean toggle
%
% Toggles between basic/user (=false) and detailed/developer documentation (=true) [default: true]
% 
% When false, classes and class members with private, protected, or hidden
% attributes are excluded from the output documentation. [default = true]
%
% TODO: Choose desired formating for documentation methods (e.g. documentbasic, documentfunction, etc)!
    arguments
        Dr Documentor
        detailLevel {mustBeMember(detailLevel,[0 1])} = ones(size(Dr)) ;
    end
    
    if numel(Dr)>1 && numel(detailLevel)==1
        detailLevel = repmat( detailLevel, [numel(Dr) 1] ) ; 
    elseif numel(detailLevel)~=numel(Dr) 
        error( ['Input assignment for `detailLevel` must be a scalar ' ...
                'or an array with the same number of elements as the given `Documentor` object' ] )
    end

    fprintf( strcat("Generating doc content (x/", num2str(numel(Dr)), "): " ) ) ;
    
    for iM = 1 : numel(Dr)
        fprintf( [ num2str(iM) '...' ] ) ;

        Info = Informer( Dr(iM).mFile ) ;
        Info = Info.Attributes ;

        switch Info.mType{:}
            case 'script'
                Dr(iM).docContent = Documentor.documentbasic( Info ) ; 
            case 'function'
                Dr(iM).docContent = Documentor.documentfunction( Info ) ;
            case 'classdef'
                Dr(iM).docContent = Documentor.documentclassdef( Info, detailLevel(iM) ) ;
        end
    end
        
    fprintf( '\n' ) ;

end

