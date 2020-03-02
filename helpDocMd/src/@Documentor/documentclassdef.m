function [docStr] = documentclassdef( Dr )
%DOCUMENTCLASSDEF Return string vector of class documentation 
%
% DOCUMENTCLASSDEF documents basic class attributes followed by class member
% documentation (courtesy of calls to Documentor.documentclassproperties and
% Documentator.documementclassmethods)

Info = Dr.Info.Attributes ;

assert( strcmp(Info.mType, "classdef"), 'mFile is not a class' ) ;

docStr = documentbasic( Dr ) ;

fields = string( fieldnames( Info ) ) ;

% remove fields included in documentbasic 
fields( fields=="mType" )               = [] ;
fields( fields=="Name" )                = [] ;
fields( fields=="Description" )         = [] ;
fields( fields=="DetailedDescription" ) = [] ;

docStr = [docStr ; "" ; "### Attributes ###"] ;

for iField = 1 : numel(fields)
    
    field = char( fields( iField ) ) ;
    
    if isempty( Info.(field) )
        docStr(end+1) = strcat( "- ", field, " : [N/A] " ) ;
    elseif strcmp( field, 'SuperclassList' )
        docStr(end+1) = strcat( "- Superclasses: ", strjoin( Info.SuperclassList, ", " ) ) ;
    elseif ~isstruct( Info.(field) ) % Property + MethodList etc. will be structs (handle them separately)
        docStr(end+1) = strcat( "- ", fields(iField), " : ", string( Info.( field ) ) ) ;
    end
end

docStr = [ docStr ; Dr.documentclassproperties ] ;    
docStr = [ docStr ; Dr.documentclassmethods ] ;    

end
    % names = strcat( "_", join(names, ", "), "_" ) ; % italicize
    % info  = [ info ; strcat( "- Parent classes: ", names ) ] ;
    %
    % if isempty( Mc.InferiorClasses ) 
    %     names = "(None)" ;
    % else 
    %     % NOTE: 'Mc.InferiorClasses' is a cell array whereas Mc.SuperclassList
    %     % is an object array, hence the for-loop: 
    %     names = string( Mc.InferiorClasses{1}.Name ) ;
    %     for iClass = 2 : numel(Mc.InferiorClasses)
    %         names = [ names string( Mc.InferiorClasses{ iClass }.Name ) ] ;
    %     end
    % end
    %
    % names = strcat( "_", join(names, ", "), "_" ) ; % italicize
    % info  = [ info ; strcat( "- Child classes: ", names ) ] ;
    %
    % if isempty( Mc.ContainingPackage ) 
    %     pkg = "N/A" ;
    % else
    %     pkg = string( Mc.ContainingPackage.Name ) ;
    % end
    %
    % info = [ info ; strcat( "- Containing Package: ", "_", pkg, "_" ) ] ;

% function [info] = documentclassattributes( Mc )
% %DOCUMENTCLASSATTRIBUTES        
%     arguments
%         Mc(1,1) meta.class ;
%     end
%
%     info = [ "" ; "#### Attributes ####" ; "" ] ; 
%     
%     fields = fieldnames( Mc ) ;
%
%     %% Basic attributes:
%     for iField = 1 : length( fields )
%         if islogical( Mc.(fields{iField}) ) 
%             % Applies to: "Hidden", "Sealed", "Abstract", "Enumeration",
%             % "ConstructOnLoad", "HandleCompatible", "RestrictsSubclassing"
%             value = join( string( Mc.( fields{iField} ) ), ", " ) ;
%             info  = [ info ; strcat( "- ", fields{iField}, ": ", value ) ] ;
%         end
%     end
%
%     %% Inheritances and package info:
%     %
%     %TODO when possible: add links to packages, parents, subclasses
%
%     if isempty( Mc.SuperclassList ) 
%         names = "(None)" ;
%     else
%         names = string( {Mc.SuperclassList.Name} ) ;
%     end
%
%     names = strcat( "_", join(names, ", "), "_" ) ; % italicize
%     info  = [ info ; strcat( "- Parent classes: ", names ) ] ;
%
%     if isempty( Mc.InferiorClasses ) 
%         names = "(None)" ;
%     else 
%         % NOTE: 'Mc.InferiorClasses' is a cell array whereas Mc.SuperclassList
%         % is an object array, hence the for-loop: 
%         names = string( Mc.InferiorClasses{1}.Name ) ;
%         for iClass = 2 : numel(Mc.InferiorClasses)
%             names = [ names string( Mc.InferiorClasses{ iClass }.Name ) ] ;
%         end
%     end
%     
%     names = strcat( "_", join(names, ", "), "_" ) ; % italicize
%     info  = [ info ; strcat( "- Child classes: ", names ) ] ;
%
%     if isempty( Mc.ContainingPackage ) 
%         pkg = "N/A" ;
%     else
%         pkg = string( Mc.ContainingPackage.Name ) ;
%     end
%
%     info = [ info ; strcat( "- Containing Package: ", "_", pkg, "_" ) ] ;
%
% end
