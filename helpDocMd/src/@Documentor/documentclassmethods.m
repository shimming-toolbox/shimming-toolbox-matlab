function docStr = documentclassmethods( Info )
%DOCUMENTCLASSMETHODS Return string vector of class method documentation
    arguments
        Info struct ;
    end
% Info = Dr.Info.Attributes ;

% assert( strcmp(Info.mType, "classdef"), 'mFile is not a class' ) ;

docStr = [ "---" ; "## Methods" ; "" ] ;

if isempty( Info.MethodList )
    docStr = [ docStr ; "" ; "[No Methods]"; "" ] ;
    return ;
else
    for iMthd = 1 : numel( Info.MethodList ) 
        docStr = [ docStr ; "" ; "---" ; "" ;] ;

        Mthd = Info.MethodList( iMthd ) ; 
         
        if contains( which( Mthd.DefiningClass ), 'built-in' ) 
        % MATLAB built-in method: include name and defining class only:
            docStr(end+1) = strcat( "### ", Mthd.Name ) ;
            docStr(end+1) = [ "[ _built-in method derived from class_ **" + Mthd.DefiningClass + "** ]" ];
            docStr(end+1) = strcat( "For more info, see MATLAB documentation]" ) ;
        else
            
            docStr = [ docStr ; "" ; strcat( "### ", Mthd.Name ) ; "" ] ; 
            
            if ~isempty( Mthd.Description )
                docStr = [ docStr ; strcat( "**Synopsis**: _", Mthd.Description, "_ " ) ; "" ] ;
            end

            if ~isempty( Mthd.DetailedDescription )
                % NOTE: could change Informer so it doesn't fill DetailedDescription with the copy of Description...
                if ~strcmp( Mthd.DetailedDescription, Mthd.Description )                         
                    docStr = [ docStr ; Mthd.DetailedDescription ; "" ] ;
                end
            end
    
            Mthd = rmfield( Mthd, {'Name'; 'Description' ; 'DetailedDescription'} ) ;
            mthdFields = fieldnames( Mthd ) ;

            %% Place basic (logical) attributes into an HTML table
            basicFields = mthdFields( structfun( @islogical, Mthd ) ) ;
            
            for iF = 1 : numel( basicFields )
               BasicAttributes.( basicFields{iF} ) = Mthd.( basicFields{iF} ) ;
            end
            
            docStr = [ docStr ; "" ; "#### Attributes:" ; "" ; Documentor.tableattributes( BasicAttributes ) ; "" ] ;
            
            clear BasicAttributes ;
            Mthd = rmfield( Mthd, basicFields ) ;
            
            %% 
            mthdFields = fieldnames( Mthd ) ;
     
            for iField = 1 : length( mthdFields ) 
            
                field = string( mthdFields( iField ) ) ;

                if isempty( Mthd.(field) )
                    docStr(end+1) = strcat( "- ", field, " : [N/A] " ) ;
                else 
                    docStr(end+1) = strcat( "- ", field, " : ", strjoin( string( Mthd.( field ) ), ", " ) ) ;
                end
            end
        end
    end
end

end
