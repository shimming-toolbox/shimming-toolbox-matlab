function [tableStr] = tableattributes( Attributes )  
%TABLEATTRIBUTES Return html-table of class/classmember attributes 
%
% tableStr = GETHELPTEXT( Attributes )

fields = string( fieldnames( Attributes ) ) ;

titles = "<table border=1><tr>" ; % title row
values = "<tr>" ;

for iField = 1 : numel( fields )
    field = fields( iField ) ;
    titles = titles + "<th>"+ field + "</th>" ;
    values =  values+ "<td>"+ string( Attributes.( field ) ) + "</td>" ;
end

titles = titles + "</tr>" ;
values = values+ "</tr>"  ;

tableStr = [ "<table>" ; titles; values ; "</table>" ] ;

end
