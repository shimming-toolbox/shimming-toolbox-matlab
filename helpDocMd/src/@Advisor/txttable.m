function [txt] = txttable( isCopying, nSize, nCharPerCell, prefix ) 
%TXTTABLE Create and [return|copy|print] a character table
%
%   [~] = txttable( 1, nSize, nSizeCell, "%" ) 
% 
% Copies an empty character table with `nSize(1)` rows and `nSize(2)` columns
% to the clipboard and prints to the standard output. nSizeCell(2) sets the
% number of ' ' spaces each cell is allotted in the horizontal dimension.
% (Note: nSizeCell(2) is not yet used/implented.)
% 
% 1. Example 
%
% >> Helper.txttable(1,[3 5],[1 10])% outputs:
% '
% |          |          |          |          |          |
% |----------|----------|----------|----------|----------|
% |          |          |          |          |          |
% |          |          |          |          |          |
%
%
% When called without any input or output arguments, TABLE prints the
%  to the standard output. 
% When the first input is `[ 1 | true ]`,`txt` is copied to the clipboard.
%
% __DESCRIPTION__
%
% To suppress the print display, call the function with a return argument, 
% e.g. using the dummy argument `[~] = txttable();`.
    arguments
        isCopying(1,1) {mustBeBoolean}                     = false ;
        nSize(1,2) {mustBePositive, mustBeInteger}         = [3 2] ;
        nCharPerCell(1,2) {mustBePositive, mustBeInteger}  = [1 10] ;
        prefix {mustBeStringScalarOrCharVector}            = "%    " ; 
        % fName {mustBeStringScalarOrCharVector}             = tempname ;
    end

%------------------------------ 
%% Create table

cell1 = [ '|' repmat( ' ', [1 nCharPerCell(2)] ) ] ;
cell1 = repmat( cell1, [ nCharPerCell(1) 1 ] ) ;

row1  = [ repmat( cell1, [1 nSize(2) ] ) '|'] ;
row2  = replace( row1, ' ', '-' ) ;

cTable = [ row1 ; row2 ; repmat(row1,[nSize(1)-1 1] ) ] ;

% sandwich between empty lines
cTable = [ prefix ; repmat( prefix,[size(cTable,1) 1 ] ) + string(cTable) ; prefix ] ;

cTable = compose( cTable + '\n' ) ;

%------------------------------ 
%% Output
txt = sprintf( '%s', cTable ) ;

if isCopying
    clipboard( 'copy', txt ) ;
end

if nargout == 0
    fprintf(txt) ;
    clear txt ; 
    return
end

end

% function [DEFAULTS]   = getdefaults( )
%     
%     DEFAULTS.nSize        = [3 3] ;
%     DEFAULTS.nSizeCell    = [1 10] ;
%     DEFAULTS.isCommented  = [true] ; % precede 1st column with commentChar
%     DEFAULTS.commentChar  = '%' ; % character to precede each row ;
%     DEFAULTS.nCharIndent  = 1 ;
%     DEFAULTS.isCopying    = [true] ;
%
% end
