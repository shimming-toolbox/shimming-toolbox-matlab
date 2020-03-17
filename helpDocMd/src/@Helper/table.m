function [txt] = table( isCopying, nSize, nCharPerCell ) 
%TABLE Create and [return|copy|print] a character table
%
%   [~] = table( 1, nSize, nSizeCell ) 
% 
% Copies an empty character table with `nSize(1)` rows and `nSize(2)` columns
% to the clipboard and prints to the standard output. nSizeCell(2) sets the
% number of ' ' spaces each cell is allotted in the horizontal dimension.
% (Note: nSizeCell(2) is not yet used/implented.)
% 
% 1. Example 
%
% >> Helper.table(1,[3 5],[1 10])% outputs:
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
% e.g. using the dummy argument `[~] = txt();`.
    arguments
        isCopying(1,1) {mustBeBoolean}                     = false ;
        nSize(1,2) {mustBePositive, mustBeInteger}         = [3 3] ;
        nCharPerCell(1,2) {mustBePositive, mustBeInteger}  = [1 10] ;
        % fName {mustBeStringScalarOrCharVector}             = tempname ;
    end

%TODO implement nCharsPerCell(1) >1 spacing (not comptable w/Markdown, but could still look good

%------------------------------ 
%% Create table
cell1 = [ '|' repmat( ' ', [1 nCharPerCell(2)] ) ] ;
row1  = [ repmat( cell1, [1 nSize(2) ] ) '|'] ;
row2  = replace( row1, ' ', '-' ) ;

cTable = [ row1 ; row2 ; repmat(row1,[nSize(1)-1 1] ) ] ;
cTable = compose( string(cTable) + '\n' ) ;


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
