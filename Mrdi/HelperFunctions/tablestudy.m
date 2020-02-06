function [ studyDirs ] = tablestudy( sortedDicomDir )
%TABLESTUDY Return a cell array of MaRdI-loadable images series within the given directory
%
%   [ studyDirs ] = MaRdI.tablestudy( sortedDicomDir ) 
%
% e.g. Protocol to load MaRdI-object :
%
%   % omit semi-colon to display the full cell array (i.e. with the row indices)
%   [ studyDirs ] = MaRdI.tablestudy( sortedDicomDir ) 
%
%   % determine the row index of the series you want to load (e.g. 10):
%   Img = MaRdI( studyDirs{ 10, 2 } ) ;
% 
% (The 1st column of studyDirs is merely the row index.)
%
% NOTE: Possibly deprecated.

assert( nargin == 1, 'Function requires sortedDicomDirectory as input argument.' ) ;

if ~strcmp( sortedDicomDir(end), '/' ) 
    sortedDicomDir(end+1) = '/' ;
end

studyDirs = cell( 0 ) ;

Tmp      = dir( [ sortedDicomDir ] );
Tmp      = Tmp( 3:end ) ; % ignore self ('.') and parent ('..') dirs
nEntries = length( Tmp ) ;

for iEntry = 1 : nEntries 

   if Tmp( iEntry ).isdir
   
       tmpSeriesSubdir = [ Tmp( iEntry ).name '/'] ;
    
        TmpEchoSubdirs = dir( [ sortedDicomDir tmpSeriesSubdir 'echo*' ] ) ;
        nEchoSubdirs   = length( TmpEchoSubdirs ) ;

        if nEchoSubdirs ~= 0

            for iEchoSubdir = 1 : nEchoSubdirs

                studyDirs{end+1, 2} = strcat( sortedDicomDir, tmpSeriesSubdir, TmpEchoSubdirs(iEchoSubdir).name )  ;
                iSeries = size( studyDirs, 1 ) ;
                studyDirs{ iSeries, 1 } = iSeries ;
            end

        % check if tmpSeriesSubdir itself contains images
        elseif length( dir( [ sortedDicomDir tmpSeriesSubdir '/*.dcm'] ) ) ~= 0 || ...
                length( dir( [ sortedDicomDir tmpSeriesSubdir '/*.IMA'] ) ) ~= 0 

           studyDirs{end+1, 2} = strcat( sortedDicomDir,tmpSeriesSubdir )  ;
            iSeries = size( studyDirs, 1 ) ;
            studyDirs{ iSeries, 1 } = iSeries ;

        end

   end
   
end

end

