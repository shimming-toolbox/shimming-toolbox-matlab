function [] = SortData( unsortedDicomDir, sortedDicomDir)
% SortData : 
%
% Sort raw dicom files from a scanner in a folder called sorted_data
% 
% [] = SORTDATA( unsortedDicomDirectory ) 
% [] = SORTDATA( unsortedDicomDirectory, sortedDicomDirectory ) 
% 
% If sortedDicomDirectory is unspecified, a subdirectory ('sorted') is created
% within unsortedDicomDirectory to contain the sorted images.
%  
%--------------------------------------------------------------------------
DEFAULT_SORTEDDICOMDIR = [ unsortedDicomDir '/sorted/' ] ;

%Dicom files inside the directory -----------------------------------------
listOfImages = dir( [ unsortedDicomDir '/*.dcm'] );
nImages      = length(listOfImages);

if length(listOfImages) == 0
    % try .IMA
    listOfImages = dir( [imgDirectory '/*.IMA'] ) ;
end

assert( length(listOfImages) ~= 0, 'No .dcm or .IMA files found in given directory' ) ;

%Declare initial parameters -----------------------------------------------

if nargin < 2 || isempty( sortedDicomDir )
    sortedDicomDir = DEFAULTSORTEDDICOMDIR ;
end

if ~exist( sortedDicomDir,'dir') 
    mkdir( sortedDicomDir ); 
end

seriesDescription = [];

nImagesnew=0;

%Read Dicom files header and extract series names and numbers -------------

Hdr = dicominfo( [unsortedDicomDir '/' listOfImages(1).name]);
    
    
folderName = [ num2str(Hdr.SeriesNumber) '_' Hdr.SeriesDescription ];

if ( num2str(Hdr.SeriesNumber) < 10)
    folderName = [ '0' folderName ] ;
end

iSeries=1;
 
seriesDescription{iSeries}=folderName;

sortedFilename = fullfile( sortedDicomDir, folderName );

% Create directories for each series
if ~exist( sortedFilename,'dir' )
  mkdir(sortedDicomDir,folderName);
end

% Copy the images into the new directory
copyfile(Hdr.Filename, sortedFilename);
    

while(1)   
    
    if( nImagesnew~=nImages)
       display('Sorting... (in progress)');   %Feedback
   
   
    %Check the presence for new files in the raw data directory----------------
    newFiles=abs( nImagesnew - nImages);
       
    listOfImages=dir( [ unsortedDicomDir '/*.dcm'] );
            if(nImagesnew==0)
                a=2;              %1st file is already sorted
                b=nImages;
            else
                a=nImages+1;
                b=nImages+newFiles;
            end
        nImages =length(listOfImages);
        
        for iscan=a:b
            newHdr = dicominfo( [unsortedDicomDir '/' listOfImages(iscan).name] ) ;
        
            if strcmp(newHdr.SeriesNumber,Hdr.SeriesNumber)==1
                  copyfile(newHdr.Filename,sortedFilename);
            else
            newseriesName=newHdr.SeriesDescription;
            newfolderNumber=newHdr.SeriesNumber;
            
                  if (newfolderNumber < 10)
                    folderName= strcat('0',num2str(newfolderNumber),'_',newseriesName); %Add 0 for consistency between all folders.
                  else
                     folderName= strcat(num2str(newfolderNumber),'_',newseriesName);
                  end
                
                  if ~any(strcmp(seriesDescription,folderName));
                        iSeries=iSeries+1;
                        seriesDescription{iSeries}= folderName;
                        sortedFilename = fullfile (sortedDicomDir, folderName);
                    
                   if ~exist(sortedFilename,'dir')                 
                            mkdir(sortedDicomDir,folderName);
                   end
                   end
                   
                   copyfile(newHdr.Filename,sortedFilename);        
                   Hdr = newHdr;
             end
     
        end





%Sort files inside series folder ------------------------------------------
    for i= 1:iSeries
        sortedFilenametoseries = fullfile (sortedDicomDir, seriesDescription{i});
%Dicom files inside series folder -----------------------------------------
        newDicomfilesinseries=dir( [sortedFilenametoseries '/*.dcm'] );
        nFilesinseries=length(newDicomfilesinseries);
    
        if ~isempty(newDicomfilesinseries) 
        
%Read Dicom files header and extract imageType and Echonumber -------------
        filesinseriesHeader = dicominfo( [sortedFilenametoseries '/' newDicomfilesinseries(1).name]);
        imageType=filesinseriesHeader.ImageType;
        name=fullfile(imageType,num2str(filesinseriesHeader.EchoNumber));
        sortedFilenametoseriestype = fullfile (sortedFilenametoseries,name);
        
%Create new directories for each imageType and Echo number-----------------
        if ~exist(sortedFilenametoseriestype,'dir')
            mkdir(sortedFilenametoseries,name);
        end
        
        movefile(filesinseriesHeader.Filename,sortedFilenametoseriestype);
    
            for j=2:nFilesinseries
                 newfilesinseriesHeader = dicominfo( [sortedFilenametoseries '/' newDicomfilesinseries(j).name]);
            
                if strcmp(newfilesinseriesHeader.ImageType,Hdr.ImageType)==1
                    movefile(newfilesinseriesHeader.Filename,sortedFilenametoseriestype);
                else
                    foldertype=newfilesinseriesHeader.ImageType;
                    name=fullfile(foldertype,num2str(filesinseriesHeader.EchoNumber));
                    sortedFilenametoseriestype = fullfile (sortedFilenametoseries, name);
                        if ~exist(sortedFilenametoseriestype,'dir')
                        mkdir(sortedFilenametoseries,name);
                        end
                    movefile(newfilesinseriesHeader.Filename,sortedFilenametoseriestype);
                end
                filesinseriesHeader = newfilesinseriesHeader;
            end
        end
       
     end

 listOfImagesnew=dir( [ unsortedDicomDir '/*.dcm'] );
 nImagesnew =length(listOfImagesnew);
 display('Files in raw data directory are sorted');  

else    

  listOfImagesnew=dir( [ unsortedDicomDir '/*.dcm'] );
  nImagesnew =length(listOfImagesnew);  

end

end
