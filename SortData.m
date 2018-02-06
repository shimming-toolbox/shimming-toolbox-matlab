%==========================================================================
%Function to sort the Raw data from the scan 
%==========================================================================

function SortData(rawdataDirectory,sortdataDirectory)

%Dicom file inside the directory ------------------------------------------
dicomFiles=dir( [ rawdataDirectory '/*.dcm'] );
nDicomfiles =length(dicomFiles);
%Declare initial parameters -----------------------------------------------
seriesDescription=[];
nDicomfilesnew=0;
%Writing a feedback in a textfile background------------------------------- 
diary('background');

if ~exist(sortdataDirectory,'dir'), mkdir(sortdataDirectory); end

%Read Dicom files header and extract series name and number ---------------
    filesHeader = dicominfo( [rawdataDirectory '/' dicomFiles(1).name]);
    
    seriesName=filesHeader.SeriesDescription;
    
    folderNumber=filesHeader.SeriesNumber;
    
    nseriesDescription=1;
    
    if (folderNumber < 10)
    folderName= strcat('0',num2str(folderNumber),'_',seriesName);
    else
    folderName= strcat(num2str(folderNumber),'_',seriesName);
    end
    seriesDescription{nseriesDescription}=folderName;
  
    newPath = fullfile (sortdataDirectory, folderName);

   if ~exist(newPath,'dir')
      %Create new directories for each series name-------------------------
      mkdir(sortdataDirectory,folderName);
   end
      %Copy the dicom files in the new directory---------------------------
    copyfile(filesHeader.Filename,newPath);
    

while(1)   
    
if(nDicomfilesnew~=nDicomfiles)
   display('Sorting in process');
   diary('background');
   
   
   %Check the presence for new files in the raw data directory-------------
   newFiles=abs(nDicomfilesnew-nDicomfiles);
   dicomFiles=dir( [ rawdataDirectory '/*.dcm'] );
        if(nDicomfilesnew==0)
            a=2;
            b=nDicomfiles;
        else
            a=nDicomfiles+1;
            b=nDicomfiles+newFiles;
        end
    nDicomfiles =length(dicomFiles);
    for iscan=a:b
        newfilesHeader = dicominfo( [rawdataDirectory '/' dicomFiles(iscan).name] ) ;
    
        if strcmp(newfilesHeader.SeriesNumber,filesHeader.SeriesNumber)==1
              copyfile(newfilesHeader.Filename,newPath);
        else
        newseriesName=newfilesHeader.SeriesDescription;
        newfolderNumber=newfilesHeader.SeriesNumber;
        
              if (newfolderNumber < 10)
                folderName= strcat('0',num2str(newfolderNumber),'_',newseriesName);
              else
                 folderName= strcat(num2str(newfolderNumber),'_',newseriesName);
              end
            
              if ~any(strcmp(seriesDescription,folderName));
                    nseriesDescription=nseriesDescription+1;
                    seriesDescription{nseriesDescription}= folderName;
                    newPath = fullfile (sortdataDirectory, folderName);
                
               if ~exist(newPath,'dir')                 
                        mkdir(sortdataDirectory,folderName);
               end
               end
               
               copyfile(newfilesHeader.Filename,newPath);
               
               
               filesHeader = newfilesHeader;
         end
 
    end

    for i= 1:nseriesDescription
        newPathtoseries = fullfile (sortdataDirectory, seriesDescription{i});
        %Dicom files inside the new directories----------------------------
        newDicomfilesinseries=dir( [newPathtoseries '/*.dcm'] );
        nFilesinseries=length(newDicomfilesinseries);
    
        if ~isempty(newDicomfilesinseries) 
        
        filesinseriesHeader = dicominfo( [newPathtoseries '/' newDicomfilesinseries(1).name]);
        imageType=filesinseriesHeader.ImageType;
        name=fullfile(imageType,num2str(filesinseriesHeader.EchoNumber));
        newPathtoseriestype = fullfile (newPathtoseries,name);
        if ~exist(newPathtoseriestype,'dir')
        mkdir(newPathtoseries,name);
        end
        movefile(filesinseriesHeader.Filename,newPathtoseriestype);
    
            for j=2:nFilesinseries
                 newfilesinseriesHeader = dicominfo( [newPathtoseries '/' newDicomfilesinseries(j).name]);
            
                if strcmp(newfilesinseriesHeader.ImageType,filesHeader.ImageType)==1
                    movefile(newfilesinseriesHeader.Filename,newPathtoseriestype);
                else
                    foldertype=newfilesinseriesHeader.ImageType;
                    name=fullfile(foldertype,num2str(filesinseriesHeader.EchoNumber));
                    newPathtoseriestype = fullfile (newPathtoseries, name);
                        if ~exist(newPathtoseriestype,'dir')
                        mkdir(newPathtoseries,name);
                        end
                    movefile(newfilesinseriesHeader.Filename,newPathtoseriestype);
                end
                filesinseriesHeader = newfilesinseriesHeader;
            end
        end
       
     end
 dicomFilesnew=dir( [ rawdataDirectory '/*.dcm'] );
 nDicomfilesnew =length(dicomFilesnew);
 display('Files in raw data directory are sorted');
 diary('background');
else    
  dicomFilesnew=dir( [ rawdataDirectory '/*.dcm'] );
  nDicomfilesnew =length(dicomFilesnew);  
end
end
