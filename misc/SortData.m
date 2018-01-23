%==========================================================================
%Function to sort the Raw data from the scan 
%==========================================================================

function SortData(rawdataDirectory,sortdataDirectory)

%Dicom file inside the directory ------------------------------------------
dicomfiles=dir( [ rawdataDirectory '/*.dcm'] );
ndicomfiles =length(dicomfiles);
%Declare initial parameters -----------------------------------------------
SeriesDescription=[];
ndicomfilesnew=0;
%Writing a feedback in a textfile background------------------------------- 
diary('background');

if ~exist(sortdataDirectory,'dir'), mkdir(sortdataDirectory); end

%Read Dicom files header and extract series name and number ---------------
    files_hdr_prec = dicominfo( [rawdataDirectory '/' dicomfiles(1).name]);
    series_name=files_hdr_prec.SeriesDescription;
    foldernumber=files_hdr_prec.SeriesNumber;
    nseriesdescription=1;
    if (foldernumber < 10)
    foldername= strcat('0',num2str(foldernumber),'_',series_name);
    else
    foldername= strcat(num2str(foldernumber),'_',series_name);
    end
    SeriesDescription{nseriesdescription}=foldername;
  
    new_path = fullfile (sortdataDirectory, foldername);

   if ~exist(new_path,'dir')
      %Create new directories for each series name-------------------------
      mkdir(sortdataDirectory,foldername);
   end
      %Copy the dicom files in the new directory---------------------------
    copyfile(files_hdr_prec.Filename,new_path);
    

while(1)   
if(ndicomfilesnew~=ndicomfiles)
   display('Sorting in process');
   diary('background');
   %Check the presence for new files in the raw data directory-------------
   newfiles=abs(ndicomfilesnew-ndicomfiles);
   dicomfiles=dir( [ rawdataDirectory '/*.dcm'] );
        if(ndicomfilesnew==0)
            a=2;
            b=ndicomfiles;
        else
            a=ndicomfiles+1;
            b=ndicomfiles+newfiles;
        end
    ndicomfiles =length(dicomfiles);
    for iscan=a:b
        files_hdr = dicominfo( [rawdataDirectory '/' dicomfiles(iscan).name] ) ;
    
        if strcmp(files_hdr.SeriesNumber,files_hdr_prec.SeriesNumber)==1
              copyfile(files_hdr.Filename,new_path);
        else
            series_name2=files_hdr.SeriesDescription;
            foldernumber2=files_hdr.SeriesNumber;
                if (foldernumber2 < 10)
                foldername= strcat('0',num2str(foldernumber2),'_',series_name2);
                else
                 foldername= strcat(num2str(foldernumber2),'_',series_name2);
                end
            
              if ~any(strcmp(SeriesDescription,foldername));
                    nseriesdescription=nseriesdescription+1;
                    SeriesDescription{nseriesdescription}= foldername;
                    new_path = fullfile (sortdataDirectory, foldername);
                
                        if ~exist(new_path,'dir')
                        
                        mkdir(sortdataDirectory,foldername);
                        end
               end
               copyfile(files_hdr.Filename,new_path);
               %Comparison with the precedent dicom files sorted-----------
               files_hdr_prec = files_hdr;
         end
 
    end

    for i= 1:nseriesdescription
        new_pathtoseries = fullfile (sortdataDirectory, SeriesDescription{i});
        %Dicom files inside the new directories----------------------------
        new_dicomfilesinseries=dir( [new_pathtoseries '/*.dcm'] );
        nfilesinseries=length(new_dicomfilesinseries);
    
        if ~isempty(new_dicomfilesinseries) 
        
        filesinseries_hdr_prec = dicominfo( [new_pathtoseries '/' new_dicomfilesinseries(1).name]);
        ImageType=filesinseries_hdr_prec.ImageType;
        name=fullfile(ImageType,num2str(filesinseries_hdr_prec.EchoNumber));
        new_pathtoseriestype = fullfile (new_pathtoseries,name);
        if ~exist(new_pathtoseriestype,'dir')
        mkdir(new_pathtoseries,name);
        end
        movefile(filesinseries_hdr_prec.Filename,new_pathtoseriestype);
    
            for j=2:nfilesinseries
                 filesinseries_hdr = dicominfo( [new_pathtoseries '/' new_dicomfilesinseries(j).name]);
            
                if strcmp(filesinseries_hdr.ImageType,files_hdr_prec.ImageType)==1
                    movefile(filesinseries_hdr.Filename,new_pathtoseriestype);
                else
                    foldertype=filesinseries_hdr.ImageType;
                    name=fullfile(foldertype,num2str(filesinseries_hdr_prec.EchoNumber));
                    new_pathtoseriestype = fullfile (new_pathtoseries, name);
                        if ~exist(new_pathtoseriestype,'dir')
                        mkdir(new_pathtoseries,name);
                        end
                    movefile(filesinseries_hdr.Filename,new_pathtoseriestype);
                end
                filesinseries_hdr_prec = filesinseries_hdr;
            end
        end
       
     end
 dicomfilesnew=dir( [ rawdataDirectory '/*.dcm'] );
 ndicomfilesnew =length(dicomfilesnew);
 display('Files in raw data directory are sorted');
 diary('background');
else    
  dicomfilesnew=dir( [ rawdataDirectory '/*.dcm'] );
  ndicomfilesnew =length(dicomfilesnew);  
end
end
