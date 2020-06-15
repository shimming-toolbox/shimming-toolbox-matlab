function dicom_to_nifti( unsortedDicomDir, niftiPath )

%TODO: Output single multi-echo series as 4d nifti

mkdir(niftiPath);
disp(unsortedDicomDir);
disp(niftiPath);

% Check for which "which" to use depending on OS
if ispc == 1
    which = 'where';
else
    which = 'which';
end

% Make sur dcm2niix is installed
if system([which ' dcm2niix']) ~= 0
    error 'dcm2niix is not installed.'
end
% Make sure dcm2bids is installed
if system([which ' dcm2bids']) ~= 0
    error 'dcm2bids is not installed.'
end

% Create bids structure for data
participant = '';
if system(['dcm2bids_scaffold -o ' niftiPath]) ~= 0
    error 'dcm2bids_scaffold'
end

% Add original data to niftiPath/sourcedata
if system(['cp -r ' unsortedDicomDir ' ' fullfile(niftiPath,'sourcedata')]) ~= 0
    error 'copy'
end

% Call the dcm2bids_helper
if system(['dcm2bids_helper -d ' unsortedDicomDir ' -o ' niftiPath]) ~= 0
    error 'dcm2bids_helper'
end

% Check if there is data
helperPath = fullfile(niftiPath,'tmp_dcm2bids','helper');
if ~isfolder(helperPath)
    error 'dcm2bids_helper could not create directory helper'
end

% Make sure there is data in niftiPath/tmp_dcm2bids/helper
helperfileList = dir(helperPath);
if isempty(helperfileList)
    error 'No data to process'
end

% Create list of acquisitions
acquisitionNames = {};
acquisitionNumbers = {};
iAcq = 0;

% Create list containing all files
helperfileList =  helperfileList(~ismember({helperfileList.name},{'.','..'}));
for iFile = 1:length(helperfileList)
    
    % If it's a json file
    [filePath,name,ext] = fileparts(fullfile(helperfileList(iFile).folder, helperfileList(iFile).name));
    if strcmp(ext,'.json')
        % Check for both .gz and .nii
        niftiFile =  helperfileList(ismember({helperfileList.name},{[name '.nii'],[name '.nii.gz']}));
        if length(niftiFile) ~= 1
            error('Did not find linked nii file')
        end
        
        % Read json file
        [~,~,jsonInfo] = imutils.read_nii(fullfile( filePath , niftiFile.name ));
        iAcq = iAcq + 1;
        
        % Create future folder name
        acquisitionNumbers{iAcq,1} = sprintf( '%03d', jsonInfo.SeriesNumber ) ;
        acquisitionNames{iAcq,1} = jsonInfo.SeriesDescription;
        % *******Modality could be acquisition name ********
        modality{iAcq,1} = jsonInfo.Modality;
    end
end

% Remove duplicates (stable is specified to keep the same order)
[acquisitionNumbers,ia] = unique(acquisitionNumbers, 'stable');
acquisitionNames = acquisitionNames(ia);
modality = modality(ia);

% For every acqs,
outputDir = fullfile(niftiPath, 'code');

clear iAcq
for iAcq = 1:length(acquisitionNames)
%     create config file, place in niftiPath/code
    configFilePath = createConfig(outputDir, acquisitionNumbers{iAcq}, acquisitionNames{iAcq}, modality{iAcq});

%     call dcm2bids
    if system(['dcm2bids -d "' unsortedDicomDir '"' ' -o '  '"' niftiPath '"' ' -p '  '"' participant '"' ' -c '  configFilePath]) ~= 0
      error 'dcm2bids failed'
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filePath = createConfig(outputDir, acquisitionNumber, acquisitionName, modality)
    
% Create names
name = [acquisitionNumber '_' acquisitionName];
ext = '.json';
filePath = fullfile(outputDir,[name ext]);
    
% Create file
fid = fopen(filePath,'w');
if fid == -1
    error 'could not create config file'
end  

% Write to file
fprintf(fid,'{\n');
fprintf(fid,'    "searchMethod": "fnmatch",\n');
fprintf(fid,'    "defaceTpl": "pydeface --outfile {dstFile} {srcFile}",\n');
fprintf(fid,'    "descriptions": [\n');
fprintf(fid,'        {\n');
fprintf(fid,'            "dataType": "%s",\n',name);
fprintf(fid,'            "SeriesDescription": "%s",\n',name);
fprintf(fid,'            "modalityLabel": "%s",\n', modality);
fprintf(fid,'            "criteria": {\n');
fprintf(fid,'                "SeriesDescription": "%s",\n', acquisitionName);
fprintf(fid,'                "SeriesNumber": "%s"\n', num2str(str2num(acquisitionNumber)));
fprintf(fid,'              }\n');
fprintf(fid,'        }\n');
fprintf(fid,'    ]\n');
fprintf(fid,'}\n');

% TODO: Add criteria in each file to add phase and magnitude seperation?

% Close file
fclose(fid); 

end
