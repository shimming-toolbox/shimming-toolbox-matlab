function dicom_to_nifti( sortedDicomDir, pathNifti )

%TODO: Output single multi-echo series as 4d nifti

mkdir(pathNifti);
disp(sortedDicomDir);
disp(pathNifti);

folders = dir(sortedDicomDir);
for iFolder = 3:length(folders)
    % BEWARE: shell injection attacks here
    if system(['which dcm2niix']) == 1
      error 'dcm2niix is not installed.'
    end
    subFolderNifti = fullfile(pathNifti, folders(iFolder).name);
    subFolderDicom = fullfile(sortedDicomDir, folders(iFolder).name);
    mkdir(subFolderNifti);
    system(['dcm2niix -b y -a y -o "' subFolderNifti '" "' subFolderDicom '"']);
    
end

end


