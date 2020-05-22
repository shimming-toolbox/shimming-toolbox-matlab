disp('Hello');

% hack: install the shimming-toolbox package
addpath(genpath('..'))

data = 'data_testing/'
% TODO: check if this folder exists and prompt to download (or maybe just download it directly)

tmp = tempname
mkdir(tmp)
dicom_sorted = fullfile(tmp, 'dicom_sorted')
sortdicoms(fullfile(data, 'dicom_unsorted'), dicom_sorted)

disp(['Put results in ' dicom_sorted])
ls(dicom_sorted)
disp(['-----'])

tmp = tempname
mkdir(tmp)
nifti_path = fullfile(tmp, 'niftis')
dicom_to_nifti(dicom_sorted, nifti_path)
disp(['-----'])
disp(nifti_path)
ls(nifti_path)

list = dir(fullfile(nifti_path, ['**' filesep '*.nii']));
imgs  = cell(length(list),1) ;
infos = cell(length(list),1) ;
jsons = cell(length(list),1) ;
for iImg = 1 : length( list )
   [imgs{iImg}, infos{iImg}, jsons{iImg}] = img.read_nii( ...
       fullfile( list(iImg).folder , list(iImg).name ) );
end

    
    
% exit;
