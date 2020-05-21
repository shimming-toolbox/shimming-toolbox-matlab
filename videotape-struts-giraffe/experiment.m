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

exit;
