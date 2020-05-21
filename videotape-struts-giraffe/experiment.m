disp('Hello');

% hack: install the shimming-toolbox package
addpath(genpath('shimming-toolbox'))

data = 'data_testing/'

tmp = tempname
mkdir(tmp)
dicom_sorted = fullfile(tmp, 'dicom_sorted')
sortdicoms(fullfile(data, 'dicom_unsorted'), dicom_sorted)

disp(['Put results in ' dicom_sorted])
ls(dicom_sorted)
disp(['-----'])

exit;
