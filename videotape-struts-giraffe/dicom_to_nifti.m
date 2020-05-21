function dicom_to_nifti( unsortedDicomDir, pathNifti )

mkdir(pathNifti);
disp(unsortedDicomDir);
disp(unsortedDicomDir);
disp(unsortedDicomDir);
disp(pathNifti);
% BEWARE: shell injection attacks here
if system(['which dcm2niix']) == 1
  error 'dcm2niix is not installed.'
end
system(['dcm2niix -o ' pathNifti ' ' unsortedDicomDir]);

end


