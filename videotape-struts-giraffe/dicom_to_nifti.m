function dicom_to_nifti( unsortedDicomDir, pathNifti )

mkdir(pathNifti);
disp(unsortedDicomDir);
disp(unsortedDicomDir);
disp(unsortedDicomDir);
disp(pathNifti);
% BEWARE: shell injection attacks here
system(['dcm2niix -o ' pathNifti ' ' unsortedDicomDir]);

end


