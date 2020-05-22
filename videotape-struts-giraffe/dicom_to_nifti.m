function dicom_to_nifti( unsortedDicomDir, pathNifti )

%TODO: Output single multi-echo series as 4d nifti

mkdir(pathNifti);
disp(unsortedDicomDir);
disp(pathNifti);

% BEWARE: shell injection attacks here
if ispc == 1
    if system(['where dcm2niix']) == 1
        error 'dcm2niix is not installed.'
    end
else
    if system(['which dcm2niix']) == 1
        error 'dcm2niix is not installed.'
    end
end
participant = '';
system(['dcm2bids -d "' unsortedDicomDir '"' ' -o '  '"' pathNifti '"' ' -p '  '"' participant '"' ' -c '  'config.json']);

end


