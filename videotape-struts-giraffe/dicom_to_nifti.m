function dicom_to_nifti( unsortedDicomDir, niftiPath )

%TODO: Output single multi-echo series as 4d nifti

mkdir(niftiPath);
disp(unsortedDicomDir);
disp(niftiPath);

if ispc == 1
    which = 'where'
else
    which = 'which'
end
if system([which ' dcm2niix']) ~= 0
    error 'dcm2niix is not installed.'
end
if system([which ' dcm2bids']) ~= 0
    error 'dcm2bids is not installed.'
end

% BEWARE: shell injection attacks here
participant = '';
if system(['dcm2bids -d "' unsortedDicomDir '"' ' -o '  '"' niftiPath '"' ' -p '  '"' participant '"' ' -c '  'config.json']) ~= 0
  error 'dcm2bids failed'
end

end


