% Test functions dealing with file management, copy, dicom conversion, etc.

function test_suite=test_file_management
    test_functions=localfunctions();
    initTestSuite;

function test_dicomsorts
%     TODO: the function below requires dicominfo, which can be installed
%      via Octave's package: pkg install -forge dicom.
%      For now, the lines below are commented until the required package is
%      installed.
%     sortdicoms('tests/data_testing/dicom_unsorted/', 'tests/data_testing/dicom_sorted/');
%     assert(isfile('tests/data_testing/dicom_sorted/06_a_gre_DYNshim/echo_11.5/acdc_95p-HC7;NC1,2-0001-0001.dcm'));
