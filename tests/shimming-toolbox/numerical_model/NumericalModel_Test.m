classdef (TestTags = {'Simulation', 'Unit'}) NumericalModel_Test < matlab.unittest.TestCase

    properties
    end

    methods (TestClassSetup)
    end

    methods (TestClassTeardown)
    end

    methods (Test)
        function test_initiate_object_is_expected_class(testCase)

            testObj = NumericalModel();

            expectedClass = 'NumericalModel';
            actualClass = class(testObj);
            testCase.verifyEqual(actualClass, expectedClass)
        end
 
        function test_empty_initialization_returns_expected_starting_volume(testCase)
            testObj = NumericalModel();

            expectedVolume = zeros(128, 128);
            actualVolume = testObj.starting_volume;
            testCase.verifyEqual(actualVolume, expectedVolume)
        end
 
        function test_shepplogan_initialization_returns_expected_starting_volume(testCase)
            testObj = NumericalModel('Shepp-Logan');

            expectedVolume = phantom(128, 128);
            actualVolume = testObj.starting_volume;
            testCase.verifyEqual(actualVolume, expectedVolume)
        end
    end

end