%% Run unit test script

startup;

[testsResults] = runTestSuite('Unit')

for iTest = 1:length(testsResults) 
    if(~testsResults(iTest).Passed)
       quit(1); 
    end
end

quit(0);