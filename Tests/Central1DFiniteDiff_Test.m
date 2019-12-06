classdef Central1DFiniteDiff_Test < matlab.unittest.TestCase
    
    properties
        Property1
    end
    
    methods (Test)
        function testPositiveSlope(testCase)
            actSolution = Central1DFiniteDiff(1, 1, 1, 3, 2, 1, 1);
            expSolution = 3;
            testCase.verifyEqual(actSolution,expSolution);
    end
end

