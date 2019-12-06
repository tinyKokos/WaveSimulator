%% SineInput(Cycles, Xfinal, Xinitial, NumberOfPoints)
%   Cycles - How many full cycles the output array must complete in the 
%   given amount of space
%   Xfinal - the value of the final point in space of the array
%   Xinitial - the value of the initial point in the space of the array
%   NumberOfPoints - how many discrete points the output array must have
%    
%   Creates a input for the Finite Difference update equation to use as the
%   initial condition of the function. This input is a simple sine wave
%   that is fixed to be 0 at the boundaries of the given space
function output = SineInput(Cycles, Xfinal, Xinitial, NumberOfPoints)

%Check for valid arguments
arguments
   Cycles (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
   Xfinal (1,:) {mustBeNumeric, mustBeFinite}
   Xinitial (1,:) {mustBeNumeric, mustBeFinite}
   NumberOfPoints (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
end

    x = linspace(Xinitial,Xfinal,NumberOfPoints);
    output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
end