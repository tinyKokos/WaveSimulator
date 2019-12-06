%% Central1DFiniteDiff(speed, deltaT, deltaX, funcAheadX, func, funcBehindX, funcBehindT)
%   speed - velocity of the wave 
%   deltaT - change in time between the points of the function
%   deltaX - change in space between the points of the function 
%   funcAheadX - f(x+1,t) 
%   func - f(x,t) 
%   funcBehindX - f(x-1,t) 
%   funcBehindT - f(x,t-1)
%
%   This is the update equation for the 1D standing wave with the Finite
%   Difference Method. The equation is based on the acceleration of the
%   wave in discrete points. This equation itself would have to be adjusted
%   for the initial run of the function from t0 to t1; for that process the
%   funcBehindT will be set to 0 and the output of this function must be
%   multiplied by 0.5 
function output2 = Central1DFiniteDiff(speed, deltaT, deltaX, ...
    funcAheadX, func, funcBehindX, funcBehindT)

%Check for valid arguments
arguments
   speed (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
   deltaT (1,:) {mustBeNumeric, mustBeFinite}
   deltaX (1,:) {mustBeNumeric, mustBeFinite}
   funcAheadX {mustBeNumeric, mustBeFinite}
   func {mustBeNumeric, mustBeFinite}
   funcBehindX {mustBeNumeric, mustBeFinite}
   funcBehindT {mustBeNumeric, mustBeFinite}
end

    output2 = (((speed^2)*(deltaT^2))/(deltaX^2))*(funcAheadX ...
    - 2*func + funcBehindX) + 2*func - funcBehindT; 
end