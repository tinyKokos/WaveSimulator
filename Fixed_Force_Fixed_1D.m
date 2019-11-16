%Fixed_Force_Fixed_1D by Adam Gleichman
% Fixed Force Fixed 1D String setup: created on 11/16/2019
% Test Ideas I should think about
%       Is there a way to break the source point?

close all;

initialX = 0;
finalX = 3;
Npoints = 100; %number of points between the initial and final X
sourcePt = floor((Npoints)/2)
frequency = 3; %frequency of sine for input of source point

%set inputs for the time
initialTime = 0;
finalTime = 1;
NtimePoints = 50; 

PropagationSpeed = 300;
tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;
x = linspace(initialX,finalX,Npoints);

func = zeros(Npoints);
nextFunc = zeros(Npoints);
pastFunc = zeros(Npoints);
func(sourcePt) = sin(frequency*2*pi*tDelta);

CFL = (PropagationSpeed*tDelta)/xDelta;

if CFL > 1
    fprintf('Your inputs will create an unstable system. Speed will be automatically adjusted for stability\n\n');
    prompt = 'Enter desired CFL: ';
    desiredCFL = input(prompt);
    if desiredCFL <= 1
        PropagationSpeed = (desiredCFL*xDelta)/tDelta;
    else
        fprintf('Not a valid CFL. Exiting Program\n');
        return
    end
end

for t = 1:(NtimePoints)
    %Note to self: the central Finite Diff function use in the loops
    % will be replaced with the Power Law Equation Function I made
        for n = 2:(sourcePt-1)
            if t == 1
                nextFunc(n)= 0.5*Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), 0);
            else
                nextFunc(n) = Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), pastFunc(n));
            end
        end
        for m = (sourcePt+1):(Npoints-1)
            if t == 1
                nextFunc(m)= 0.5*Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(m+1), func(m), func(m-1), 0);
            else
                nextFunc(m) = Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(m+1), func(m), func(m-1), pastFunc(m));
            end
        end
    nextFunc(sourcePt) = sin(frequency*2*pi*t*tDelta);
    
    %force the fixed boundary conditions
    % this might ultimately be redundant
    nextFunc(1) = 0;
    nextFunc(end) = 0;

    pastFunc = func; %update the past function to equal the old present function
    func = nextFunc; %update the present function to equal the old future function
    %plot the new output
    plot(x, func);
    ylim([-2 2]);
    pause(0.1); % I might want to modify this to be a different value
end

function output2 = Central1DFiniteDiff(speed, deltaT, deltaX, ...
    funcAheadX, func, funcBehindX, funcBehindT)
%Central1DFiniteDiff: solves Utt = (c^2)Uxx equation
%   This function solves Utt = (c^2)Uxx wave equation for each individual
%   point so this will require a for loop setup to calculate all points
%   numerically in the function. For the calculating the first set of
%   points after the initial function is given use funcBehindT = 0 and
%   half the result of the function
output2 = (((speed^2)*(deltaT^2))/(deltaX^2))*(funcAheadX ...
    - 2*func + funcBehindX) + 2*func - funcBehindT; 
end

%This is for later in the project
function output3 = PML_1dFiniteDiff(speed, alpha, deltaT, ...
    deltaX, funcAheadX, func, funcBehindX, funcBehindT, time)
%Central1DFiniteDiffPML: <summary>
%   This is using the numeric approxiation of the first derivative of the
%   function; (d/dt)f(x,t) -> (f(x,t) - f(x,t-tDELTA))/tDELTA
constant1 = ((speed^2)*(deltaT^2))/(deltaX^2);
constant2 = alpha*speed*deltaT;
if time == 0 %this funciton only works for t == initial time
    output3 = constant1*funcAheadX + ...
        (1 - 2*constant1 - constant2^2)*func + ...
        constant1*funcBehindX;
else %this function only works for t != initial time
    output3 = constant1*funcAheadX + ...
        (2 - 2*constant1 - constant2^2 - 2*constant2)*func + ...
        constant1*funcBehindX + (2*constant2 - 1)*funcBehindT;
end
end

