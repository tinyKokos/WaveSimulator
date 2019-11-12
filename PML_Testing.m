%PML Testing by Adam Gleichman
% started from 11/7/2019
% started by copying the FiniteDifference.m file on 11/7/2019

initialX = 0;
finalX = 3;
Npoints = 10; %number of points between the initial and final X
Ncycles = 1; %number of full sine cycles

%set inputs for the time
initialTime = 0;
finalTime = 1;
NtimePoints = 50; 

PropagationSpeed = 300;
tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;

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

%Must make it be an array that matches the size of the data set
alphaPoints = 10;
alphaMax = 1; %can't allow negatives

for n = 1:alphaPoints

if (2*alphaPoints) > Npoints
   fprintf("Thickness of PML is too large\n");
   return
elseif (2*alphaPoints) == Npoints
    fprintf("Warning: PML covers entire space\n");
end

%plot the initial function set
func = SineInput(Ncycles, finalX, initialX, Npoints);
x = linspace(initialX,finalX,Npoints);
plot(x, func);
ylim([-2 2]);

nextFunc = zeros(length(func));
pastFunc = zeros(length(func));

for t = 1:(NtimePoints-1)
        for n = 2:(Npoints-1)
            if t == 1
                nextFunc(n)= 0.5*Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), 0);
            else
                nextFunc(n) = Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), pastFunc(n));
            end
        end
 
    %force the fixed boundary conditions
    nextFunc(1) = 0;
    nextFunc(end) = 0;

    pastFunc = func; %update the past function to equal the old present function
    func = nextFunc; %update the present function to equal the old future function
    %plot the new output
    plot(x, func);
    ylim([-2 2]);
    pause(0.1);
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

function output3 = PML_1dFiniteDiff(speed, alpha, deltaT, ...
    deltaX, funcAheadX, func, funcBehindX, funcBehindT)
%Central1DFiniteDiffPML: <summary>
%   This is using the numeric approxiation of the first derivative of the
%   function; (d/dt)f(x,t) -> (f(x,t) - f(x,t-tDELTA))/tDELTA
constant1 = ((speed^2)*(deltaT^2))/(deltaX^2);
constant2 = alpha*speed*deltaT;
%this function only works for t != initial time
output3 = constant1*funcAheadX + ...
    (2 - 2*constant1 - constant2^2 - 2*constant2)*func + ...
    constant1*funcBehindX + (2*constant2 - 1)*funcBehindT;
end

function output = SineInput(Cycles,Xfinal, Xinitial, NumberOfPoints)
%SineInput Summary of this function goes here
%   Detailed explanation goes here
x = linspace(Xinitial,Xfinal,NumberOfPoints);
output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
end
