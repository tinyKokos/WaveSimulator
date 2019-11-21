%PML Testing by Adam Gleichman
% started from 11/7/2019
% started by copying the FiniteDifference.m file on 11/7/2019

%Notes to self:
% 1st - switch this to a Fixed Force Fixed String setup
% 2nd - create the alpha array thing and finish the PML
% 3rd - modify the automatic CFL adjust to be better
close all;

initialX = 0;
finalX = 3;
Npoints = 100; %number of points between the initial and final X
sourcePt = floor((Npoints)/2);
%Ncycles = 1; %number of full sine cycles for fixed string
frequency = 3; %frequency of sine for input of source point

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

%Alpha Stuff
%Must make it be an array that matches the size of the data set
alphaPoints = 10; %this is the amount of points on each side of the graph
alphaMax = 1; %can't allow negatives
alphaDummy = zeros(alphaPoints);
alpha = zeros(Npoints);
limit1 = length(alpha) - alphaPoints;

if (2*alphaPoints) >= Npoints
   fprintf("Thickness of PML is too large\n");
   return
else
    for e = 1:alphaPoints
       alphaDummy(e) = 0.25*e; 
    end
    
    for r = 1:(alphaPoints-1)
        alpha(r) = alphaDummy(end-r);
    end
    alpha(limit:Npoints) = alphaDummy;
    alpha(sourcePt) = 0; %force alpha to not affect the source point
end


%plot the initial function set
%func = SineInput(Ncycles, finalX, initialX, Npoints);
x = linspace(initialX,finalX,Npoints);
%plot(x, func);
%ylim([-2 2]);

func = zeros(Npoints);
nextFunc = zeros(Npoints);
pastFunc = zeros(Npoints);
func(sourcePt) = sin(frequency*2*pi*tDelta);
for t = 1:(NtimePoints)
    %func(sourcePt) = sin(frequency*2*pi*t*tDelta);
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
    deltaX, funcAheadX, func, funcBehindX, funcBehindT, time)
%Central1DFiniteDiffPML: <summary>
%   This is using the numeric approxiation of the first derivative of the
%   function; (d/dt)f(x,t) -> (f(x,t) - f(x,t-tDELTA))/tDELTA
constant1 = ((speed^2)*(deltaT^2))/(deltaX^2);
constant2 = alpha*speed*deltaT;
if time == 0 %this function only works for t == initial time
    output3 = constant1*funcAheadX + ...
        (1 - 2*constant1 - constant2^2)*func + ...
        constant1*funcBehindX;
else %this function only works for t != initial time
    output3 = constant1*funcAheadX + ...
        (2 - 2*constant1 - constant2^2 - 2*constant2)*func + ...
        constant1*funcBehindX + (2*constant2 - 1)*funcBehindT;
end
end

function output = SineInput(Cycles,Xfinal, Xinitial, NumberOfPoints)
%SineInput Summary of this function goes here
%   Detailed explanation goes here
x = linspace(Xinitial,Xfinal,NumberOfPoints);
output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
end
