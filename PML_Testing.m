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

%PropagationSpeed = physconst('LightSpeed');
tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;

CFL = 0.75;

%check CFL condition 
%CFL = (PropagationSpeed*tDelta)/xDelta;

%Ideally I should just check the CFL and adjust some input to make it work
%but right now I will be lazy and just adjust the PropagationSpeed to fix
%my CFL

PropagationSpeed = (CFL*xDelta)/tDelta;

%New constants that I do not know what they are for the PML
alpha0 = 1;
c0 = 1;

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

function output3 = Central1DFiniteDiffPML(speed, c0, alpha, deltaT, ...
    deltaX, funcAheadX, func, funcBehindX, funcBehindT)
%Central1DFiniteDiffPML: <summary>
%   <function details>
constant1 = (deltaT^2*speed^2*c0)/(deltaX^2*c0 + (deltaX^2)*deltaT*(speed^2)*alpha);
constant2 = (alpha*deltaT*speed^2 - c0)/(c0+deltaT*alpha*speed^2);
constant3 = 2*c0*deltaX^2;
constant4 = -1*(deltaT^2)*(deltaX^2)*(speed^2)*(alpha^2)*c0;
constant5 = -2*(deltaT^2)*(speed^2)*(c0^2);
constant6 = c0*deltaX^2 + alpha*deltaT*(deltaX^2)*(speed^2);
output3 = constant1*funcAheadX + ...
    ((constant3+constant4+constant5)/constant6)*func + ...
    constant1*funcBehindX + constant2*funcBehindT;
end

function output = SineInput(Cycles,Xfinal, Xinitial, NumberOfPoints)
%SineInput Summary of this function goes here
%   Detailed explanation goes here
x = linspace(Xinitial,Xfinal,NumberOfPoints);
output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
end
