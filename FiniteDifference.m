%% FiniteDifference.m
% This function is written to simulate a standing wave in 1 dimensional
% space with the Finite Difference Method. 

%% Declared Input Values
close all; 

% Spacial Inputs
initialX = 0;
finalX = 3;
Npoints = 50; %number of points between the initial and final X
Ncycles = 1.5; %number of full sine cycles

%set inputs for the time
initialTime = 0;
finalTime = 5;
NtimePoints = 100; 

% Group Velocity of the Wave
PropagationSpeed = 300;

tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;

%filename of gif output of the program
filename = 'output.gif';

%delay value between each frame of the gif
FrameDelay = 0;

%% Stability Checking
% This is meant to check the stability of the Finite Difference Function
% across every update in time. This is done by using the
% Courant-Friedrichs-Lewy (CFL) Condition. Because this function is in 1
% dimension then the CFL number must be <= 1 or else the function is proven
% to be unstable. When the program finds out that the given values are
% unstable then the program will ask for a CFL constant and if a value is
% valid (i.e. <= 1) then the velocity of the function will be scaled to fit
% the given CFL from the user.
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

%% Calculate and Plot Function
% 
h = figure;
x = linspace(initialX,finalX,Npoints);

func = SineInput(Ncycles, finalX, initialX, Npoints);

nextFunc = zeros(length(func));
pastFunc = zeros(length(func));

for t = 1:(NtimePoints)
    %plot the new output
    plot(x, func);
    %This would have to be adjusted if amplitude is adjusted
    ylim([-2 2]);    
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if t == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end
    
    %Finite Difference Update Equation
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

    %update the past function to equal the old present function
    pastFunc = func; 
    %update the present function to equal the old future function
    func = nextFunc; 

    pause(FrameDelay);
end

function output2 = Central1DFiniteDiff(speed, deltaT, deltaX, ...
    funcAheadX, func, funcBehindX, funcBehindT)
%% Central1DFiniteDiff(speed, deltaT, deltaX, funcAheadX, func, funcBehindX, funcBehindT)
%   This is the update equation for the 1D standing wave with the Finite
%   Difference Method. The equation is based on the acceleration of the
%   wave in discrete points. This equation itself would have to be adjusted
%   for the initial run of the function from t0 to t1; for that process the
%   funcBehindT will be set to 0 and the output of this function must be
%   multiplied by 0.5 
%
%   speed - velocity of the wave 
%   deltaT - change in time between the points of the function
%   deltaX - change in space between the points of the function 
%   funcAheadX - This is the point that is one point ahead in space, but
%   not in time; this is typically represented as f(x+1,t) 
%   func - The current point that updated in the update equation;
%   represented as f(x,t) 
%   funcBehindX - This is the function that is behind one point in space
%   but not in time; represented as f(x-1,t) 
%   funcBehindT - This is the point in the function in the previous time 
%   iteration; represented as f(x,t-1)
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

function output = SineInput(Cycles, Xfinal, Xinitial, NumberOfPoints)
%% SineInput(Cycles, Xfinal, Xinitial, NumberOfPoints)
%   Creates a input for the Finite Difference update equation to use as the
%   initial condition of the function. This input is a simple sine wave
%   that is fixed to be 0 at the boundaries of the given space
%
%   Cycles - How many full cycles the output array must complete in the 
%   given amount of space
%   Xfinal - the value of the final point in space of the array
%   Xinitial - the value of the initial point in the space of the array
%   NumberOfPoints - how many discrete points the output array must have
%    
%   Example: SineInput(1, 6, 1, 4) = [0    0.9848   -0.3420   -0.8660]
arguments
   Cycles (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
   Xfinal (1,:) {mustBeNumeric, mustBeFinite}
   Xinitial (1,:) {mustBeNumeric, mustBeFinite}
   NumberOfPoints (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
end
x = linspace(Xinitial,Xfinal,NumberOfPoints);
output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
end