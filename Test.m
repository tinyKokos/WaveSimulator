%% Fixed_Force_Fixed_1D.m
% This is meant to model a wave traveling from a source point to a fixed
% boundary. I realized that my old way of adjusting the CFL when the CFL
% condition was not met by changing the speed of the propagating wave was
% creating the situation where the wave was never meeting the boundary so I
% do not want to publish this as complete software yet.

close all;

%% Declare Constants
%name for the gif file output
filename = 'FixedForceFixed1D.gif';
%create variable for the figure for saving the output into a gif
h = figure;
%frame delay between each frame of the simulation t points
FrameDelay = 0.025;
%this is meant to be a default for the CFL check prompt
ScaleChoice = 1;

%set inputs for the time
initialTime = 0;
finalTime = 8;
NtimePoints = 800;

initialX = 0;
finalX = 3;
Npoints = 200; %number of points between the initial and final X

%wave properties
PropagationSpeed = 1;
frequency = 1; %frequency of sine for input of source point
%% Create Intermediate Constants

sourcePt = floor((Npoints)/2);

%percent of the total time points that the source is creating a point
SourceDuration = 0.125*NtimePoints;

tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;
x = linspace(initialX,finalX,Npoints);

func = zeros(1,Npoints);
nextFunc = zeros(1,Npoints);
pastFunc = zeros(1,Npoints);

func(sourcePt) = sin(frequency*2*pi*tDelta);
%% Create Attenuation Array
alpha = zeros(1,Npoints);
endLeft = sourcePt-1;
beginRight = sourcePt+2;

alpha_x = x(1:endLeft);
alpha_model = alpha_x;
alpha(1:endLeft) = flip(alpha_model);
alpha(beginRight:end) = alpha_model;
%% CFL Checking
% This is currently a work in progress. The program corrects the CFL by
% changing the propagation speed of the wave to meet the CFL condition but
% I noticed that because the propagation speed is changed such that the
% wave never meets the boundary so I am not sure if the code has the
% correct behaivor at the boundary.
CFL = (PropagationSpeed*tDelta)/xDelta;
if CFL > 1
    fprintf('Your inputs will create an unstable system. Speed will be automatically adjusted for stability\n\n');
    prompt = 'Enter desired CFL: ';
    desiredCFL = input(prompt);
    if desiredCFL <= 1
        prompt = 'Choose to scale Speed (1), Time(2), or X(3): ';
        ScaleChoice = input(prompt);
        switch ScaleChoice
            case 1
                PropagationSpeed = (desiredCFL*xDelta)/tDelta;
            case 2
                %remember that smaller tDelta is the goal for stability
                tDelta = (desiredCFL*xDelta)/PropagationSpeed;
                NtimePoints = ceil((finalTime-initialTime)/tDelta);
            case 3
                %remember that bigger xDelta is the goal for stability
                xDelta = (PropagationSpeed*tDelta)/desiredCFL;
                Npoints = floor((finalX-initialX)/xDelta);
            otherwise
                %input is not valid so exit the program
                printf('This input is not valid. Please read the prompt\n');
                return
        end
    else
        %Given CFL number is incorrect; exit the program
        fprintf('Not a valid CFL. Exiting Program\n');
        return
    end
end
%% Simulation Loop
for t = 1:NtimePoints
    %% Plot Function
    %plot the new output
    plot(x, func, x, alpha);
    ylim([-2 2]);
    %% Save Result for .GIF
    %{
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
    %}
    %% Update Source Point
    if t <= SourceDuration 
        %% Update Equation Left of Source 
        for n = 2:(sourcePt-1)
            %NOTE TO SELF: the central Finite Diff function use in the loops
            % will be replaced with the Power Law Equation Function I made
            if t == 1
                nextFunc(n)= 0.5*Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), 0);
            else
                nextFunc(n) = Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), pastFunc(n));
            end
        end
        %% Update Equation Right of Source
        for m = (sourcePt+1):(Npoints-1)
            if t == 1
                nextFunc(m)= 0.5*Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(m+1), func(m), func(m-1), 0);
            else
                nextFunc(m) = Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(m+1), func(m), func(m-1), pastFunc(m));
            end
        end
        nextFunc(sourcePt) = sin(frequency*2*pi*t*tDelta);
    else
        %nextFunc(sourcePt) = 0;
        for u = 2:(Npoints-1)
           if t == 1
               nextFunc(u) = 0.5*Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(u+1), func(u), func(u-1), 0);
           else
               nextFunc(u) = Central1DFiniteDiff(PropagationSpeed, tDelta, xDelta, func(u+1), func(u), func(u-1), pastFunc(u));
           end
        end
    end
    
    %% Force Boundary Condition
    nextFunc(1) = 0;
    nextFunc(end) = 0;

    %% Swap Function Arrays
    pastFunc = func; %update the past function to equal the old present function
    func = nextFunc; %update the present function to equal the old future function
    pause(FrameDelay);
end
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
   funcAheadX (:,:) {mustBeNumeric, mustBeFinite}
   func (:,:) {mustBeNumeric, mustBeFinite}
   funcBehindX (:,:) {mustBeNumeric, mustBeFinite}
   funcBehindT (:,:) {mustBeNumeric, mustBeFinite}
end
output2 = (((speed^2)*(deltaT^2))/(deltaX^2))*(funcAheadX ...
    -2*func+funcBehindX)+2*func-funcBehindT; 
end