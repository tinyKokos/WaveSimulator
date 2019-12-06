%% Fixed_Force_Fixed_1D.m
% This is meant to model a wave traveling from a source point to a fixed
% boundary. I realized that my old way of adjusting the CFL when the CFL
% condition was not met by changing the speed of the propagating wave was
% creating the situation where the wave was never meeting the boundary so I
% do not want to publish this as complete software yet.

close all;

filename = 'FixedForceFixed1D.gif';
FrameDelay = 0;
ScaleChoice = 1;

initialX = 0;
finalX = 3;
Npoints = 200; %number of points between the initial and final X
sourcePt = floor((Npoints)/2);
frequency = 3; %frequency of sine for input of source point

%set inputs for the time
initialTime = 0;
finalTime = 20;
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

%Modifying the propagation speed to fix the CFL makes working with this
%simulation pretty frustrating. I should try some other way.
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
                
            case 3
                
            otherwise
                
        end
    else
        fprintf('Not a valid CFL. Exiting Program\n');
        return
    end
end
h = figure;
for t = 1:NtimePoints
    %plot the new output
    plot(x, func);
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
      for n = 2:(sourcePt-1)
            %NOTE TO SELF: the central Finite Diff function use in the loops
            % will be replaced with the Power Law Equation Function I made
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
    %(this might ultimately be redundant)
    nextFunc(1) = 0;
    nextFunc(end) = 0;

    pastFunc = func; %update the past function to equal the old present function
    func = nextFunc; %update the present function to equal the old future function
    pause(FrameDelay);
end

function output2 = Central1DFiniteDiff(speed, deltaT, deltaX, ...
    funcAheadX, func, funcBehindX, funcBehindT)
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