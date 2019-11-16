% FiniteDifference2D by Adam Gleichman
% Lastest Version saved on 11/16/2019
% Desciption: this is the general 2D simulation of a wave that is the
% product of a sine in the x and in the y direction.
close all;

initialX = 0;
finalX = 3;
initialY = 0;
finalY = 3;
Npoints = 10; %number of points between the initial and final X, Y
Ncycles = 1.5; %number of full sine cycles

%set inputs for the time
initialTime = 0;
finalTime = 1;
NtimePoints = 50; 

PropagationSpeed = 300; %this is just meant as a place holder value
PropagationSpeedX = PropagationSpeed; %I assume that speed in Y & X
PropagationSpeedY = PropagationSpeed; % must be equal to each other
tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;
yDelta = (finalY-initialY)/Npoints;

x = linspace(initialX,finalX,Npoints);
y = linspace(initialY,finalY,Npoints);

func = SineInput2D(Ncycles, finalX, initialX, finalY, initialY, Npoints);
mesh(x, y, func);
zlim([-2 2]);

nextFunc = zeros(length(x), length(y));
pastFunc = zeros(length(x), length(y));

%check CFL condition 
CFL = (PropagationSpeed*tDelta)/xDelta;

if CFL > 1
    fprintf('Your inputs will create an unstable system. Speed will be automatically adjusted for stability\n\n');
    prompt = 'Enter desired CFL: ';
    desiredCFL = input(prompt);
    if desiredCFL <= 1
        PropagationSpeedX = (desiredCFL*xDelta)/tDelta;
        PropagationSpeedY = (desiredCFL*yDelta)/tDelta;
    else
        fprintf('Not a valid CFL. Exiting Program\n');
        return
    end
end

for t = 1:(NtimePoints-1)
   for n = 2:(Npoints-1) %correspond to X
      for m = 2:(Npoints-1) %correspond to Y
          if t == 1
              nextFunc(n,m) = 0.5*Central2DFiniteDiff(PropagationSpeedX, PropagationSpeedY, ...
                  tDelta, xDelta, yDelta, func(n+1,m), func(n,m+1), func(n,m), ...
                  func(n-1,m), func(n,m-1), 0);
          else
              nextFunc(n,m) = Central2DFiniteDiff(PropagationSpeedX, PropagationSpeedY, ...
                  tDelta, xDelta, yDelta, func(n+1,m), func(n,m+1), func(n,m), ...
                  func(n-1,m), func(n,m-1), pastFunc(n,m));
          end
          nextFunc(1,m) = 0;
          nextFunc(end,m) = 0;
      end
      nextFunc(n,1) = 0;
      nextFunc(n,end) = 0;
   end
   pastFunc = func; %update the past function to equal the old present function
   func = nextFunc; %update the present function to equal the old future function
   
   mesh(x, y, func);
   zlim([-2 2]);
   pause(0.1);
end

function output4 = Central2DFiniteDiff(speedX, speedY, deltaT, deltaX, ...
    deltaY, funcAheadX, funcAheadY, func, funcBehindX, funcBehindY, funcBehindT)
%Central1DFiniteDiff Summary of this function goes here
%   Detailed explanation goes here

%Note to self: this equation might be correct
output4 = (((speedX^2)*(deltaT^2))/(deltaX^2))*(funcAheadX ...
    - 2*func + funcBehindX) + (((speedY^2)*(deltaT^2))/(deltaY^2))*(funcAheadY ...
    - 2*func + funcBehindY) + 2*func - funcBehindT; 

end

function output3 = SineInput2D(Cycles, Xfinal, Xinitial, Yfinal, Yinitial, NumberOfPoints)
%SineInput2D Summary of this function goes here
%   NOTE: this does assume that the number of points is equal
%       in the y direction and the x direction
x = linspace(Xinitial,Xfinal,NumberOfPoints);
y = linspace(Yinitial,Yfinal,NumberOfPoints);
output3 = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial)).*sin(((Cycles*2*pi)/Yfinal)*(y'-Yinitial));
end

%Archive
%{

%Test of central difference 1D function
% used in the main loop
% (this is archived just in case)
%{
n= 6; % this to crunch the normal equation

%for time = 0 equation
%nextFunc(n) = (((PropagationSpeed^2)*(tDelta^2))...
%                /(2*(xDelta^2)))*(func(n+1)-2*func(n)+func(n-1))...
%                + func(n);

%for all the other times
nextFunc(n) = (((PropagationSpeed^2)*(tDelta^2))...
                    /(xDelta^2))*(func(n+1)-2*func(n)+func(n-1))...
                    + 2*func(n) - pastFunc(n);

%for time = 0                
%something = 0.5*CentralFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), 0);

%for time not = 0
something = CentralFiniteDiff(PropagationSpeed, tDelta, xDelta, func(n+1), func(n), func(n-1), pastFunc(n));

if something == nextFunc(n)
    disp('Function checks out');
else
    disp('Function doesnt checks out');
    disp('Function =%d\n', something);
    disp('Equation =%d\n', nextFunc(n));
end
%}
%}