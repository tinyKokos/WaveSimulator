% test different speeds, boundary settings, and cycles
% try to undo the CFL setting the speed of the wave

% Then go back to 1D to add PML (Vaughn and McGough have the code) to 
%the 1d code
% idea is to replace my boundary with a PML (perfectly matched layer)
%(this is with the power law wave equation; it will be #12 in the paper
%that McGough will email you)

%Objective: start writing this code to work for 2D. 
% new input for the 2D (Amplitude)sin(Kx)sin(Ky)
% Kx = (n pi) / (x max)
% Ky = (n pi) / (ymax)

initialX = 0;
finalX = 3;
initialY = 0;
finalY = 3;
Npoints = 10; %number of points between the initial and final X, Y
Ncycles = 1.5; %number of full sine cycles

%set inputs for the time
initialTime = 0;
finalTime = 1;
NtimePoints = 200; 

%PropagationSpeed = physconst('LightSpeed');
tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;
yDelta = (finalY-initialY)/Npoints;

CFL = 0.75;

%check CFL condition 
%CFL = (PropagationSpeed*tDelta)/xDelta;

%Ideally I should just check the CFL and adjust some input to make it work
%but right now I will be lazy and just adjust the PropagationSpeed to fix
%my CFL

%for 2d the speed of X and Y should be scaled depending on
% delta in X compared with delta in Y

PropagationSpeedX = 0.5*(CFL*xDelta)/tDelta;
PropagationSpeedY = 0.5*(CFL*yDelta)/tDelta;

x = linspace(initialX,finalX,Npoints);
y = linspace(initialY,finalY,Npoints);

func = SineInput2D(Ncycles, finalX, initialX, finalY, initialY, Npoints);
mesh(x, y, func);
zlim([-2 2]);

nextFunc = zeros(length(x), length(y));
pastFunc = zeros(length(x), length(y));

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

%rest of the code
%{
%{
desiredCFL = 0.75;
if CFL > 1
    %this is my fun trying to make a dialog log box setup to correct the
    %CFL
    %{
    prompt = {'Your inputs will create an unstable system. It is not recommended to proceed with your setup. Would you prefer to run anyways (1), exit (2), or automatically set a new space interval to avoid this (3)?'};
    dlgtitle = 'Potential Unstability Detected';
    dims = [1 35];
    definput = {'2'};
    answer = inputdlg(prompt,dlgtitle,dims,definput)
    switch(answer)
        case '1'
            continue;
        case '3'
            prompt2 = {'Set your desired CFL (must be less than or equal to 1)'};
            dlgtitle2 = 'Set CFL';
            definput = {'0.75'};
            answer2 = inputdlg(prompt2, dlgtitle2,dims, definput);
            
        otherwise
            quit(1); 
%}
    
    %My hack way to force the CFL to meet a desired CFL
    xDelta = desiredCFL/(PropagationSpeed*tDelta);
    Npoints = round((finalX-initialX)/xDelta);
    disp("Note: CFL was not less than or equal to 1. Spatial dimensions are adjusted to Npoints = " + num2str(Npoints));
end
%}
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
%}

function output4 = Central2DFiniteDiff(speedX, speedY, deltaT, deltaX, ...
    deltaY, funcAheadX, funcAheadY, func, funcBehindX, funcBehindY, funcBehindT)
%Central1DFiniteDiff Summary of this function goes here
%   Detailed explanation goes here

%Note to self: this equation might be correct
output4 = (((speedX^2)*(deltaT^2))/(deltaX^2))*(funcAheadX ...
    - 2*func + funcBehindX) + (((speedY^2)*(deltaT^2))/(deltaY^2))*(funcAheadY ...
    - 2*func + funcBehindY) + 2*func - funcBehindT; 

end

function output2 = Central1DFiniteDiff(speed, deltaT, deltaX, ...
    funcAheadX, func, funcBehindX, funcBehindT)
%Central1DFiniteDiff Summary of this function goes here
%   Detailed explanation goes here
output2 = (((speed^2)*(deltaT^2))/(deltaX^2))*(funcAheadX ...
    - 2*func + funcBehindX) + 2*func - funcBehindT; 
end

function output3 = SineInput2D(Cycles, Xfinal, Xinitial, Yfinal, Yinitial, NumberOfPoints)
%SineInput2D Summary of this function goes here
%   NOTE: this does assume that the number of points is equal
%       in the y direction and the x direction
x = linspace(Xinitial,Xfinal,NumberOfPoints);
y = linspace(Yinitial,Yfinal,NumberOfPoints);
output3 = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial)).*sin(((Cycles*2*pi)/Yfinal)*(y'-Yinitial));
end

function output = SineInput(Cycles,Xfinal, Xinitial, NumberOfPoints)
%SineInput Summary of this function goes here
%   Detailed explanation goes here
x = linspace(Xinitial,Xfinal,NumberOfPoints);
output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
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