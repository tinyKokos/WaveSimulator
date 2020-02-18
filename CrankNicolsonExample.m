%% CrankNicolsonExample.m
% This function is written to simulate a standing wave in 1 dimensional
% space with the implicit Finite Difference Method. This is solved in
% explicit Finite Difference. This will be my 2nd order implicit file
% solving a combination of advection and the diffusion equations.

%% Declared Input Values
close all; 

% Spacial Inputs
initialX = 0;
finalX = 3;
Npoints = 50; %number of points between the initial and final X
Ncycles = 1; %number of full sine cycles

%set inputs for the time
initialTime = 0;
finalTime = 1;
NtimePoints = 100000; 

% Group Velocity of the Wave
PropagationSpeed = 300;
% Dissapation Constant
D = 10;

%filename of gif output of the program
filename = 'output.gif';

%delay value between each frame of the gif
FrameDelay = 0.1;

%% Calculate Constants for Later

%find delta constants
tDelta = (finalTime-initialTime)/NtimePoints;
xDelta = (finalX-initialX)/Npoints;

%this is for efficiency because it is used so often
M = Npoints - 1; %end array location
N = M - 1; %matrix and vector/subvector size

%create space vector
x = linspace(initialX,finalX,Npoints);

%create solving matrix
A = zeros(N, N);
%create current function vector
b = zeros(1,N);

%Create constants for Update Equation
constant1 = (PropagationSpeed*tDelta)/(2*xDelta^2);
constant2 = 2*constant1;
constant3 = 1+constant2;
constant4 = 1-constant2;
%% Create Initial Input
% Also create arrays for the function which we use as 
% pastFunc ->   f(x,t-1)
% func ->       f(x, t)
% nextFunc ->   f(x,t+1)
func = SineInput(Ncycles, finalX, initialX, Npoints);
%nextFunc = zeros(length(func));
%pastFunc = zeros(length(func));
%% CFL Checking
% I just copied this from the Fixed_Force_Fixed_1D.m file 
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
                oldTimePoints = NtimePoints;
                NtimePoints = ceil((finalTime-initialTime)/tDelta);
                EstimateRunTime = NtimePoints*(FrameDelay+0.01);
                fprintf('This will adjust system from %d to %d time points. Estimated run time will be %.2f seconds\n'...
                    , oldTimePoints, NtimePoints, EstimateRunTime);
                prompt = 'Would this be acceptable? (y/n): ';
                Choice = input(prompt, 's');
                if strcmp(Choice, 'y') == 1
                    %empty case to continue the program
                elseif strcmp(Choice, 'n') == 1
                    return
                else
                    printf('Input is not recognized. Exiting execution.\n');
                    return
                end
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
%% Display & Update Loop in Time
% create display
h = figure;
for t = 1:NtimePoints
    %plot the new output
    plot(x(2:end-1), func(2:end-1));
    %This would have to be adjusted if amplitude is adjusted
    ylim([-2 2]);    
    %% Create present time vector (b vector)
    for r = 2:N
        b(r) = constant1*func(r+1)+constant4*func(r)+constant1*func(r-1);
    end
    %% Create Matrix
    for n = 1:N
        switch n
            case 1
                A(1,1) = constant3;
                A(1,2) = -1*constant1;
                A(1,N) = -1*constant1;
            case N
                A(N,N-1) = -1*constant1;
                A(N,N) = constant3;
                A(N,1) = -1*constant1;
                break;
            otherwise
                A(n,n-1) = -1*constant1;
                A(n,n) = constant3;
                A(n,n+1) = -1*constant1;
        end
    end
 
    %% Solve for Next Time Step
    
    opts.TRANSA = false;
    opts.RECT = true;
    
    a = transpose(b);
    c = linsolve(A,a,opts);
    d = transpose(c);
    func(2:end-1) = d;
    %% Force Function to be 0 at Boundaries
    func(1) = 0;
    func(end) = 0;

    pause(FrameDelay);
end
%% SineInput(Cycles, Xfinal, Xinitial, NumberOfPoints)
%   Cycles - How many full cycles the output array must complete in the 
%   given amount of space
%   Xfinal - the value of the final point in space of the array
%   Xinitial - the value of the initial point in the space of the array
%   NumberOfPoints - how many discrete points the output array must have
%    
%   Creates a input for the Finite Difference update equation to use as the
%   initial condition of the function. This input is a simple sine wave
%   that is fixed to be 0 at the boundaries of the given space
function output = SineInput(Cycles, Xfinal, Xinitial, NumberOfPoints)

%Check for valid arguments
arguments
   Cycles (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
   Xfinal (1,:) {mustBeNumeric, mustBeFinite}
   Xinitial (1,:) {mustBeNumeric, mustBeFinite}
   NumberOfPoints (1,:) {mustBeNumeric, mustBeFinite, mustBePositive}
end

    x = linspace(Xinitial,Xfinal,NumberOfPoints);
    output = sin(((Cycles*2*pi)/Xfinal)*(x-Xinitial));
end