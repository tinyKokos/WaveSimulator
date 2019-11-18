%should I determine the movement of the wave from velocity, stepsize of T,
%or stepsize of X (I can only pick 2 of the 3)

%Define Spacial Boundary
initialX = 0;
stepSizeX = pi/32;
finalX = 63/32*pi;

%Define Time Boundary
initialT = 0;
stepSizeT = 0.01;
finalT = 5;

%Initial state points
velocity = 5; %this will only be used to determine sign of the wave
x = (initialX-stepSizeX):stepSizeX:(finalX+stepSizeX); 
%I need to create 2 ghost cells with 1 before the initial 
% state and 1 step beyond the final state to provide enough
% data for the central difference method to compute

u = sin(x + velocity*initialT);

%Plot initial state of function
figure(1)
plot(initialX:stepSizeX:finalX, u(2:end-1))
axis([initialX,finalX,-1.5,1.5])
title(['time = ' num2str(initialT)])

FinalStepX = (finalX - initialX)/stepSizeX+1;
% this is to fix the for loop on the correct array subset

space = linspace(initialX,finalX,(FinalStepX-2));

%Create u_new array
u_new = u;

for t = initialT:stepSizeT:(finalT-stepSizeT)
   for x = 2:FinalStepX
       RightFlux = velocity*(0.5*u(x+1) - u(x));
       LeftFlux = velocity*(u(x) - 1.5*u(x-1));
       u_new(x) = u(x) + (stepSizeT/(stepSizeX))*(RightFlux - LeftFlux);
   end
   % I think this is for the boundary conditions
   % after looking at my new results I think the program is not
   % properly shifting the wave values from the right to the left
   % boundary
   u(1) = u(end-1); % This is to shift
   u(end) = u(2);   % the ghost cell values
   
   perfect = sin((initialX:stepSizeX:finalX) + velocity*(t+stepSizeT));
   figure(1);
   plot(initialX:stepSizeX:finalX, u_new(2:end-1))%, initialX:stepSizeX:finalX, perfect)
   axis([initialX,finalX,-1.5,1.5])
   title(['time = ' num2str(t + stepSizeT)])
   pause(stepSizeT);
   
   if length(u_new) ~= length(u)
       print('u_new does not match u size')
   else 
        u = u_new;
   end
  
   %numeric_peak = find(u==max(u));
   %analytical_peak = find(perfect==max(perfect));
   
end

% for Attenuation: set RightFlux = u(x+1) - u(x) and LeftFlux = u(x) - u(x-1)