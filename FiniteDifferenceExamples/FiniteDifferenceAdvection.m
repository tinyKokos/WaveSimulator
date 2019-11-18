pause on;

initialX = 0;
stepSizeX = pi/32;
finalX = 63/32*pi;

initialT = 0;
stepSizeT = 0.01;
finalT = 10;

velocity = -1;

finalStepT = (finalT-initialT)/stepSizeT;
finalStepX = (finalX-initialX)/stepSizeX+3; 

x = (initialX-2*stepSizeX):stepSizeX:(finalX+2*stepSizeX);
u = sin(x + velocity*initialT);
%I think the initial plot of u might have more points than the plots in the
%for loop

figure(1)
plot(initialX:stepSizeX, u)

u_new = u;

for t = 1:finalStepT
   x = 3:finalStepX;
   u_new(x) = u(x) - velocity*stepSizeT*(-3/2*u(x)+2*u(x+1)-0.5*u(x+2))/(stepSizeX); 
   u_new(1) = u_new(finalStepX-1);
   u_new(2) = u_new(finalStepX);
   u_new(finalStepX+2) = u_new(4);
   u_new(finalStepX+1) = u_new(3);
   
   plot(x, u_new(3:finalStepX))
   title(['time = ' num2str(t*stepSizeT+initialT)])
   pause(0.01);
   
   u = u_new;
end

% notes on forward versus backward loop
%
% Backward:     (x) - (x-1)
%               positive v
%               loop Final:-1:initial
% 
% Forward:      (x+1) - (x)
%               negative v
%               loop initial:Final