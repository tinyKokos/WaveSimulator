% This is my attempt to create a 2D finite difference code for a 2d sine
% wave with the advection equation. Currently the main problems I see is
% that after some time running the loop in time, I end up with an array
% filled with NaN (not a number) values. The output looks very unstable in
% the loop. But I know for sure that the initial function input before the
% loop is correct. 

pause on;

initialX = 0;
stepSizeX = pi/32;
finalX = 63/32*pi;

initialY = initialX;
stepSizeY = stepSizeX;
finalY = finalX;

initialT = 0;
stepSizeT = 0.01;
finalT = 10;

velocity = -1;

finalStepT = (finalT-initialT)/stepSizeT;
finalStepX = (finalX-initialX)/stepSizeX+3; 
finalStepY = (finalY-initialY)/stepSizeY+3;

x = (initialX-2*stepSizeX):stepSizeX:(finalX+2*stepSizeX);
y = (initialY-2*stepSizeY):stepSizeY:(finalY+2*stepSizeY);
u = sin(x + y.' + velocity*initialT);
%I think the initial plot of u might have more points than the plots in the
%for loop
figure(1)
contour(u)

u_new = u;

for t = 1:finalStepT
   x = 3:finalStepX;
   y = 3:finalStepY;
   Residue=(stepSizeT*((3/2*u(x,y)-2*u(x+1,y)+0.5*u(x+2,y))/stepSizeX^2 ... 
                      + (3/2*u(x,y)-2*u(x,y+1)+0.5*u(x,y+2))/stepSizeY^2) ...
                      + u(x,y))-u_new(x,y);
   u_new(x, y.') = u(x, y.') - velocity*Residue./stepSizeX;  
   %u_new(x, y) = u(x, y) - velocity*stepSizeT*
   % (-3/2*u(x)+2*u(x+1)-0.5*u(x+2))/(stepSizeX); 
   u_new(1,1:end) = u_new(finalStepX-1,1:end);
   u_new(2,1:end) = u_new(finalStepX,1:end);
   u_new(finalStepX+2,1:end) = u_new(4,1:end);
   u_new(finalStepX+1,1:end) = u_new(3,1:end);
   
   u_new(1:end,1) = u_new(1:end,finalStepX-1);
   u_new(1:end,2) = u_new(1:end,finalStepX);
   u_new(1:end,finalStepX+2) = u_new(1:end,4);
   u_new(1:end,finalStepX+1) = u_new(1:end,3);
   
   contour(u_new)
   %title(['time = ' num2str(t*stepSizeT+initialT)])
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
% 
% Forward:      (x+1) - (x)
%               negative v
%               loop initial:Final