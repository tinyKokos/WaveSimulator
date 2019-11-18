clear all;
pause on;

initialX = 0;
stepSizeX = pi/32;
finalX = 63/32*pi;

initialT = 0;
stepSizeT = .5;
finalT = 200;

velocity = -1;

h= stepSizeT/(2*stepSizeX);

finalStepT = (finalT-initialT)/stepSizeT;
finalStepX = (finalX-initialX)/stepSizeX+1; 

x = (initialX-stepSizeX):stepSizeX:(finalX+stepSizeX);
A = zeros(length(x),length(x));
u = sin(x + velocity*initialT); %my initial condition has to be the final result of the model
%I think the initial plot of u might have more points than the plots in the
%for loop

figure(1)
plot(x, u)
axis([initialX, finalX, -1, 1])
b = u';
for t = 1:finalStepT
   for xidx = 1:1:finalStepX+2 %this will construct the matrix
       switch xidx
           case 1 %1st element: 2h+1 2nd element: -h last element = -h
               A(1,1) = 2*h+1;
               A(1,2) = -h;
               A(1,finalStepX+2) = -h;
           case (finalStepX+2) % 1st element: -h last element: 2h+1 2nd to last element: -h
               A(finalStepX+2,1) = -h;
               A(finalStepX+2,finalStepX+2) = 2*h+1;
               A(finalStepX+2,finalStepX+1) = -h;
           otherwise
               A(xidx,xidx) = 2*h+1; %u(t,x) term
               A(xidx,xidx+1) = -h; %u(t,x+1) term
               A(xidx,xidx-1) = -h; %u(t,x-1) term
       end
   end
   %disp(A) % 66x66
   %disp(length(u)) %66
   %b = eye(length(u)).*u;
   
   u = linsolve(A,b);
   b = u;
   
   %u(1) = u(finalStepX+1); %this logic should not have to change
   %u(finalStepX+1) = u(2); %this logic should not have to change
   plot((initialX:stepSizeX:finalX), u(2:(length(u)-1)))
   ylim([-1.5, 1.5])
   title(['time = ' num2str(t*stepSizeT+initialT)])
   pause(0.01);
end