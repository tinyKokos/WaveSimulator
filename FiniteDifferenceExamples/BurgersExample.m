% Burgers program
% this program computes numerical solutions for viscid and inviscid Burgers
% equation: PDF founded @ http://www.bcamath.org/projects/NUMERIWAVES/Burgers_Equation_M_Landajuela.pdf

% It needs the functions df.m, f.m, nf.m, and uinit.m

clear all;

% Selection of equation and method.
type_burger = menu('Choose the equation:', ...
    'Inviscid Burgers Equation','Viscid Burgers equation');

if (type_burger == 1)
    method = menu('Choose a numerical method:', ...
        'Up-Wind Nonconservative','Up-wind conservative','Lax-Friedrichs',...
        'Lars-Wendroff','MacCormack','Godunov');
else
    method = menu('Choose a numerical method:','Parabolic Method');
end

ictype = menu('Choose the initial condition type:',...
    'Piecewise constant (shock)','Piecewise constant (expansion)',...
    'Gaussian','Piecewise continuous');

% Selection of numerical parameters (time step, grid spacing, etc...)
xend = 2; % x-axis size
tend = 2; % t-axis size
N = input('Enter the number of grid points: ');
dx = xend/N; % Grid spacing
dt = input('Enter time step dt: ');

x = 0:dx:xend;
nt = floor(tend/dt);
dt = tend/nt;

% Set up the initial solution values
u0 = uinit(x,ictype); % Call to the function "uinit"
u = u0;
unew = 0*u;

% Implementation of the numerical methods
if (type_burger == 1)
    for i = 1:nt,
        switch method
            case 1 % Up-wind nonconservative
                unew(2:end) = u(2:end) - dt/dx .* (u(2:end) - u(1:end-1));
                unew(1) = u(1); % u(3:end): subvector de u desde 3 hasta el final
            case 2 % Up-wind conservative
                unew(2:end) = u(2:end) - dt/dx * (f(u(2:end)) - f(u(1:end-1)));
                unew(1) = u(1);
            case 3 % Lax-Friedrichs
                unew(2:end-1) = 0.5*(u(3:end) + u(1:end-2)) - 0.5*dt/dx * ...
                    (f(u(3:end)) - f(u(1:end-2)));
                unew(1) = u(1);
                unew(end) = u(end);
            case 4 % Lax-Wendroff
                unew(2:end-1) = u(2:end-1)  ...
                    - 0.5*dt/dx * (f(u(3:end)) - f(u(1:end-2))) ...
                    + 0.5*(dt/dx)^2 * ...
                    (df(0.5*(u(3:end) + u(2:end-1))) .* (f(u(3:end)) - f(u(2:end-1))) - ...
                    df(0.5*(u(2:end-1) + u(1:end-2))) .* (f(u(2:end-1)) - f(u(1:end-2))) );
                unew(1) = u(1);
                unew(end) = u(end);
            case 5 % MacCormack
                us = u(1:end-1) - dt/dx * (f(u(2:end)) - f(u(1:end-1)));
                unew(2:end-1) = 0.5*(u(2:end1) + us(2:end)) - ...
                    0.5*dt/dx * (f(us(2:end)) - f(us(1:end-1)));
                unew(1) = u(1);
                unew(end) = u(end);
            case 6 % Godunov
                unew(2:end-1) = u(2:end-1) - dt/dx*(nf(u(2:end-1),u(3:end)) - nf(u(1:end-2),u(2:end-1)));
                unew(1) = u(1);
                unew(end) = u(end);
                
        end
        
        u = unew;
        U(i,:) = u(:);
    end
    
else
    for i = 1:nt,
        switch method
            case 1 % Parabolic method
                D = 0.01;
                fminus = 0.5*( f(u(2:end-1)) + f(u(1:end-2)) );
                fplus = 0.5*( f(u(2:end-1)) + f(u(3:end)) );
                unew(2:end-1) = u(2:end-1) + dt*(D*(u(3:end)-2*u(2:end-1)+u(1:end-2))/(dx)^2 ...
                    - (fplus - fminus)/dx );
                unew(1) = u(1);
                unew(end) = u(end);
                
        end
        u = unew;
        U(1,:) = u(:);
    end
end

U = [u0:U];
T = 0:dt:tend;

% Plot of the solutions
figure(1)
surf(x,T,U)
shading interp
xlabel('x'), ylabel('t'), zlabel('u(x,t)');
grid on
colormap('Gray');