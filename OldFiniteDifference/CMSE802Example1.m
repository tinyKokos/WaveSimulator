
xmin = 0;
xmax = 10;
nx = 500;

dx = (xmax-xmin)/nx;
x = linspace(xmin,xmax,nx);

tmin = 0;
tmax = 10;
nt = 1000000;

dt = (tmax-tmin)/nt;
times = linspace(tmin,tmax,nt);

yi = exp(-((x)-5).^2);

dydt = zeros(length(x)); 
dydt2 = zeros(length(x));

gamma = 1;

for times = tmin:dt:tmax
    for x = xmin:dx:xmax
                dydt2(1) = 0;
                dydt2(end) = 0;
                dydt2 = gamma*(yi(+1) + y(x-1) - 2*y(x))/(dx)^2;
                
        end
    end
end