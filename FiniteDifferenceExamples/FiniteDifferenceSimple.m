%Utt = Vx^2Uxx
N = 51;  
DX=0.1; % step size
Nx=5; 
X=0:DX:Nx; %0 to 5 step 0.1 ->  50 points
Vx = 5;

%boundary conditions
U(1:N,1:N) = 0;
U(1,1:N) = 100; 
U(N,1:N) = 0;  
U(1:N,1) = 0;  
U(1:N,N) = 0;  

%Initial Condition

%end

DTx = DX^2/(2*Vx); % time step
M=2000; % maximum number of allowed iteration

fram=0;
Ncount=0;
loop=1;
while loop==1;
    U_old = U;
    for i = 2:N-1
            Residue=(DTx*((U_old(i+1)-2*U_old(i)+U_old(i-1))/DX^2 + U_old(i)))-U(i);
            U(i)=U(i)+Residue;
        
    end
    Ncount=Ncount+1;
    if (mod(Ncount,50)==0) % displays movie frame every 50 time steps
        fram=fram+1;
        plot(U);
        axis([1 N])
        h=gca; 
        get(h,'FontSize') 
        set(h,'FontSize',12)
        colorbar('location','eastoutside','fontsize',12);
        xlabel('X','fontSize',12);
        title('Finite Difference Simulation','fontsize',12);
        fh = figure(1);
        set(fh, 'color', 'white'); 
        F=getframe;
    end
    if(Ncount>M)
            loop=0; %make the program exit the loop
            disp(['Finished ',num2str(M),...
            'time steps'])
    end
end
%display the output
movie(F,fram,1)