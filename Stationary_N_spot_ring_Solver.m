%%% -------------------------------------------------- %%%
%%% 2D Gas Discharge Model pseudo-spectral solver            %%%
%%% -------------------------------------------------- %%%


clear, close all
format longE

% ------------------------------------------------------------ %

% declaration of global variables:
global x y dy Du Dv Kx Ky K2 Lx Ly Nx Ny Nxy X Y tau k3 spots_num
% ------------------------------------------------------------ %
%%% Physical parameters:
Du = 1.1e-4;	% 
Dv=9.64*1e-4;
k3= 0.3;
tau =0.1;%1/k3+0.005;

% ------------------------------------------------------------ %

%%% Domain definition:
Lx = 1;    % Domain half-length in x-direction
Ly = 1;    % Domain half-length in y-direction

%%% Numerical parameters:
Nx = 256;  % number of Fourier modes in discrete solution x-dir
Ny = 256;	% number of Fourier modes in discrete solution y-dir
Nxy = Nx*Ny;

dx = 2*Lx/Nx;   		% distance between two physical points
x = (1-Nx/2:Nx/2)'*dx;  % physical space discretization

dy = 2*Ly/Ny;   		% distance between two physical points
y = (1-Ny/2:Ny/2)'*dy;  % physical space discretization

[X,Y] = meshgrid(x,y);	% 2D composed grid

% ------------------------------------------------------------ %

% vectors of wavenumbers in the transformed space:
kx = [0:Nx/2 1-Nx/2:-1]'*pi/Lx;
ky = [0:Ny/2 1-Ny/2:-1]'*pi/Ly;

% ------------------------------------------------------------ %

%%% Some operators arising in GasDischarge equations:
[Kx, Ky] = meshgrid(kx,ky);
K2 = Kx.^2 + Ky.^2;     % to compute the Laplace operator



% ------------------------------------------------------------ %

fftw('planner', 'hybrid');

% ------------------------------------------------------------ %

%%% Time-stepping parameters:
t = 0.0;           	% the discrete time variable
Tf = 100000;          	% final simulation time
ds = 1;			% write time of the results
tol= 1e-4;   % Set the stop criterion
ops = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);



% ------------------------------------------------------------ %
%%% Set initial conditions.  
spots_num=2; % number of spots 
     thetas=linspace(0,2*pi,spots_num+1);
     init_r=0.08; % initial radius of the ring for the first binding radius.    For the second binding radius, we chooose 0.16;
     centers_x=  init_r/sin(pi/spots_num).*cos(thetas); 
     centers_y=  init_r/sin(pi/spots_num).*sin(thetas);
     u=-0.27+0*X;
     v=-0.27+0*X;
     for k=1:spots_num % place N spots on a ring
        u =u+2* exp(-((X-centers_x(k)).^2+(Y-centers_y(k)).^2)*500)+k*0.001;

        v =v+ 2*exp(-((X-centers_x(k)).^2+(Y-centers_y(k)).^2)*500);

     end
   
    initial_u=reshape(u, Nxy, 1);
    initial_v=reshape(v, Nxy, 1);
    Om_vec=[initial_u;initial_v];


%% Plot initial condition
FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 1200, 500]);

frame = 0;
xcs=[];  % record spots' locations.
ycs=[];
ts=[];
 [x_c,y_c]=FindSpotLocation(X,Y,u,spots_num); % Get spots' locations.
 xcs=[xcs;x_c];
 ycs=[ycs;y_c];
 ts=[ts,t];

 plot_Contour(X,Y,u,t,x_c,y_c);% we plot the initial condition
 pause(0.1);

%% Run simulation
 
while (t < Tf) % main loop in time
    
	[~, sol] = ode113(@RHS, [0:ds:ds], Om_vec, ops); % Solve discrete PDE system
    
    err=max(abs(sol(end,:)-sol(1,:))); % Set tolerance.
    
	Om_vec = sol(end,:);
    
    u=Om_vec(1:Nxy);
    v=Om_vec(Nxy+1:end);
	
    Om_u = reshape(u, Nx, Ny);
    Om_v= reshape(v, Nx, Ny);

	t = t + ds; frame = frame + 1;
    
    [x_c,y_c]=FindSpotLocation(X,Y,Om_u,spots_num);
    

    subplot(1,2,1)
    contourf(X,Y,Om_u-Om_u(end,end),'LineColor','none');
    title ([' t = ',num2str(t,'%4.4f'), ' $\tau$ = ',num2str(tau),' error= ', num2str(err) ], 'interpreter', 'latex', 'fontsize', 12);
    colorbar;
    grid on
    grid minor
	shading interp;
    hold on;
    for k=1:spots_num
        plot(x_c(k),y_c(k),'ro','markersize',12)
    end
    axis equal
    hold off
    
    subplot(1,2,2)
    xcs=[xcs;x_c];
    ycs=[ycs;y_c];
    ts=[ts,t];
    
    plot(x,Om_u(Nx/2,:)','-') 
    ylim([-0.5,0.5])
    grid minor


    drawnow;
    
    if err<tol;
        break;
    end
    
    
    
end 
 data_name= sprintf("%dSpots.mat",spots_num );
 save( data_name,'X','Y','Om_u','Om_v'); % save stationary N-spot ring profile.




