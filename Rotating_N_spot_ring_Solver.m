Stationary_N_spot_ring_Solver;
% Get Stationary N-spot ring first; If you have already got the stationary steady state, ignore this line.

tau=1/k3+ 0.01; % increase the bifurcation parameter beyond the threshold 1/k_3;

u_hat=fft2(Om_u);
dx_u=real(ifft2(1i*Kx.*u_hat));
dy_u=real(ifft2(1i*Ky.*u_hat));


Om_v=Om_v-0.005*Y.*(dx_u)+0.005.*X.*(dy_u); % add a rotating pertubation to the stationary state.
v=reshape(Om_v,1,Nxy);
Om_vec=[u,v];
close all;

ds=10;
imag_index=0;
t=0;

FigHandle = figure(2);
set(FigHandle, 'Position', [100, 100, 600, 600]);
Tf=100000;
tt=4000;

%% Run simulation
while (t < Tf) % main loop in time
    
     % if t>tt
     %     tt=tt+4000;  %slowly increasing tau if necessary.
     %     tau=tau+0.01;  
     % end
    
	[~, sol] = ode113(@RHS, [0:ds:ds], Om_vec, ops);
    
    err=max(sol(end,:)-sol(1,:));
    
	Om_vec = sol(end,:);
    
    u=Om_vec(1:Nxy);
    v=Om_vec(Nxy+1:end);
	
    Om_u = reshape(u, Nx, Ny);
    Om_v= reshape(v, Nx, Ny);

	t = t + ds; frame = frame + 1;
    
    [x_c,y_c]=FindSpotLocation(X,Y,Om_u,spots_num);
    
    radius=mean( sqrt(x_c.^2+y_c.^2) );
    contourf(X,Y,Om_u-Om_u(end,end),'LineColor','none');
    title ([' t = ',num2str(t,'%4.3f'), ' $\tau$ = ',num2str(tau,'%4.3f'),' $r_0$ = ', num2str(radius) ], 'interpreter', 'latex', 'fontsize', 12);
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
    imag_name=sprintf("Rotating%dSpot%d.png",spots_num,imag_index)
    saveas(gcf,imag_name)
    imag_index=imag_index+1;
    xcs=[xcs;x_c];
    ycs=[ycs;y_c];
    ts=[ts,t];



    drawnow;
end 

AjustTrajectoriesm; %
save('Rotating_trajecotories.mat','ts','newxc','newyc')
figure(3)
plot(ts,newxc); % plot x-cooridiantes of N spots w.r.t time .
 

