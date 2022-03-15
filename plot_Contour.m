function plot_Contour(X,Y,u,t,x_c,y_c)
global x y tau spots_num Nx Ny


subplot(1,2,1)
    contourf(X,Y,u,'LineColor','none');
    title ([' t = ',num2str(t,'%4.2f'), ' $\tau$ = ',num2str(tau,'%4.2f')], 'interpreter', 'latex', 'fontsize', 12);
    colorbar;
    grid on
    grid minor
	shading interp;
    hold on;
    for i=1:spots_num
        plot(x_c(i),y_c(i),'ro','markersize',12)
    end
    axis equal
    hold off
    
    subplot(1,2,2)
    plot(x,u(Nx/2,:),'o-') 
    ylim([-0.5,0.5])