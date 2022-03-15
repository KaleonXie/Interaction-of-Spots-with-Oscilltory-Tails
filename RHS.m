

function rhs = RHS (t, Om_vec)

	global Du Dv Kx Ky K2 Nx Ny Nxy tau k3
    
    u = Om_vec(1:Nxy);
    v=Om_vec(Nxy+1:end);
	uu = reshape(u, Nx, Ny);
    vv= reshape(v, Nx, Ny);
	u_hat = fft2(uu);
    v_hat = fft2(vv);
    w_hat= u_hat./(1+ Dv*K2);
    ww=real( ifft2(w_hat));
    
    
	rhs1 = reshape(real(ifft2(-Du*K2.*u_hat)) + 1.01*uu - uu.^3 - k3*vv - ww -0.1, Nxy, 1);
    rhs2 = reshape( uu - vv , Nxy, 1);
    rhs=[rhs1;rhs2./tau];
end % RHS ()