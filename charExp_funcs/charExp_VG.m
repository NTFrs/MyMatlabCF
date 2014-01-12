function z=charExp_VG(u,theta,sigma,kappa)
%CHAREXP_VG characteristic exponent of VG process with subordinator
%notation
%
%Usage charExp_VG(u,theta,sigma,kappa)
%
%theta drift of brownian motion
%sigma diffusion of brownian motion
%kappa subordinator variance

z=-1/kappa*log(1+u.^2*sigma^2*kappa/2-1i*theta*u*kappa);

end