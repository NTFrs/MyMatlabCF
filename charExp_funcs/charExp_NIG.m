function z=charExp_NIG(u,theta,sigma,kappa)
%CHAREXP_NIG characteristic exponent of NIG process with subordinator
%notation
%
%Usage charExp_NIG(u,theta,sigma,kappa)
%
%theta drift of brownian motion
%sigma diffusion of brownian motion
%kappa subordinator variance

z=1/kappa-1/kappa*sqrt(1+u.^2*sigma^2*kappa-2i*theta*u*kappa);

end
