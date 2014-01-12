function z=charExp_Merton(u,gamma,sigma,lambda,jump_mu,jump_delta)
%CHAREXP_MERTON characteristic exponent of Merton process with subordinator
%notation
%
%Usage charExp_Merton(u,gamma,sigma,lambda,jump_mu,jump_delta)
%
%gamma drift part
%sigma diffusion part
%lambda jump intensity
%jump_mu jumps mean
%jump_delta jumps variance


z=-sigma^2*u.^2/2+1i*gamma*u+lambda*(exp(-jump_delta*u.^2/2+1i*jump_mu*u)-1);

end
