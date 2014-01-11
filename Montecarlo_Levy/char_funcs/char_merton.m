function z=char_merton(u,gamma,sigma,lambda,jump_mu,jump_delta)

z=-sigma^2*u.^2/2+1i*gamma*u+lambda*(exp(-jump_delta*u.^2/2+1i*jump_mu*u)-1);

end
