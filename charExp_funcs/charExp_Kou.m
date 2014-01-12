function z=charExp_Kou(u,gamma,sigma,lambda,lambda_u,lambda_d,prob)
%CHAREXP_KOU characteristic exponent of Kou process with subordinator
%notation
%
%Usage charExp_Kou(u,gamma,sigma,lambda,lambda_u,lambda_d,prob)
%
%gamma drift part
%sigma diffusion part
%lambda jump intensity
%lambda_u positive jumps parameter
%lambda_d negative jumps parameter
%prob probability of a positive jump

z=-sigma^2*u.^2/2+1i*gamma*u+1i*u*lambda.*...
    (prob./(lambda_u-1i*u)-(1-prob)./(lambda_d+1i*u));

end
