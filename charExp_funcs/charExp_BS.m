function z=charExp_BS(u,mu,sigma)
%CHAREXP_BS characteristic exponent of a BS model
% (mu-sigma^2/2)*1i*u-sigma^2
%
%Usage charExp_BS(u,mu,sigma)
%
%mu drift part
%sigma diffusion part


z=(mu-sigma^2/2)*1i*u-sigma^2/2*u.^2;
end