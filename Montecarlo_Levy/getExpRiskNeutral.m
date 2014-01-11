function S=getExpRiskNeutral(S0,r,T,X,char)
%GETEXPRISKNEUTRAL returns the underlyng risk neutral dynamic
%
%Usage getExpRiskNeutral(S0,X,char)
%
%Computes the values of an underlyng S for some model where X are it's
%realizations, whose caracteristic function is specified by char
%
%S0 starting value
%r risk free interest rate
%X vector of NxM values with N number of observed or simulated and M-1 time
%   steps
%char function handle with the characteristic function of X
%
%Output exp((r-char(-i))*t+X(t)) ie risk neutral exponential dynamic

Nsteps=size(X,2)-1;
Nsims=size(X,1);
t=linspace(0,T,Nsteps+1);

S=S0*exp((r-char(-1i))*ones(Nsims,1)*t+X);


end