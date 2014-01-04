function Price=EuroBS_FD_2(S0,K,r,d,sigma,T,M,N,FlagCP,theta)
%EUROBS_FD_D Computes price of an european option using DF
%
%Computes price of a vanilla option  in B&S framework using a FD method. By
%default it uses Crank-Nicolson algorithm for time discretization, thus
%being unconditionnaly stable.
%
%Output:
%Price at time 0 and S value S0
%
%Inputs:
%S0: Initial stock value
%K: Strike
%r: risk free constant interest rate
%d: dividend rate
%sigma: volatility
%T: maturity
%M: Number of time intervals
%N: Number of "space" intervals
%FlagCP: 1 for Call, -1 (or anything else) for Put
%theta: (OPTIONAL) theta coefficient:
%        theta=0 Implicit Euler
%        theta=1 Explicit Euler

if (nargin<10)
    theta=0.5;
end

ctheta=1-theta;

Smin=S0*exp( (r-sigma^2/2)*T-sigma*sqrt(T)*6)
Smax=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*6)

S=linspace(Smin,Smax,N+1)'; dS=S(2)-S(1);
t=linspace(0,T,M+1)'; dt=t(2)-t(1);

nodes=S(2:end-1);

lambda=dt/dS^2;
if theta > 0.5
if lambda>0.5
    disp('May have stability issues')
end
end

if FlagCP==1
    r1=@(j) 0;
    r2=@(j) Smax-K*exp(-r*(T-t(j)));
    V=max(nodes-K,0);
else
    r1=@(j) K*exp(-r*(T-t(j)))-Smin;
    r2=@(j) 0;
    V=max(K-nodes,0);
end

%Implicit part matrix
underdiag=-nodes*(r-d)/(2*dS)+0.5*sigma^2*nodes.^2/dS^2;
BCund=-underdiag(1);

diag=-sigma^2*nodes.^2/dS^2-r;

upperdiag=nodes*(r-d)/(2*dS)+0.5*sigma^2*nodes.^2/dS^2;
BCupp=-upperdiag(end);

tridiag=[ctheta*underdiag, -1/dt+ctheta*diag, ctheta*upperdiag];
A=spdiags(tridiag,1:-1:-1,N-1,N-1)';

tridiag=[-theta*underdiag,-1/dt-theta*diag,-theta*upperdiag];
B=spdiags(tridiag,1:-1:-1,N-1,N-1)';
clear tridiag underdiag upperdiag diag

BC=sparse(N-1,1);

for j=M:-1:1
    BC(1)=ctheta*BCund*r1(j)+theta*BCund*r1(j+1);
    BC(end)=ctheta*BCupp*r2(j)+theta*BCupp*r2(j+1);
    V=A\(B*V+BC);
end

V=[r1(1);V;r2(1)];

plot(S,V);
title('Price of option as function of S');
xlabel('Underlying S');
ylabel('Option price');

Price=interp1(S,V,S0,'spline');
end
