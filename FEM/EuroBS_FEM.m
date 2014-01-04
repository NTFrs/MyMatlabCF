function Price=EuroBS_FEM(S0,K,r,d,sigma,T,M,N,FlagCP,theta)
%EUROBS_FEM Computes price of an european option using FEM
%
%Computes price of a vanilla option  in B&S framework using a FEM method.
%By default it uses Crank-Nicolson algorithm for time discretization, thus
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

Smin=S0*exp( (r-sigma^2/2)*T-sigma*sqrt(T)*6);
Smax=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*6);

x=linspace(log(Smin),log(Smax),N+1)'; dx=x(2)-x(1);
t=linspace(0,T,M+1)'; dt=t(2)-t(1);
nodes=x(2:end-1);

lambda=dt/dx^2;
if theta > 0.5
if lambda>0.5
    disp('May have stability issues')
end
end

if FlagCP==1
    r1=@(j) 0;
    r2=@(j) Smax-K*exp(-r*(T-t(j)));
    V=max(exp(nodes)-K,0);
else
    r1=@(j) K*exp(-r*(T-t(j)))-Smin;
    r2=@(j) 0;
    V=max(K-exp(nodes),0);
end

tridiag=[ones(N-1,1)*dx/6,ones(N-1,1)*2*dx/3,ones(N-1,1)*dx/6];
B=spdiags(tridiag,-1:1,N-1,N-1);
tridiag=[ones(N-1,1)*(-1)/dx,ones(N-1,1)*2/dx,ones(N-1,1)*(-1)/dx];
A=spdiags(tridiag,-1:1,N-1,N-1);
clear tridiag

bidiag=[ones(N-1,1)*(-0.5),ones(N-1,1)*0.5];
C=spdiags(bidiag,[-1 1],N-1,N-1);

M1=B*(1/dt+r*ctheta)-(r-sigma^2/2)*ctheta*C+ctheta*sigma^2/2*A;
M2=B*(1/dt-r*theta)+(r-sigma^2/2)*theta*C-theta*sigma^2/2*A;

BC1_0=M1(2,1); BC2_0=M2(2,1); BC1_N=M1(1,2); BC2_N=M2(1,2);

BC=sparse(N-1,1);

for k=M:-1:1
    BC(1)=-BC1_0*r1(k)+BC2_0*r1(k+1);
    BC(end)=-BC1_N*r2(k)+BC2_N*r2(k+1);
    V=M1\(M2*V+BC);
end

Price=interp1(nodes,V,log(S0),'spline');
V=[r1(1);V;r2(1)];
close all
S=exp(x);
figure
plot(S,V)
title('Prezzo')

error=abs(Price-Exact_Euro_BS(S0,K,r,0,T,sigma,FlagCP))

end