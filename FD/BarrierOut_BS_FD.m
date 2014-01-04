function Price=BarrierOut_BS_FD(S0,K,r,d,sigma,T,M,N,FlagCP,varargin)
%EUROBS_FD_D Computes price of a knock out barrier option using DF
%
%Computes price of a knock out barrier option  in B&S framework using a FD
%method. By default it uses Crank-Nicolson algorithm for time
%discretization, thus being unconditionnaly stable.
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
%
%Optional Inputs:
% ..,'Theta',value,... theta coefficient:
%                       theta=0 Implicit Euler
%                       theta=1 Explicit Euler
%...,'Low',Number,... to specify lower Barrier
%...,'Upp',Number,... to specify Upper Barrier

indexU=find(strcmp('Upp',varargin));
indexL=find(strcmp('Low',varargin));
indexT=find(strcmp('Theta',varargin));

if(isempty(indexT)==0)
    theta=varargin{indexT+1};
else
    theta=0.5;
end

ctheta=1-theta;

if(isempty(indexU)==0)
    U=varargin{indexU+1};
end

if(isempty(indexL)==0)
    L=varargin{indexL+1};
end

if(exist('L','var'))
    Smin=L;
    r1=@(j) 0;
else
    Smin=S0*exp( (r-sigma^2/2)*T-sigma*sqrt(T)*6);
end

if(exist('U','var'))
    Smax=U;
    r2=@(j) 0;
else
    Smax=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*6);
end

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
    if(~exist('U','var'))
    r2=@(j) Smax-K*exp(-r*(T-t(j)));
    end
else
    if(~exist('L','var'))
    r1=@(j) K*exp(-r*(T-t(j)))-Smin;
    end
    r2=@(j) 0;
end

V=max(FlagCP*(nodes-K),0);
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

figure
plot(S,V);
title('Price of option as function of S');
xlabel('Underlying S');
ylabel('Option price');

figure
plot(nodes, (V(3:end)-V(1:end-2))/(2*dS))
title('Delta');
figure
plot(nodes, (V(3:end)-2*V(2:end-1)+V(1:end-2))/(dS^2))
title('Gamma');

Price=interp1(S,V,S0,'spline');
end