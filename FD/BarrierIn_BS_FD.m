function Price=BarrierIn_BS_FD(S0,K,r,d,sigma,T,M,N,FlagCP,varargin)
%BARRIERIN_BS_FD Computes price of a knock in barrier option using DF
%
%Computes price of a knock in barrier option  in B&S framework using a FD
%method to compute the coresponding out option price. By default it uses
%Crank-Nicolson algorithm for time discretization, thus being
%unconditionnaly stable.
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

Price=Exact_Euro_BS(S0,K,r,0,T,sigma,FlagCP);

if(isempty(indexU)==0)
    U=varargin{indexU+1};
    if(isempty(indexL)==0)
        L=varargin{indexL+1};
        Price=Price-BarrierOut_BS_FD(S0,K,r,d,sigma,T,M,N,FlagCP,'Upp',U,'Low',L,'Theta',theta);
    else
        Price=Price-BarrierOut_BS_FD(S0,K,r,d,sigma,T,M,N,FlagCP,'Upp',U,'Theta',theta);
    end
elseif(isempty(indexL)==0)
        L=varargin{indexL+1};
        Price=Price-BarrierOut_BS_FD(S0,K,r,d,sigma,T,M,N,FlagCP,'Upp',U,'Low',L,'Theta',theta);
else
    disp('No Barrier specified');
    Price=0;
end

end