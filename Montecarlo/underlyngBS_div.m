function [S,S_AV]=underlyngBS(S0,r,d,sigma,T,Nsim,Nsteps)
%UNDERLYNGBS Simulation of the dynamic of a BS underlyng
%
%Generates Nsim simulations with Nsteps time steps. If the second argument
%is present, generates also the antitethic variable corresponding to S
%
%INPUTS:
%S0 spot price
%r risk free constant interest rate
%sigma volatility
%T Time horizon
%Nsim number of simulations
%Nsteps number of time steps
%
%OUTPUTS:
%S simulated underlying
%S_AV antithetic simulated underlying

dt=T/(Nsteps);
temp=randn(Nsim,Nsteps);
S=zeros(Nsim,Nsteps+1); S(:,1)=S0;
S_AV=zeros(Nsim,Nsteps+1); S_AV(:,1)=S0;

if nargout==2
    for i=1:Nsteps
        S(:,i+1)=S(:,i).*exp( (r-d-sigma^2/2)*dt+...
            sigma*sqrt(dt)*temp(:,i) );
        S_AV(:,i+1)=S_AV(:,i).*exp( (r-d-sigma^2/2)*dt-...
            sigma*sqrt(dt)*temp(:,i) );
    end
else
    
    for i=1:Nsteps
        S(:,i+1)=S(:,i).*exp( (r-d-sigma^2/2)*dt+...
            sigma*sqrt(dt)*temp(:,i) );
    end
end

end