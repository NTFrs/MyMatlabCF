function [X,t]=MCL_NIG(theta,sigma,kappa,T,Nsteps,Nsims)
%MCL_NIG Generates paths of a NIG levy process
%
%usage MCL_NIG(S0,sigma,eta,kappa,T,Nsteps,Nsims)
%
%S0 asset value at time 0
%theta brownian drift
%sigma brownian volatility
%kappa subordinator volatility
%T time horizon
%Nsteps Number of time steps
%Nsims number of simulations

t=linspace(0,T,Nsteps+1); dt=t(2)-t(1);

lambda=dt^2/kappa; mu=dt;

Y=randn(Nsims,Nsteps).^2;
U=rand(Nsims,Nsteps);
dS=mu+mu^2*Y/(2*lambda)-mu/(2*lambda)*sqrt(4*mu*lambda*Y+mu^2*Y.^2);

for i=1:Nsims
    for j=1:Nsteps
        if U(i,j)>mu/(mu+dS(i,j))
            dS(i,j)=mu^2/dS(i,j);
        end
    end
end



X=zeros(Nsims,Nsteps+1);
gauss=randn(Nsims,Nsteps);

for j=1:Nsteps
    X(:,j+1)=X(:,j)+theta*dS(:,j)+sigma*gauss(:,j).*sqrt(dS(:,j));
end

end