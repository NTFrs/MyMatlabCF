function [X,V]=Heston_Sim(x0,v0,r,T,kappa,theta,eps,rho,Nsims,Nsteps,psicrit)
%HESTON_SIM Simulates the dynamics of a Heston Stochastic volatility model
%
%Usage [X,V]=Heston_Sim(x0,v0,r,T,kappa,theta,eps,rho,Nsims,Nsteps,psicrit)

if nargin < 11
    psicrit=1.5;
end

V=zeros(Nsims,Nsteps+1);
logX=zeros(Nsims,Nsteps+1);
logX(:,1)=ones(Nsims,1)*log(x0);
V(:,1)=ones(Nsims,1)*v0;
t=linspace(0,T,Nsteps+1); dt=t(2)-t(1);

inv_cumulated=@(u,p,b) log((1-p)./(1-u))./b.*(p<u);

gam1=0.5;gam2=0.5;
K1=gam1*dt*(kappa*rho/eps-0.5)-rho/eps;
K2=gam1*dt*(kappa*rho/eps-0.5)+rho/eps;
K3=gam1*dt*(1-rho^2); K4=gam2*dt*(1-rho^2);

bigA=K2+0.5*K4;
gaus=randn(Nsims,Nsteps);


for i=1:Nsteps
    M=zeros(Nsims,1);
    m=theta+(V(:,i)-theta)*exp(-kappa*dt);
    ssquared=V(:,i)*eps^2*exp(-kappa *dt)/kappa*(1-exp(-kappa*dt))+...
        theta*eps^2/(2*kappa)*(1-exp(-kappa*dt))^2;
    psi=ssquared./m.^2;
    
    low=find(psi<=psicrit);
    high=find(psi>psicrit);
    
    if  ~isempty(low)
    bsqu=2./psi(low)-1+sqrt(2./psi(low)).*sqrt(2./psi(low)-1);
    a=m(low)./(1+bsqu);
    V(low,i+1)=a.*(sqrt(bsqu)+randn(length(low),1)).^2;
    M(low)=exp((bigA*bsqu.*a)./(1-2*bigA*a))./sqrt(1-2*bigA*a);
    end
    
    if ~isempty(high)
    p=(psi(high)-1)./(psi(high)+1);
    beta=(1-p)./m(high);
    V(high,i+1)=inv_cumulated(rand(length(high),1),p,beta);
    M(high)=p+beta.*(1-p)./(beta-bigA);
    end
    
    
    %%% logX part
    
    K0=-log(M)-(K1+0.5*K3)*V(:,i);
    
    logX(:,i+1)=logX(:,i)+r*dt+K0+K1*V(:,i)+K2*V(:,i+1)...
        +sqrt(K3*V(:,i)+K4*V(:,i+1)).*gaus(:,i);
%    logX(:,i+1)=logX(:,i)+K0+K1*V(:,i)+K2*V(:,i+1)...
%        +sqrt(K3*V(:,i)+K4*V(:,i+1)).*gaus(:,i);
    
end

X=exp(logX);

end
