function Price=Merton_Call(S0,K,r,T,sigma,lambda,mu_ju,delta_ju,M,N,theta)

ctheta=1-theta;
Smin=0.5*S0*exp((r-sigma^2/2)*T-6*sigma*sqrt(T));
Smax=2*S0*exp((r-sigma^2/2)*T+6*sigma*sqrt(T));
ymin=log(Smin/S0);
ymax=log(Smax/S0);

y=linspace(ymin,ymax,N+1)'; dy=y(2)-y(1);
t=linspace(0,T,M+1)'; dt=t(2)-t(1);
nodes=y(2:end-1);

k=@(z) lambda/(delta_ju*sqrt(2*pi))*exp(-(z-mu_ju).^2/(2*delta_ju^2));

[cap_alpha,cap_lam,Bl,Bu]=get_alph_lam(k,min(ymin,mu_ju),max(ymax,mu_ju),2*N);
cap_alpha=0
cap_lam=0
under=sigma^2/(2*dy^2)-(r-sigma^2/2-cap_alpha)/(2*dy);
mid_nodt=sigma^2/dy^2+cap_lam;
upper=sigma^2/(2*dy^2)+(r-sigma^2/2-cap_alpha)/(2*dy);

un=-ctheta*under*ones(N-1,1);
mid=1/dt+ctheta*mid_nodt*ones(N-1,1);
up=-ctheta*upper*ones(N-1,1);

A=spdiags([un mid up],-1:1,N-1,N-1);

un=theta*under*ones(N-1,1);
mid=1/dt-theta*mid_nodt*ones(N-1,1);
up=theta*upper*ones(N-1,1);

B=spdiags([un mid up],-1:1,N-1,N-1);
V=exp(-r*T)*max(S0*exp(nodes)-K,0);
payoff1=@(S) 0;
payoff2=@(S,t) exp(-r*(T-t))*max(S0*exp(S)-K,0);

V=exp(-r*T)*max(S0*exp(nodes)-K,0);
BC=sparse(N-1,1);

J=sparse(N-1,1);

for n=1:M
    BC(end)=(ctheta*upper+theta*upper)*payoff2(ymax,t(n));
    %J=rhs_integral(k,Bl,Bu,y,[0;V;0],2*N,payoff2,t(n));
    
    V=A\(B*V+J+BC);

end

V=[0;V;payoff2(ymax,t(end))];
S=S0*exp(y);
C=V;

plot(S,C)

Price=interp1(S,C,S0,'spline');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cap_alpha,cap_lam,Bl,Bu]=get_alph_lam(k,Bl,Bu,NN)

tol=1e-8;
step=0.5;
while k(Bl)>tol
    Bl=Bl-step;
end

while k(Bu)>tol
    Bu=Bu+step;
end

dB=(Bu-Bl)/(NN-1);

z=Bl+dB*(0:NN-1)';

w=ones(1,NN);w(1)=0.5;w(end)=0.5;
w=w*dB;
cap_lam=w*k(z);
cap_alpha=w*((exp(z)-1).*k(z));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=rhs_integral(k,Bl,Bu,y,V,NN,pay,t)

dB=(Bu-Bl)/(NN-1);

z=Bl+dB*(0:NN-1)';

w=ones(1,NN);w(1)=0.5;w(end)=0.5;
w=w*dB;
J=zeros(length(y)-2,1);
for i=1:length(y)-2
    temp=cut_f(y,V,y(i+1)+z,pay,t);
    J(i)=w*(temp.*k(z));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value=cut_f(y,V,val,pay,t)

value=zeros(size(val));
index_in=val>=y(1) & val<=y(end);
index_up=val>y(end);

value(index_in)=interp1(y,V,val(index_in));
value(index_up)=pay(val(index_up),t);
end


