function Price=American_BS_FD(S0,K,r,d,sigma,T,FlagCP,M,N,theta)

if (nargin<10)
    theta=0.5;
end

ctheta=1-theta;

Smin=S0*exp( (r-d-sigma^2/2)*T-sigma*sqrt(T)*6);
Smax=S0*exp( (r-d-sigma^2/2)*T+sigma*sqrt(T)*6);

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
omega=1.5;
MAXITER=25;
tol=1e-4;
for j=M:-1:1
    BC(1)=ctheta*BCund*r1(j)+theta*BCund*r1(j+1);
    BC(end)=ctheta*BCupp*r2(j)+theta*BCupp*r2(j+1);
    
    rhs=B*V+BC;
    Vold=V;
    
    for k=1:MAXITER
        for i=1:length(V)
            if i==1
                z=(rhs(i)-A(i,i+1)*Vold(i+1))/A(i,i);
            elseif i==length(V)
                z=(rhs(i)-A(i,i-1)*V(i-1))/A(i,i);
            else
                
                
                z=(rhs(i)-A(i,i-1)*V(i-1)-A(i,i+1)*Vold(i+1))/A(i,i);
            end
            V(i)=max(FlagCP*(nodes(i)-K),Vold(i)+omega*(z-Vold(i)));
        end
        if norm(V-Vold,Inf)<tol
            break;
        else
            Vold=V;
        end
    end
    if isempty(find(V==FlagCP*(nodes-K),1,'last'))==0
        F(j)=find(V==FlagCP*(nodes-K),1,'last');
    end
end

V=[r1(1);V;r2(1)];
figure
plot(S,V);
title('Price of option as function of S');
xlabel('Underlying S');
ylabel('Option price');

if exist('F','var')
figure
plot(t(1:end-1),nodes(F));
title('Frontiera');
end
Price=interp1(S,V,S0,'spline');


end
