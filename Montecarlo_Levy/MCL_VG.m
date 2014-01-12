function [X,t]=MCL_VG(theta,sigma,kappa,T,Nsteps,Nsims)
%MCL_VG Generates paths of a VG levy process
%
%usage MCL_VG(S0,sigma,eta,kappa,T,Nsteps,Nsims)
%
%S0 asset value at time 0
%theta brownian drift
%sigma brownian volatility
%kappa subordinator volatility
%T time horizon
%Nsteps Number of time steps
%Nsims number of simulations

t=linspace(0,T,Nsteps+1); dt=t(2)-t(1);

a=dt/kappa;
dS=zeros(Nsims,Nsteps);

% if a<=1 % Johnk's
%     
%     for i=1:Nsims
%         
%         for j=1:Nsteps
%             A=rand(1,1)^(1/a);B=rand(1,1)^(1/(1-a));
%             while(A+B>1)
%                 A=rand(1,1)^(1/a);B=rand(1,1)^(1/(1-a));
%             end
%             dS(i,j)=A*icdf('Exp',rand(1,1),1)/(A+B);
%         end
%     end
%     
% else %Best's
%     b=a-1;c=3*a-3/4;
%     for i=1:Nsims
%         for j=1:Nsteps
%             neg=0;
%             
%             U=rand(1,1);
%             W=U*(1-U);
%             Y=sqrt(c/W)*(U-0.5);
%             A=b+Y;
%             if A<0
%                 neg=1;
%             else
%                 neg=0;
%             end
%             Z=64*W^3*rand(1,1)^3;
%             
%             while (log(Z)>2*(b*log(A/b)-Y) || neg==1)
%                 U=rand(1,1);
%                 W=U*(1-U);
%                 Y=sqrt(c/W)*(U-0.5);
%                 A=b+Y;
%                 if A<0
%                     neg=1;
%                 else
%                     neg=0;
%                 end
%                 Z=64*W^3*rand(1,1)^3;
%             end
%             
%             dS(i,j)=A;
%             
%         end
%     end
%     
%     
% end


%is this right? According to manual seems to be
dS=icdf('Gamma',rand(Nsims,Nsteps),a,1);

dS=kappa*dS;

X=zeros(Nsims,Nsteps+1);
gauss=randn(Nsims,Nsteps);

for j=1:Nsteps
    X(:,j+1)=X(:,j)+theta*dS(:,j)+sigma*gauss(:,j).*sqrt(dS(:,j));
end

end
