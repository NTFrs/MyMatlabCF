function [X,t]=MCL_Kou(gamma,sigma,lambda,lambda_u,lambda_d,p_up,T,Nsteps,Nsims)
%MCL_KOU simulates jump trajectory for a Kou JD Process till time T
%
%
%Usage [S,t]=MCL_Kou(S0,gamma,sigma,lambda,lambda_u,lambda_d,p_up,T,Nsteps,Nsims)
%gamma JD drift
%sigma diffusion coefficien
%lambda jump intensity
%jump_mu mean of the jumps
%jump_delta variance of the jumps
%T horizon time
%Nsteps number of time discretization steps
%Nsims number of simulations

t=linspace(0,T,Nsteps+1)'; dt=t(2)-t(1);


NT=icdf('Poisson',rand(Nsims,1),lambda*T);

X=zeros(Nsims,Nsteps+1);
X(:,1)=0;
gauss=randn(Nsims,Nsteps);
for i_sims=1:Nsims        
    JumpTimes=sort(rand(NT(i_sims),1)*T);
    indexs=zeros(NT(i_sims),1);
    for i_jumps=1:NT(i_sims)
        [dum,indtemp]=min(abs(JumpTimes(i_jumps)-t(2:end)));
        indexs(i_jumps)=indtemp(1);
    end
    
    shift_jump_time=t(indexs+1);

    u_prob=rand(NT(i_sims),1); u_intens=rand(NT(i_sims),1);
    for i_time=1:Nsteps
        X(i_sims,i_time+1)=X(i_sims,i_time)+gamma*dt+sigma*sqrt(dt)*gauss(i_sims,i_time);
    
        for i_jumps=1:NT(i_sims)
            if shift_jump_time(i_jumps)==t(i_time)
                if u_prob(i_jumps)<=p_up
                    Y=icdf('Exp',u_intens(i_jumps),1/lambda_u);
                else
                    Y=-icdf('Exp',u_intens(i_jumps),1/lambda_d);
                end
                X(i_sims,i_time+1)=X(i_sims,i_time+1)+Y;
            end
        end
    end
    
end


end