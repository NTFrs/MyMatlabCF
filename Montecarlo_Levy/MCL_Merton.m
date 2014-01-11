function [S,t]=MCL_Merton(S0,gamma,sigma,lambda,jump_mu,jump_delta,T,Nsteps,Nsims)
%MCL_MERTON simulates jump trajectory for a Kou JD Process till time T
%
%Usage [S,t]=MCL_Merton(S0,gamma,sigma,lambda,jump_mu,jump_delta,T,Nsteps,Nsims)
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

    for i_time=1:Nsteps
        X(i_sims,i_time+1)=X(i_sims,i_time)+gamma*dt+sigma*sqrt(dt)*gauss(i_sims,i_time);
    
        for i_jumps=1:NT(i_sims)
            if shift_jump_time(i_jumps)==t(i_time)
                X(i_sims,i_time+1)=X(i_sims,i_time+1)+jump_mu+jump_delta*randn(1,1);
            end
        end
    end
    
end

S=S0*exp(X);

end