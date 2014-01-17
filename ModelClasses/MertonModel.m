classdef MertonModel < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        trend
        sigma
        jumpIntensity
        jumpMean
        jumpVar
        
    end
    
    methods
        function MM=MertonModel(gam,sig,lam,mu,del)
            MM.trend=gam;
            MM.sigma=sig;
            MM.jumpIntensity=lam;
            MM.jumpMean=mu;
            MM.jumpVar=del;
        end
        
        function z=charExp(Mod,u)
            z=-Mod.sigma^2*u.^2/2+1i*Mod.trend*u+Mod.jumpIntensity*(exp(-Mod.jumpVar*u.^2/2+1i*Mod.jumpMean*u)-1);
        end
        
        function z=charExpNeutral(Mod,u)
            z=Mod.charExp(u)-Mod.charExp(-1i)*1i*u;
        end
        
        function makeNeutral(Mod)
            Mod.trend=-Mod.sigma^2/2-Mod.jumpIntensity*(exp(Mod.jumpVar^2/2+Mod.jumpMean)-1);
        end
        
        function [X,t]=simPaths(Mod,T,Nsims,Nsteps)
            
            t=linspace(0,T,Nsteps+1)'; dt=t(2)-t(1);
            
            
            NT=icdf('Poisson',rand(Nsims,1),Mod.jumpIntensity*T);
            
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
                    X(i_sims,i_time+1)=X(i_sims,i_time)+Mod.trend*dt+Mod.sigma*sqrt(dt)*gauss(i_sims,i_time);
                    
                    for i_jumps=1:NT(i_sims)
                        if shift_jump_time(i_jumps)==t(i_time)
                            X(i_sims,i_time+1)=X(i_sims,i_time+1)+Mod.jumpMean+Mod.jumpVar*randn(1,1);
                        end
                    end
                end
                
            end
            
        end
        
    end
    
end

