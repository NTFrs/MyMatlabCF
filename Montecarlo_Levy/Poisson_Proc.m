function JumpTimes=Poisson_Proc(lambda,T,FlagPlot)
%POISSON_PROC simulates jump times for a Poisson Process till time T
%
%lambda jump intensity
%T time to simulate
%FlagPlot determines if process is plotted
NT=icdf('Poisson',rand,lambda*T);

JumpTimes=sort(rand(NT,1)*T);

if FlagPlot
    t=linspace(0,T,252*T)';
    S=zeros(length(t),1);
    
    for i=1:NT
        [~,index]=min(abs(JumpTimes(i)-t));
        ind=index(1);
        S(ind+1:end)=S(ind+1:end)+1;
    end
plot(t,S);
end


end