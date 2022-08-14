function ret=Obj_Heckman_full(param)

global data data2
global simdraws

x=data2(:,1);
y=data2(:,2);
d=data2(:,3);

alpha_wt=zeros(3,1);
alpha_wt(1)=param(1);
alpha_wt(2)=param(2);
alpha_wt(3)=param(3);

delta=zeros(3,1);
delta(1)=param(4);
delta(2)=param(5);
delta(3)=param(6);

k=param(7);

sigma2 = max(abs(param(8)),0.00001);
sigma_xi=max(abs(param(9)),0.00001);

vector_ones=ones(size(data2,1),1);
latent=alpha_wt(1)*vector_ones+alpha_wt(2)*x+alpha_wt(3)*d; %mean utility WT
latent2=delta(1)*vector_ones+delta(2)*x+delta(3)*y;  %mean utility competitors

pd = makedist('Normal','mu',0,'sigma',1); 
%pd2 = makedist('Normal','mu',0,'sigma',sigma2);  
pd3 = makedist('Normal','mu',0,'sigma',sigma_xi);

prob=log(pdf(pd,(data2(:,4)-latent2(:,1))/sigma_xi)+vector_ones*0.00000001)+log(cdf(pd,(data2(:,4)-latent(:,1)-k*(data2(:,4)-latent2(:,1)))/sigma2)+vector_ones*0.00000001);
loglike=sum(prob); 
loglike=loglike+size(data2,1)*log(1/sigma_xi); 

vector_ones=ones(size(data,1),1);
x=data(:,1);
y=data(:,2);
d=data(:,3);
latent=alpha_wt(1)*vector_ones+alpha_wt(2)*x+alpha_wt(3)*d; %mean utility WT
latent2=delta(1)*vector_ones+delta(2)*x+delta(3)*y;  %mean utility competitors

prob_epsilon_unc=zeros(size(data,1),1);
for i=1:size(data,1)
    if (data(i,5))==0 %(uncensored)
       if cdf(pd3,data(i,4)-latent2(i))==0
       loglike=loglike+log(0.00000001)+log(1/(sigma2));
       else
       
       simdraws_aux=simdraws*sigma_xi;
       draws_xi=simdraws_aux(simdraws_aux(:,1)<data(i,4)-latent2(i));
       n_xi=size(draws_xi,1);
       
       if n_xi>=5
       column_ones=ones(n_xi,1);
       prob_epsilon_unc(i,1)=log(sum(pdf(pd,((data(i,4)*column_ones-latent(i)*column_ones-k*draws_xi))/sigma2))/n_xi*cdf(pd,(data(i,4)-latent2(i))/sigma_xi));
       loglike=loglike+prob_epsilon_unc(i,1);
       loglike=loglike+log(1/(sigma2)); 
       else
       loglike=loglike+log(0.00000001)+log(1/(sigma2));
       end
      
       end
    end
end

% Alternative (a bit slower)
% sim=100;
% column_ones=ones(1,sim);
% draws_xi=zeros(size(data,1),sim);
% prob_epsilon_unc=zeros(size(data,1),1);
% for i=1:size(data,1)
%     if (data(i,5))==0 %(uncensored)
%        if cdf(pd3,data(i,4)-latent2(i))==0
%        loglike=loglike+log(0.00000001)+log(1/(sigma2));
%        else
%        t_up = truncate(pd3,-inf,data(i,4)-latent2(i)); %for uncensored
%        draws_xi(i,:) = transpose(random(t_up,sim,1)); %for uncensored
%        prob_epsilon_unc(i,1)=log(sum(pdf(pd,((data(i,4)*column_ones-latent(i)*column_ones-k*draws_xi(i,:)))/sigma2))/sim*cdf(pd,(data(i,4)-latent2(i))/sigma_xi));
%        loglike=loglike+prob_epsilon_unc(i,1);
%        loglike=loglike+log(1/(sigma2)); 
%        end
%     end
% end


ret=-loglike;
