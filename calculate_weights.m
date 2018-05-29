function gamma1=calculate_weights(X,mu_true,S,pi)
k=length(pi);
for i=1:k
    
    %         if sum(eig(Sigma{i})<0)>0
    %             keyboard
    %         end
    gamma1(:,i)=log(mvnpdf(X,mu_true(i,:),S)*pi(i));
end
%keyboard
%gamma=diag(1./(k*eps+sum(gamma,2)))*(gamma+eps);
gamma_not_normalize=gamma1;
for i=1:k
    gamma1(:,i)=1./sum(exp(gamma_not_normalize-gamma_not_normalize(:,i)*ones(1,k)),2);
end
gamma1(isnan(gamma1)) = 0 ;

end