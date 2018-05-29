function [gamma,mu,Sigma,pi,v,v1,difdif,mu_record]=EM(X,mu,sigma,mu_true_samp,pi,max_iter)
% [gamma,mu,Sigma,pi,v,v1,difdif,mu_record]=EM(X,mu,sigma,mu_true_samp,pi,max_iter)
% mu_true_samp is ground truth
% mu is initial

isplot = 1;
n = size(X,1);
s = 2/(max(pi)+min(pi));
[k,d] = size(mu_true_samp);

%mu=eye(d);
for i=1:k
    Sigma{i}=sigma*eye(d);
end
gamma=zeros(n,k);

%allperm = perms(1:k);
%dist = zeros(size(allperm,1),1);
dist_truth = sum(diag(pdist2(mu,mu_true_samp)));
dist_truth_old = dist_truth;

%mu_0=mu_true+1*randn(size(mu_true));

mu_record=cell(max_iter+1, 1) ;
mu_record{1}=mu;
for iter=1:max_iter
    %s=1/sqrt(iter);
    %y=zeros(n,k);
    muold=mu;
    %expectation
    for i=1:k
        gamma(:,i)=log(mvnpdf(X,mu(i,:),Sigma{i})*pi(i));  %pi is CONSIDERED
    end
    
    gamma_not_normalize=gamma;
    for i=1:k
        gamma(:,i)=1./sum(exp(gamma_not_normalize-gamma_not_normalize(:,i)*ones(1,k)),2);
    end
    gamma(isnan(gamma)) = 0 ;
    %gamma=exp(log(gamma+1e-300)-log(sum(gamma,2)*ones(1,k)+k*1e-300));
    
    
    %mu=gamma'*X/sum(gamma(:,i));
    if isplot,
        plot(mu(2,1),mu(2,2),'x-','color',[179,35,39]/255,'markersize',20,'linewidth',3);
        hold on;
        plot(mu(1,1),mu(1,2),'^-','color',[40,137,41]/255,'markersize',20,'linewidth',3)
        plot(mu(3,1),mu(3,2),'o-','color',[2,13,130]/255,'markersize',20,'linewidth',3)
        pause(0.001)
    end
    g = gradient_gmm(X,mu,sigma,gamma);
    
    for i=1:k
        
        mu(i,:)=mu(i,:)+s*g(i,:);
        
    end
    
    
%     for i = 1:size(allperm,1),
%         %dist(i) = max(max(abs(mu(allperm(i,:),:)-muold)));
%         dist(i) = sum(diag(pdist2(mu(allperm(i,:),:),muold)));
%         dist_truthV(i) = sum(diag(pdist2(mu(allperm(i,:),:),mu_true_samp)));
%     end
    dist_last = norm(diag(pdist2(mu,muold)));
    
    dist_truth = norm(diag(pdist2(mu,mu_true_samp)));
    
    %[dist_truth,I]=min(dist_truthV);
    popu_grad = 1/k*(mu_true_samp-mu);
    difdif(iter) = norm(reshape(g,k*d,1)-reshape(popu_grad,k*d,1));
    if (dist_last < 1e-6 && norm(g) < 1e-6),
        break;
    end
    v(iter) = dist_truth/dist_truth_old;
    v1(iter) = log(dist_truth);
    dist_truth_old = dist_truth;
    
    mu_record{iter+1} = mu;
end

end

