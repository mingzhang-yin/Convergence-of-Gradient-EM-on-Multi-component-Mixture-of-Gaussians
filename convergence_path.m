% plot convergence path
close all;
t = 1;
R = 5; Rmax = 5;
gamma = Rmax/R;
Pi = [.4,.3,.3];
%Pi=ones(3,1)/3;
alpha = R/(2*sqrt(3));
prior = 0;%no prior

mu = [0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension
[k,d] = size(mu);
n = 1000; % #total points
t0 = [1,1,1];     %t is sd
[X,ytrue] = gen_mixture(n,mu,t0,Pi);
figure;plot(X(:,1),X(:,2),'.','color',[173,170,173]/255,'markersize',20)
hold on
plot(mu(:,1),mu(:,2),'ks','markersize',20,'linewidth',3)
plot(mu(:,1),mu(:,2),'kx','markersize',20,'linewidth',3)
r = 0.5;
dev = r*R;
dev_direction1 = (mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
mu0 = zeros(size(mu));
mu0(3,:) = mu(3,:)+dev*dev_direction1;
mu0(2,:) = mu(2,:)-dev*dev_direction1;
mu0(1,:) = mu(1,:);
%plot(mu0(:,1),mu0(:,2),'s','markersize',4)
%pause(0.01)
max_iter = 200;
[ga,mu_infer,Sigma,weights,v,v1,difdif,mu_record] = EM(X,mu0,1/sqrt(t),mu,Pi,max_iter);
dist_vec = diag(pdist2(mu_infer,mu));
xlim([-5,5]);ylim([-5,5])
set(gca,'visible','off')
set(gca,'Position',[0 0 1 1])