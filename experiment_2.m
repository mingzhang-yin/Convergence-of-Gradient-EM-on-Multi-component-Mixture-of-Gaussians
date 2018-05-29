%% Experiment for contration region
clear;close all;clc;

%% Rmin=Rmax=2; pi=1/3; 
% figure 1 log(optimization error-iteration) plot
% Totally symmetric case, balanced \pi
% sigma = 1, change R, get different SNR

n_t = 1;
n_circle = 10; % number of circles

t = 1;
R = 2; Rmax = 2;
gamma = Rmax/R;
Pi = 1/3*ones(1,3);%[0.6,0.2,0.2]; %balanced

alpha = R/(2*sqrt(3)); 
delta = R/2/n_circle;
radius = exp(log(0.4)+linspace(log(1),log(0.5/0.4),n_circle));
prior = 0;%no prior

mu = [0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension
[k,d] = size(mu);
n = 30000; % #total points
t0 = [1,1,1];     %t is sd
[X,ytrue] = gen_mixture(n,mu,t0,Pi);
figure;plot(X(:,1),X(:,2),'.')
hold on
plot(mu(:,1),mu(:,2),'*','markersize',10)

max_iter = 200;    
%fix alpha do n_t times
%for tt=1:length(alpha_seq)
record = zeros(n_t,length(radius));
dev_direction=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
mu0=zeros(size(mu));
for tt=1:length(radius)
    r = radius(tt);
    dev = r*R;
    dev_direction1 = (mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
    mu0 = zeros(size(mu));
    mu0(3,:) = mu(3,:)+dev*dev_direction1;
    mu0(2,:) = mu(2,:)-dev*dev_direction1;
    if tt==length(radius),
        mu0(2,:) = (mu(2,:)+mu(3,:))/2;
        mu0(3,:) = (mu(2,:)+mu(3,:))/2;
    end
    mu0(1,:) = mu(1,:);
    %plot(mu0(:,1),mu0(:,2),'s','markersize',4)
    %pause(0.01)
    [ga,mu_infer,Sigma,weights,v,v1,difdif,mu_record] = EM(X,mu0,1/sqrt(t),mu,Pi,max_iter);
    record(tt)=exp(v1(end));
end

figure;plot(radius,record,'-o','LineWidth',1.5,'markersize',6)
ylabel('$\|\hat{\mu}-\mu^*\|$','Interpreter','latex','FontSize',16)
xlabel('a/Rmin','FontSize',16)

hold on;
% %error bar
% figure;
% mlog=mean(record);
% err=std(record);
% err( find( mod( 1:length(mlog), 5 ) > 0 ) ) = NaN;
% hb=errorbar(mlog,err,'bo-');
% legend(hb,'asymmetric','Location','best');
% xlabel('Iteration')
% ylabel('$\|\mu-\mu^*\|$','Interpreter','latex')

%% Exp 2
% Rmax = Rmin = 5;

R = 5; Rmax = 5;
gamma = Rmax/R;
Pi = 1/3*ones(1,3);%[0.6,0.2,0.2]; %balanced

alpha = R/(2*sqrt(3)); 
delta = R/2/n_circle;
mu=[0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension
[k,d] = size(mu);

t0 = [1,1,1];     %t is sd
[X,ytrue] = gen_mixture(n,mu,t0,Pi);
%figure;plot(X(:,1),X(:,2),'c.')
%hold on
%plot(mu(:,1),mu(:,2),'s','markersize',5)

max_iter = 200;    
%fix alpha do n_t times
%for tt=1:length(alpha_seq)
record=zeros(n_t,length(radius));
dev_direction=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
mu0=zeros(size(mu));
for tt=1:length(radius)
    r=radius(tt);
    dev=r*R;
    dev_direction1=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
    mu0=zeros(size(mu));
    mu0(3,:)=mu(3,:)+dev*dev_direction1;
    mu0(2,:)=mu(2,:)-dev*dev_direction1;
    mu0(1,:)=mu(1,:);
        if tt==length(radius),
        mu0(2,:) = (mu(2,:)+mu(3,:))/2;
        mu0(3,:) = (mu(2,:)+mu(3,:))/2;
    end
    [ga,mu_infer,Sigma,weights,v,v1,difdif,mu_record] = EM(X,mu0,1/sqrt(t),mu,Pi,max_iter);
    record(tt)=exp(v1(end));
end

plot(radius,record,'--s','LineWidth',1.5,'markersize',6);

%% Exp 3
% Rmax = 5; Rmin = 2;

R = 2; Rmax = 5;
gamma = Rmax/R;
Pi = 1/3*ones(1,3);%[0.6,0.2,0.2]; %balanced

alpha = R/(2*sqrt(3)); 
delta = R/2/n_circle;
prior = 0;%no prior

mu=[0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension
[k,d] = size(mu);

t0 = [1,1,1];     %t is sd
[X,ytrue] = gen_mixture(n,mu,t0,Pi);
%figure;plot(X(:,1),X(:,2),'c.')
%hold on
%plot(mu(:,1),mu(:,2),'s','markersize',5)

max_iter = 200;    
%fix alpha do n_t times
%for tt=1:length(alpha_seq)
record = zeros(n_t,length(radius));
dev_direction = (mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
mu0=zeros(size(mu));
for tt=1:length(radius)
    r=radius(tt);
    dev=r*R;
    dev_direction1=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
    mu0=zeros(size(mu));
    mu0(3,:)=mu(3,:)+dev*dev_direction1;
    mu0(2,:)=mu(2,:)-dev*dev_direction1;
    mu0(1,:)=mu(1,:);
    %plot(mu0(:,1),mu0(:,2),'s','markersize',4)
    %pause(0.01)
    [ga,mu_infer,Sigma,weights,v,v1,difdif,mu_record] = EM(X,mu0,1/sqrt(t),mu,Pi,max_iter);
    record(tt)=exp(v1(end));
end
record3 = record;
plot(radius,record3,'-.x','color',[40,137,41]/255,'LineWidth',1.5,'markersize',6);

%% Exp 3
% Rmax = 5; Rmin = 2;

R = 5; Rmax = 12.5;
gamma = Rmax/R;
Pi = 1/3*ones(1,3);%[0.6,0.2,0.2]; %balanced

alpha = R/(2*sqrt(3)); 
delta = R/2/n_circle;

mu = [0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension
[k,d] = size(mu);

t0 = [1,1,1];     %t is sd
[X,ytrue] = gen_mixture(n,mu,t0,Pi);
%figure;plot(X(:,1),X(:,2),'c.')
%hold on
%plot(mu(:,1),mu(:,2),'s','markersize',5)

max_iter = 200;    
%fix alpha do n_t times
%for tt=1:length(alpha_seq)
record=zeros(n_t,length(radius));
dev_direction=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
mu0=zeros(size(mu));
for tt=1:length(radius)
    r=radius(tt);
    dev=r*R;
    dev_direction1=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
    mu0=zeros(size(mu));
    mu0(3,:)=mu(3,:)+dev*dev_direction1;
    mu0(2,:)=mu(2,:)-dev*dev_direction1;
    mu0(1,:)=mu(1,:);
    %plot(mu0(:,1),mu0(:,2),'s','markersize',4)
    %pause(0.01)
    [ga,mu_infer,Sigma,weights,v,v1,difdif,mu_record] = EM(X,mu0,1/sqrt(t),mu,Pi,max_iter);
    record(tt)=exp(v1(end));
end

plot(radius,record,'-.*','LineWidth',1.5,'markersize',6);

%% Exp 5
% Rmax = Rmin = 2;
%{
R = 2; Rmax = 2;
gamma = Rmax/R;
Pi = [0.9,0.05,0.05]; 

alpha = R/(2*sqrt(3)); 
delta = R/2/n_circle;
mu=[0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension
[k,d] = size(mu);

t0 = [1,1,1];     %t is sd
[X,ytrue] = gen_mixture(n,mu,t0,Pi);
max_iter = 200;    
%fix alpha do n_t times
%for tt=1:length(alpha_seq)
record=zeros(n_t,length(radius));
dev_direction=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
mu0=zeros(size(mu));
for tt=1:length(radius)
    r=radius(tt);
    dev=r*R;
    dev_direction1=(mu(2,:)-mu(3,:))/norm(mu(2,:)-mu(3,:));
    mu0=zeros(size(mu));
    mu0(3,:)=mu(3,:)+dev*dev_direction1;
    mu0(2,:)=mu(2,:)-dev*dev_direction1;
    mu0(1,:)=mu(1,:);
    if tt==length(radius),
        mu0(2,:) = (mu(2,:)+mu(3,:))/2;
        mu0(3,:) = (mu(2,:)+mu(3,:))/2;
    end
    [ga,mu_infer,Sigma,weights,v,v1,difdif,mu_record] = EM(X,mu0,1/sqrt(t),mu,Pi,max_iter);
    record(tt)=exp(v1(end));
end

plot(radius,record,'-.v','LineWidth',1.5,'markersize',6);
%}
%% legend
leg = legend('Rmin=2; Rmax=2','Rmin=5; Rmax=5','Rmin=2; Rmax=5','Rmin=5; Rmax=12.5');

%leg = legend('Rmin=2; Rmax=2','Rmin=5; Rmax=5','Rmin=2; Rmax=5','Rmin=5; Rmax=12.5','Imbalance Rmin=2; Rmax=2');
set(leg,'FontSize',16,'Location','northwest');
set(gca,'FontSize',14);
ylim([-0.5,6])
xlim([.4,.5])