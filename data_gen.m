function [X]=data_gen(alpha,gamma)
mu=[0,2*gamma;-sqrt(3),-1;sqrt(3),-1]*alpha; % #row->#mixtures, #col->#dimension

n=1500; % #total points
t=1;    %covariance matrix diag
[X,X_sym,ytrue]=gen_mixture(n,mu,t);
X=[X;X_sym];
end