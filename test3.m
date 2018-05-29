function [log_opt_err,mu_infer]=test3(X,ytrue,mu,M_D,Pi,max_iter,mu0)
% [log_opt_err,mu_infer]=test3(X,ytrue,mu,M_D,Pi,max_iter,mu0)
[k,d]=size(mu);
n=size(X,1);
k=size(mu,1); % #mixtures
mu_true_samp=mu;

t=1;
opts.max_iter = 200;
opts.tol = 1e-7;
opts.steptype = 'g';
opts.isplot=0;
[gamma,mu_infer,Sigma,pi,v,v1,difdif,mu_record]=EM(X,mu0,1/sqrt(t),mu_true_samp,Pi,max_iter);
log_opt_err=v1;

end