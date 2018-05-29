function [ lik ] = expect_loglik(X,mu)
loglik=0;
n=size(X,1);
for i=1:n
    loglik=loglik+log(mvnpdf(X(i,:),mu(1,:),eye(2))/3+mvnpdf(X(i,:),mu(2,:),eye(2))/3+...
    mvnpdf(X(i,:),mu(3,:),eye(2))/3);
end
lik=loglik/n;
end

