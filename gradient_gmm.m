function g=gradient_gmm(X,mu,sigma, wts)

n=size(X,1);
k=size(mu,1);
d=size(X,2);
g=zeros(k,d);

for i=1:k
   g(i,:)=mean(spdiags([zeros(n,1) wts(:,i) zeros(n,1)], -1:1, n, n)*(X-ones(n,1)*mu(i,:))/sigma^2,1);
end

end