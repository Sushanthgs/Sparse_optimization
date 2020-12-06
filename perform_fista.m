function [cost_f,x]=perform_fista(fm,dat,lambda,iter)
%my implementation of the Fast Iterative Shrinkage Thresholding algorihthm (FISTA)
%A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems 
% Amir Beck, Marc Teboulle, SIAM journal on imaging sciences (2009)

init_guess=zeros(size(fm,2),1);
lr=1/max(max(abs(eigs(fm'*fm))));
disp('Doing fista');
sl=@(x,lam)(sign(x).*(max(abs(x)-lam,0)));

x=init_guess;
y=x;
tk=1;
for i=1:iter
    x=sl(y-lr*(fm'*(fm*y-dat)),lambda*1*lr);
      cost_f(i)=sum((fm*x-dat).^2)+lambda*lr*sum(abs(x));
%      cost_f(i)=sum((fm*x-dat).^2);
     tk1=(1+sqrt(1+4*tk.^2))/2;

     x=x+((tk-1)/tk1)*(x-y);
     y=x;
     tk=tk1;
     
end
x=sl(x,lambda*lr);
end