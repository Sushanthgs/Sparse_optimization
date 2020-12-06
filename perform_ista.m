function [cost,x]=perform_ista(fm,dat,lambda,iter)
% my implementation for iterative shrinkage thresholding algorithm.
%A New TwIST: Two-Step Iterative Shrinkage/. Thresholding Algorithms for Image Restoration. 
%José M. Bioucas-Dias, Member, IEEE, and Mário A. T. Figueiredo
 %note : this is not TwIST
lr=1/max(max(abs(eigs(fm'*fm))));
init_guess=zeros(size(fm,2),1);
x=init_guess;
sl=@(x,lam)(sign(x).*(max(abs(x)-lam,0)));
disp('Doing ista');
for i=1:iter
   
    x=sl(x-lr*(fm'*(fm*x-dat)),lambda*1*lr);
   cost(i)=sum((fm*x-dat).^2);
   
end
end