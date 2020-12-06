function [cost_s,u]=perform_salsa(fm,dat,lambda,iter,mu)
% my implementation of the SALSA method ( Split Augmented Lagrangian Shrinkage Algorithm)
%Afonso, Figuerido, Bioucas dias
%% IEEE transactions on image processing 2010
init_guess=zeros(size(fm,2),1);
% lr=1/max(max(abs(eigs(fm'*fm))));
disp('Doing salsa');
sl=@(x,lam)(sign(x).*(max(abs(x)-lam,0)));
% s_vec=ones(1,size(dat,1));
% c_fun=@(u,lam,mu,fm,dat)(s_vec*((fm*u-dat).^2)+(lam/mu)*s_vec*(abs(u)));

% u=init_guess;
i_mat=inv(fm'*fm+mu*eye(size(fm,2)));
x=init_guess;
xp=fm\dat;
cost_s=zeros(iter,1);

for i=1:iter
    u=sl(x+xp,lambda/mu);
    x=i_mat*(fm'*dat+mu*(u-xp));
    xp=xp-u+x;
     cost_s(i)=sum((fm*u-dat).^2)+lambda/mu*sum(abs(u));
%     cost_s(i)=c_fun(u,lambda,mu,fm,dat);

end
end