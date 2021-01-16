clc;
clear all;
close all;
%%my implementation of recovery using low rank minimization
% Cai, Candes, Shen 2005 (Singular value thresholding), Robust PCA

I=im2double(imread('peppers.png'));
I=I(:,:,1);
iter=800;
prj=rand(size(I))>0.75;
Ip=prj.*I;
Y=0*I;
sl=@(x,lam)(sign(x).*(max(abs(x)-lam,0)));
lam=200;
lr=1;
disp_flag=1;
err=zeros(iter,1);
for i =1:iter
    [U,S,V]=svd(Y);
    St=sl(S,lam*lr);
    rec=U*St*V';
    Y=Y+lr*prj.*(Ip-rec);
    err(i)=sum(sum(prj.*(Ip-rec).^2))/sum(sum(prj));
    if(disp_flag==1)
        subplot(1,2,1),imagesc(rec);axis('square');colormap('gray');
        subplot(1,2,2),semilogy(err,'LineWidth',3);axis('square');grid on;
        pause(0.03);
    end
end