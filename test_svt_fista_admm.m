clc;
clear all;
close all;
I=im2double(imread('peppers.png'));
I=I(:,:,1);
prj=(rand(size(I))>0.75);
Ip=I.*prj;
iter=200;
Y=0*Ip;
sl=@(x,lam)(sign(x).*(max(abs(x)-lam,0)));
lam=200;
lr=1;
for i=1:iter
    [U,S,V]=svd(Y);
    St=sl(S,lam*lr);
    rec=U*St*V';
    Y=Y-lr*prj.*(rec-Ip);
    err_orig(i)=sum(sum((prj.*(rec-Ip).^2)))/sum(sum(prj));
    subplot(1,3,1),imagesc(I);axis('square');colormap('gray');
    subplot(1,3,2),imagesc(rec);axis('square');colormap('gray');
    subplot(1,3,3),semilogy(err_orig);axis('square');
    pause(0.03);
end
%
tk=1;
Y_old=0*rec;
Y=Y_old;
rec=Y;
for i=1:iter
    [U,S,V]=svd(Y);
    St=sl(S,lam*lr);
    tk1=(1+sqrt(1+4*tk^2))/2;
    rec=U*St*V';
    Y=Y-lr*prj.*(rec-Ip);
    Y=Y+((tk-1)/tk1)*prj.*(Y-Y_old);
    tk=tk1;
    Y_old=Y;
    err_fista(i)=sum(sum((prj.*(rec-Ip).^2)))/sum(sum(prj));
    subplot(1,3,1),imagesc(I);axis('square');colormap('gray');
    subplot(1,3,2),imagesc(rec);axis('square');colormap('gray');
    subplot(1,3,3),semilogy(err_orig);hold on;semilogy(err_fista);hold off;axis('square');
    pause(0.03);
end
%

Y=0*rec;
Z=Y;
rec=Y;
L=rec;
mu=0.1;
lam=0.1;
X=0*Ip;
Y=Ip;
for i=1:iter
        [U,S,V]=svd(X+L);
        St=sl(S,lam/mu);
        Z=U*St*V';
        X=(Y+mu*(Z-L))./(mu+prj);
        L=L+X-Z;
 err_admm(i)=sum(sum((prj.*(Z-Ip).^2)))/sum(sum(prj));
    subplot(1,3,1),imagesc(I);axis('square');colormap('gray');
    subplot(1,3,2),imagesc(Z);axis('square');colormap('gray');
    subplot(1,3,3),semilogy(err_orig);hold on;semilogy(err_fista);
    semilogy(err_admm);hold off;axis('square');
    pause(0.03);
end