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
lam=20;
lr=1;
lev=3;
wav='sym4';
for i=1:iter
%     [U,S,V]=svd(Y);
    [C,S]=wavedec2(Y,lev,wav);
    St=sl(C,lam*lr);
    rec=waverec2(St,S,wav);
%     rec=U*St*V';
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
%     [U,S,V]=svd(Y);
   [C,S]=wavedec2(Y,lev,wav);
    St=sl(C,lam*lr);
    tk1=(1+sqrt(1+4*tk^2))/2;
    rec=waverec2(St,S,wav);
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
mu=0.01;
lam=0.01;
X=0*Ip;
Y=Ip;
for i=1:iter
       [C,S]=wavedec2(X+L,lev,wav);

%         [U,S,V]=svd();
        St=sl(C,lam/mu);
        Z=waverec2(St,S,wav);

%         Z=U*St*V';
        X=(Y+mu*(Z-L))./(mu+prj);
        L=L+X-Z;
 err_admm(i)=sum(sum((prj.*(Z-Ip).^2)))/sum(sum(prj));
    subplot(1,3,1),imagesc(I);axis('square');colormap('gray');
    subplot(1,3,2),imagesc(Z);axis('square');colormap('gray');
    subplot(1,3,3),semilogy(err_orig);hold on;semilogy(err_fista);
    semilogy(err_admm);hold off;axis('square');
    pause(0.03);
end
%%
semilogy(abs(diff(err_orig))./err_orig(1:end-1));hold on;
semilogy(abs(diff(err_fista))./err_fista(1:end-1));hold on;
semilogy(abs(diff(err_admm))./err_admm(1:end-1));
