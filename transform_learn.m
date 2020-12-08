clc;
clear all;
close all;
%% my demo file for transform learning 
%Ravishankar and Bresler IEEE Transactions on Signal Processing 2015 
I=im2double(imread('peppers.png'));
I=I(:,:,1);
patchsize=[15,15];
Ip=im2col(I,patchsize,'sliding');
rp=randperm(size(Ip,2));
train_data=Ip(:,rp(1:10000));
% d_init=dctmtx(patchsize(1)*patchsize(2));

iter=100;
spp=4;
coef1=[1e-4,1e-3,1e-2,1e-1,1,10,100];
%  for i=1
%      for j=1

%  train_data=train_data-ones(size(train_data,1),1)*mean(train_data);
 l3=0.01*(mean(mean(train_data.^2)));
l2=0.001*(mean(mean(train_data.^2)));
 d_init=randn(size(train_data,1));
 d_init=d_init./sqrt(sum(d_init.^2));
for m=1:iter
    proj=d_init'*train_data;
    proj_s=sort(abs(proj),1,'descend');
    proj_t=proj_s(spp,:);
    proj(abs(proj)<proj_t)=0;
        err(m)=mean2(abs(d_init*proj-train_data).^2);

    dat_cov=train_data*train_data';
    [U,S,V]=svd(dat_cov+l3*eye(size(train_data,1)));
    LL=U*(S^(0.5))*V';
    LL2=pinv(LL);
    [Qi,Si,R]=svd(LL2*train_data*proj');
    sdiag=diag(Si);
    gam_val=0.5*(sdiag+sqrt((sdiag.^2)+2*l2));
    d_init=((R*diag(gam_val)*Qi')*LL2)';
      d_init=d_init./(sqrt(sum(d_init.^2)));
    k=1;
    for p=1:patchsize(1)
        for q=1:patchsize(2)
            c{p,q}=reshape(d_init(:,k)./sqrt(sum(d_init(:,k).^2)),patchsize);
            k=k+1;
        end
    end
    subplot(1,2,1),imagesc(cell2mat(c));colormap('gray');caxis([-0.25,0.25]);
    subplot(1,2,2),plot(err,'linewidth',3);
    title(['L2 val:',num2str(l2),'L3 val:',num2str(l3)]);
    pause(0.03);
end

%      end
%  end
