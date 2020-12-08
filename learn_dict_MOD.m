function [ld]=learn_dict_MOD(dict_init,Ip,spp,iter,meanflag,disp_flag,patchsize,h,w)
%my implementation of the MOD (method of optimal directions) Engan et al
%or iterative least squares (ILS) dictionary learning algorithm ICASSP 1999
dict=dict_init;
%  iter=20;
err=zeros(iter,1);
 for i=1:iter
 if(i==1)
  A=OMP_Par(dict,Ip,spp);
 else
   if(meanflag==1)
     Imn=Ip-dict(:,1)*A(1,:);
  Ap=OMP_Par(dict(:,2:end),Imn,spp);
  A=[A(1,:);Ap];
   else
         A=OMP_Par(dict,Ip,spp);

   end
 end
    err(i)=mean(mean((dict*A-Ip).^2));
     if(disp_flag==1)
      subplot(1,2,1),display_dict(dict,patchsize,h,w);
      subplot(1,2,2),semilogy(err,'LineWidth',3);
      pause(0.03);
     end
%   Acov=A(:,2:end)*A(:,2:end)';
if (meanflag==1)
  D_update=Ip*pinv(full(A(2:end,:)));
  D_update=[ones(size(D_update,1),1),D_update-mean(D_update)];
  D_update=D_update./(ones(size(D_update,1),1)*sqrt(sum(D_update.^2)));
  
else
    D_update=Ip*pinv(full(A));
%   D_update=[ones(size(D_update,1),1),D_update-mean(D_update)];
  D_update=D_update./(ones(size(D_update,1),1)*sqrt(sum(D_update.^2)));
end
D_update(isnan(D_update))=0;
  dict=D_update;
  
 end
  code_update=OMP_Par(D_update,Ip,spp);
  ld.dict=D_update;
  ld.coef=code_update;
  ld.err=err;
%   err=err;
%   err_update=mean(mean((D_update*code_update-Ip).^2));
end