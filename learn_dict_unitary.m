function ld=learn_dict_unitary(training_data,dict,spp,iter)
% my implementation of unitary dictionary learning 
%lesage and gribonval bimbot and benaroya ICASSP 2005
for i=1:iter
    proj=dict'*training_data;
    proj_s=sort(abs(proj),1,'descend');
    proj=proj.*(abs(proj)>=proj_s(spp,:));
            err(i)=mean2((dict*proj-training_data).^2);

    mat=proj*training_data';
    [U,~,V]=svd(mat);
    dict=V*U';

end
ld.err=err;
ld.coef=proj;
ld.dict=dict;
end