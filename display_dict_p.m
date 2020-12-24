function []=display_dict_p(dict,patchsize,h,w)
p=1;
 nx=@(x)((x-mean(x))/(max(x)-min(x)));

for i=1:h
    for j=1:w
                c{i,j}=reshape(nx(dict(:,p)),patchsize);
%              c{i,j}=20*log10(abs(reshape(dict(:,p),patchsize)));
%                c{i,j}=20*log10(abs(hilbert(reshape(dict(:,p),patchsize))));
%              c{i,j}=c{i,j}-max(max(c{i,j}));
        p=p+1;
    end
end
%  c{1,1}=ones(patchsize);
imagesc(cell2mat(c));
  caxis([-0.5,0.5]);colormap('gray');
%   caxis([-35,0]);colormap('gray');
end
