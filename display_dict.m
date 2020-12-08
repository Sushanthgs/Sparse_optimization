function []=display_dict(dict,patchsize,h,w)
p=1;
 nx=@(x)((x-mean(x))/(max(x)-min(x))); % normalize between 1 and -1
%display_dictionary 
for i=1:h
    for j=1:w
                 c{i,j}=reshape(nx(dict(:,p)),patchsize);
%
        p=p+1;
    end
end
%  c{1,1}=ones(patchsize);
imagesc(cell2mat(c));
 caxis([-0.5,0.5]);colormap('gray');
end
