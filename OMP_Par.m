function [A]=OMP_Par(D,X,L)
%my implementation of OMP for sparse coding (parallelized over patches)
% Joel Tropp and Anna Gilbert, IEEE transactions on Information theory
%=============================================
[~,P]=size(X);
[~,K]=size(D);
A=zeros(K,P);

parfor k=1:P
    a=[];
    x=X(:,k);
    residual=x;
%     indx=zeros(L,1); % uncomment to run sequentially
%     for j=1:1:L,
%         proj=D'*residual;
%         [maxVal,pos]=max(abs(proj));
%         pos=pos(1);
%         indx(j)=pos;
%         a=pinv(D(:,indx(1:j)))*x;
%         residual=x-D(:,indx(1:j))*a;
%         if sum(residual.^2) < 1e-6
%             break;
%         end
%     end;
[a,indx,j]=get_codes(D,residual,L,x); % comment to run sequentially
    temp=zeros(K,1);
    temp(indx(1:j))=a;
    A(:,k)=(temp);
%     disp(['Sparse Coding Progress: ',num2str(100*k/P),'%']);
end
A=sparse(A);
end
function [a,indx,j]=get_codes(D,residual,L,x) 
indx=zeros(L,1);
for j=1:1:L
        proj=D'*residual;
        [~,pos]=max(abs(proj));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
        if sum(residual.^2) < 1e-6
            break;
        end
end
end