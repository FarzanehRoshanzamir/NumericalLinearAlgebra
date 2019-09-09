%%%%%%%%%% NLA Implementations %%%%%%%%%%
%%%%%%%%%% Farzaneh Roshanzamir %%%%%%%%%
%%%%%%%% Gauss-Jordan Ellimination %%%%%%
% more efficient algorithm

function [D,c] = GaussJordanElimination2(A,b)
   % profile on
    tic
    [m,n]=size(A);
    if(m~=n)
        disp('Error: Matrix is not square')
        return
    end
    
    A(:,n+1)=b;
    j=1;
   while(j<n)
        if(~any(A(j:n,j)))
            disp('Error: There is no unique answer')
            return
        end
        r=find(A(j:n,j),1,'first')+(j-1);
        if(r~=j)
            tmp=A(r,:);
            A(r,:)=A(j,:);
            A(j,:)=tmp;
        end
        for i=j+1:n
            A(i,:)=A(i,:)-(A(i,j)/A(j,j))*A(j,:);
        end
        j=j+1;
    end
    
    if(A(n,n)==0)
        disp('Error: There is no unique answer')
        return
    end
    
    while(j>1)
        for i=1:j-1
            l=A(i,j)/A(j,j);
            A(i,j)=A(i,j)-l*A(j,j);
            A(i,n+1)=A(i,n+1)-l*A(j,n+1);
        end
        j=j-1;
    end
    
    D = A;
    c=A(:,n+1);
    toc
   % profile off
    %profile('info')
end
