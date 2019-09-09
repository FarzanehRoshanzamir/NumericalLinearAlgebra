%%%%%%%% Farzaneh Roshanzamir %%%%%%%%
%%%%%%%%% LU Decomposition %%%%%%%%%%%
%%%%%%%%%%% Crout Alorithm %%%%%%%%%%%


function [L,U]=Crout(A)
    n=size(A,1);
    U=zeros(n,n);
    L=eye(n,n);
    
    U(1,1)=A(1,1);
    L(2:n,1)=(1/U(1,1))*(A(2:n,1));
    
    for j=2:n
        for i=1:j
          U(i,j)=A(i,j)-(L(i,1:i-1)*U(1:i-1,j));
          L(j+1:n,j)=(1/U(j,j))*(A(j+1:n,j)-(L(j+1:n,1:j-1)*U(1:j-1,j)));
        end
    end
    
end
