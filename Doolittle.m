%%%%%%%% Farzaneh Roshanzamir %%%%%%%%
%%%%%%%%% LU Decomposition %%%%%%%%%%%
%%%%%%%%% Doolittle Alorithm %%%%%%%%%

function [L,U]=Doolittle(A)

    n=size(A,1);
    U=zeros(n,n);
    L=eye(n,n);

    U(1,1:n)=A(1,1:n);
    L(2:n,1)=(1/U(1,1))*A(2:n,1);
    for k=2:n-1
        U(k,k:n)=A(k,k:n)-(L(k,1:k-1)*U(1:k-1,k:n));
        L(k+1:n,k)=(1/U(k,k))*(A(k+1:n,k)-(L(k+1:n,1:k-1)*U(1:k-1,k)));
    end

    U(n,n)=A(n,n)-(L(n,1:n-1)*U(1:n-1,n));
end

