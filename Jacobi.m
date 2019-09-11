%%%%%%%% Farzaneh Roshanzamir %%%%%%%%
%%%%%%%%%%% Jacobi method %%%%%%%%%%%%

function X = Jacobi(D,L,U,b,x0,n,k)
    X=zeros(n,n);
    for i=1:k
        x0=(inv(D))*(b+(L+U)*x0);
        X(:,i)=x0;
    end
end
