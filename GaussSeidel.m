%%%%%%%% Farzaneh Roshanzamir %%%%%%%%
%%%%%%%% Gauss–Seidel method %%%%%%%%%

function X = Gauss_Seidel(D,L,U,b,x0,k,n)
    X=zeros(n,n);
    for i=1:k
      x0=(inv(D-L))*(b+U*x0);
      X(:,i)=x0;
    
    end
end
