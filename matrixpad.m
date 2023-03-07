function [B] = matrixpad(A,x)
% pad the matrix A with boundary x 
[m,n]=size(A);
B=ones(m+2,n+2)*x;
B(2:m+1,2:n+1)=A;
end

