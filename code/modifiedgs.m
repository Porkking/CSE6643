%wei fang
%math6643 project modified gs
%instructor: Dr. Haesun Park
%April 21, 2020

function [x_mgs] = modifiedgs(A,b)

[m,n] = size(A);
L = eye(n);
Q1 = eye(n);
%generated Q1 and L with modified gram-schmidt
for k = n:-1:1
    L(k,k) = sqrt(sum(A(:,k).^2));
    for i = 1:m
        A(i,k) = A(i,k)./L(k,k);
    end
    for j = 1 :k-1
        L(k,j) = sum(A(:,k).*A(:,j));
        for i = 1:m
            A(i,j) = A(i,j) - A(i,k)*L(k,j);
        end
    end
end
Q1 = A;
d = Q1'*b;
lowtriang = dsp.LowerTriangularSolver;
x_mgs = lowtriang(L,d);

end