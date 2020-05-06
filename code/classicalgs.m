%wei fang
%math6643 project classical gs
%instructor: Dr. Haesun Park
%April 21, 2020

function [x_cgs] = classicalgs(A,b)

[m,n] = size(A);
L(n,n) = norm(A(:,n));
Q1(:,n) = A(:,n)/L(n,n);

%generated Q1 and L with classical gram-schmidt
for k = n-1:-1:1
    sum_rq = 0;
    for i = n:-1:k+1
        L(i,k) = Q1(:,i)'*A(:,k);
        sum_rq = sum_rq + L(i,k)*Q1(:,i);
    end
    z(:,k) = A(:,k) - sum_rq;
    L(k,k) = norm(z(:,k));
    Q1(:,k) = z(:,k)/L(k,k);
end

%output
lowtriang = dsp.LowerTriangularSolver;
x_cgs = lowtriang(L,Q1'*b);

end