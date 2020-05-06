%wei fang
%math6643 project givens
%instructor: Dr. Haesun Park
%April 21, 2020

function [x_givens] = givens(A,b)

[m,n] = size(A);

%construct givens method
for q = 2:m
    for p = 1:min(q-1,n)
        if A(p,p)~=0
        t = A(q,p)/A(p,p);
        c = 1/sqrt(1+t^2);
        s = c*t;
        J = eye(m);
        J(p,p) = c;
        J(q,q) = c;
        J(p,q) = s;
        J(q,p) = -1*s;
        else
        J = eye(m);
        J(p,p) = 0;
        J(q,q) = 0;
        J(p,q) = 1;
        J(q,p) = 1;
        end
        A = J*A;
        b = J*b;
    end
end

R = A(1:n,:);
%output
uptriang = dsp.UpperTriangularSolver;
%x_givens = R\b(1:size(A,2));
x_givens = uptriang(R,b(1:size(A,2)));
end