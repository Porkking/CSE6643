%wei fang
%math6643 project householder
%instructor: Dr. Haesun Park
%April 21, 2020

function [x_house] = householder(A,b)
[m,n] = size(A);

%generated a [0;L] matrix with householder
for i = 0:n-1
    x = A(1:m-i,n-i);
    v = x;
    ind = eye(length(v));
    
    %pivoting
    if v(end) == 0
        [x_max,x_ind] = max(v);
        inx = 1:length(v);
        inx(end) = x_ind;
        inx(x_ind) = length(v);
        ind = ind(inx,:);
    end
    ind_new = blkdiag(ind,eye(i));
    A = ind_new*A;
    b = ind_new*b;
    x = A(1:m-i,n-i);
    v = x;

    %household vector
    v(end) = v(end) + sign(v(end))*norm(x,2);
    p = eye(m);
    p(1:m-i,1:m-i) = eye(length(v)) - 2.*v*v'./(v'*v);
    A = p*A;
    b = p*b;
end

L = A(m-n+1:m,:);
x_house = L\b(m-n+1:end);
end