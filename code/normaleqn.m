%wei fang
%math6643 project normal equation
%instructor: Dr. Haesun Park
%April 21, 2020

function [x_ne] = normaleqn(A,b)

C = A'*A;
d = A'*b;
R = chol(C);
lowtriang = dsp.LowerTriangularSolver;
uptriang = dsp.UpperTriangularSolver;
y = lowtriang(R',d);
x_ne = uptriang(R,y);

end