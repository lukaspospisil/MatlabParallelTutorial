function fxy = evaluate_something( x,y )
% example - evaluate quadratic function

X = [x;y];
A = [3 1; 1 2];
b = [19;13];

fxy = 0.5*dot(A*X,X) - dot(b,X);


end

