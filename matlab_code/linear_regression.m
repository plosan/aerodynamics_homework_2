function [a, b, r2] = linear_regression(X, Y)

S = max(size(X,1), size(X,2));
Sx = sum(X);
Sy = sum(Y);
Sxx = X' * X;
Sxy = X' * Y;
D = S * Sxx - Sx^2;
a = (Sxx * Sy - Sx * Sxy) / D;
b = (S * Sxy - Sx * Sy) / D;
X_ = X - (Sx / S);
Y_ = Y - (Sy / S);
r2 = (X_' * Y_)^2 / ((X_' * X_)*(Y_' * Y_));

end