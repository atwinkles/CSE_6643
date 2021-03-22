function Problem5(coeff,q)
% Problem5(coeff,q)
%
% A MATLAB function to solve Problem 5 of Homework 1 by Alexander Winkles
%
% The output is a series of graphs of the solutions plus average relative 
%   errors
%
% coeff : lambda value for function
% q     : pivoting option (1 == no pivot, 2 == partial pivot)


% various step sizes
vals = [200,400,800];
%vals = [200,400,800,1000,2000,4000];

% runs the gaussian eliminaiton for the step sizes
for i = 1:size(vals,2)
    tic
    x = GaussianElim(coeff,vals(i),q);
    toc
    
end

end

% x = GaussianElim(lambda,n,q)
%
% A MATLAB function to compute the solution to a differential equation
% using Gaussian elimination
%
% Returns a vector x of solutions
%
% lambda : the value of lambda for the diff eq
% n      : the step size for the computation
% q      : the pivot option

function x = GaussianElim(lambda,n,q)

% declares boundary conditions
alpha = 0;
ualpha = 0;
beta = 1;
ubeta = -2;

% creates step sizes
h = (beta - alpha)/(n+1);

% creates arrays to store values
A = zeros(n);
B = zeros(n,1);
x = zeros(n,1);

% generates the centeral difference array and b vector for Ax = b
for i=1:n
    A(i,i) = 2 + lambda * h * h;
    if (i + 1 <= n)
        A(i+1,i) = -1;
        A(i,i+1) = -1;
    end
    
    if (i == 1)
        B(i) = h*h*f(alpha + i*h) + ualpha;
    elseif (i == n)
        B(i) = h*h*f(alpha + i*h) + ubeta;
    else
        B(i) = h*h*f(alpha + i*h);
    end
end

for i=1:n
    steps(i) = alpha + i*h;
end

% Gaussian Elimination without pivoting

if (q == 1)
    
    U = A;
    L = eye(n);
    for i = 1:n-1
        for j = i+1:n
            L(j,i) = U(j,i)/U(i,i);
            for k=i:n
                U(j,k) = U(j,k) - L(j,i)*U(i,k);
            end
        end
    end
    
    B = L\B;
    
    x = U\B;

% Gaussian Elimination with partial pivoting
else
    
    U = A;
    L = eye(n);
    for k=1:n-1
        val = 0;
        count = 0;
        for i=k:n
            if abs(U(i,k)) > val
                val = abs(U(i,k));
                count = i;
            end
        end
        swap = U(count,k:n);
        U(count,k:n) = U(k,k:n);
        U(k,k:n) = swap;

        swap = L(count,1:k-1);
        L(count,1:k-1) = L(k,1:k-1);
        L(k,1:k-1) = swap;
        
        for j = k+1:n
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
        end
    end
    
    B = L\B;
    
    x = U\B;
    
    
end

% returns the average relative error and creates figure

fprintf('The average relative error for %i is %d.\n',n,error(lambda,n,x',alpha,h))

figure
plot(steps,x,'o');

end

% differential equation
function u = f(x)

u = 3*x - 1/2;

end

% xhat = error(lambda,n,vec,alpha,h)
%
% A MATLAB function that computes the relative error of a vector with its
% corresponding true value
%
% Returns a double representing the relative error
%
% lambda : the lambda value of the diff eq
% n      : the size of the vector
% vec    : the vector being studied
% alpha  : used for step size
% h      : used for step size
function xhat = error(lambda,n,vec,alpha,h)

x = zeros(n,1);
for i=1:n
    x(i) = alpha + i*h;
end
sum = 0;

if lambda == 2
    for i=1:size(x,2)
        val = (exp(-sqrt(2)*x(i))*(exp(sqrt(2)*x(i))*(1-6*x(i))-exp(2*sqrt(2)*x(i))-13*exp(sqrt(2)*(2*x(i)+1)) + exp(sqrt(2)*(x(i)+2))*(6*x(i)-1)+exp(2*sqrt(2))+13*exp(sqrt(2))))/(4*(exp(2*sqrt(2))-1));
        sum = sum + (val - vec(i))/val;
    end
    
else
    for i=1:size(x,2)
        val = 1/4*x(i)*(-2*x(i)*x(i) + x(i) - 7);
        sum = sum + (val - vec(i))/val;
    end
end

xhat = sum / n;

end