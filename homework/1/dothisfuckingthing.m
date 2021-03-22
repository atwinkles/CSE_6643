function [A,B] = dothisfuckingthing(n,lambda, alpha,beta,ualpha,ubeta)

h = (beta - alpha)/(n+1);

A = zeros(n);
B = zeros(n,1);

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

end

function u = f(x)

u = 3*x - 1/2;

end