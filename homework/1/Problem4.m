function Problem4(M,N)
% Problem4(M,N)
%
% A MATLAB function for Problem 4 of Homework 1 by Alexander Winkles
%
% Returns figures of various results
%
% M : size of the matrix
% N : the number of matrices computer per size M

% declares the starting matrix size and the step size for increasing 
% the matrix
start = 100;
step = 100;

% generates empty vectors to store average values per size M
averageOne = zeros(M/step,1);
averageTwo = zeros(M/step,1);
averageInf = zeros(M/step,1);
normRatio = zeros(M/step,1);
averageCond = zeros(M/step,1);

% loops for each size M
for i = start:step:M
    
    disp(i);
    
    % does N computations to obtain an average value per matrix size
    for j=1:N
        
        % generates a random matrix with mean of 0 and standard deviation
        % of matrixsize^(1/2)
        A = normrnd(0,i^(1/2),i,i);
        
        % computes various norms for random matrix
        averageOne(i/step) = averageOne(i/step) + norm(A,1);
        averageTwo(i/step) = averageTwo(i/step) + norm(A,2);
        averageInf(i/step) = averageInf(i/step) + norm(A,Inf);
        normRatio(i/step) = normRatio(i/step) + norm(A,2)/norm(A,Inf);
        averageCond(i/step) = averageCond(i/step) + cond(A,1);
    end
    
end

% averages out results found by dividing by N
averageOne = averageOne / N;
averageTwo = averageTwo / N;
averageInf = averageInf / N;
normRatio = normRatio / N;
averageCond = averageCond / N;


% plots results nicely
X = start:step:M;

figure
plot(X,averageOne);
title('Average One Norm Values');
figure
plot(X,averageTwo);
title('Average Two Norm Values');
figure
plot(X,averageInf);
title('Average Inf Norm Values');
figure
plot(X,normRatio);
title('Average Ratio of 2 Norm and Inf Norm');
figure
plot(X,averageCond);
title('Average Condition Number');


