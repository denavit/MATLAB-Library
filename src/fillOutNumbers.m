function numbers = fillOutNumbers(peaks,rate)
% fillOutNumbers  Evenly space numbers between peaks.
%
% numbers = fillOutNumbers(peaks, rate)
%
% Examples
% --------
% 
% >> fillOutNumbers([0, 1, -1], 0.25)
% ans =
%          0    % Peak 1
%     0.2500
%     0.5000
%     0.7500
%     1.0000    % Peak 2
%     0.7500
%     0.5000
%     0.2500
%          0
%    -0.2500
%    -0.5000
%    -0.7500
%    -1.0000    % Peak 3
%
% >> fillOutNumbers([0, 1, -1; 1, 2, -2], 0.25)
% ans =
%          0    1.0000   -1.0000
%     0.2500    1.2500   -1.2500
%     0.5000    1.5000   -1.5000
%     0.7500    1.7500   -1.7500
%     1.0000    2.0000   -2.0000
%

% if peaks is a row vector, make it a column vector
if ( size(peaks,1) == 1 )
    peaks = peaks';
end

numPeaks = size(peaks,1);
numbers(1,:) = peaks(1,:);

for i = 1:numPeaks-1
    diff = peaks(i+1,:) - peaks(i,:);
    numSteps = max([2 1+ceil(max(abs(diff./rate)))]);
    numbersToAdd = superlinspace(peaks(i,:),peaks(i+1,:),numSteps);
    numbers = vertcat(numbers,numbersToAdd(2:end,:));
end

end


function y = superlinspace(a,b,n)
y = zeros(n,size(a,2));
for i = 1:size(a,2)
    y(:,i) = linspace(a(i),b(i),n);
end
end