function numbers = fillOutNumbers(peaks,rate)

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