function letter = listLetter(i)
letters = 'abcdefghijklmnopqrstuvwxyz';
i1 = floor((i-1)/length(letters)) + 1;
if i1 >= 6
    error('Number too high for listLetter');
end
i2 = 1+rem((i-1),length(letters));
letter = char(ones(1,i1)*letters(i2));
end