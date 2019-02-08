function data = reinf_bar_lookup(bar_size,type)

if nargin < 2
    type = 'imperial';
end

all_data = [
    3  10  0.376 0.375 0.11 1.178
    4  13  0.668 0.500 0.20 1.571
    5  16  1.043 0.625 0.31 1.963
    6  19  1.502 0.750 0.44 2.356
    7  22  2.044 0.875 0.60 2.749
    8  25  2.670 1.000 0.79 3.142
    9  29  3.400 1.128 1.00 3.544
    10 32  4.303 1.270 1.27 3.990
    11 36  5.313 1.410 1.56 4.430
    14 43  7.65  1.693 2.25 5.32 
    18 57 13.60  2.257 4.00 7.09 
    20 64 16.69  2.500 4.91 7.85];

% Check input
if ischar(bar_size)
    assert(bar_size(1) == '#','bar_size should start with #')
    bar_size = str2double(bar_size(2:end));
end

if isnumeric(bar_size)
	assert(bar_size == round(bar_size),'bad input');
else
    error('bar_size should be an integer or character array')
end

% Lookup data
switch lower(type)
    case {'imperial','us'}
        ind = find(all_data(:,1)==bar_size);
    case {'soft_metric','metric'}
        ind = find(all_data(:,2)==bar_size);
    otherwise
        error('Unknown type: %s',type);
end

data = struct;
data.imperial_bar_size      = sprintf('#%i',all_data(ind,1));
data.soft_metric_bar_size   = sprintf('#%i',all_data(ind,2));
data.weight                 = all_data(ind,3);
data.weight_units           = 'lb/ft';
data.diameter               = all_data(ind,4);
data.diameter_units         = 'in';
data.area                   = all_data(ind,5);
data.area_units             = 'in^2';
data.perimeter              = all_data(ind,6);
data.perimeter_units        = 'in';

end
    