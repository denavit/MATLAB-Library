function [C,S] = csvread2(filename, delimiter)
if nargin < 2
    delimiter = ',';
end

% Open the text file
fid = fopen(filename);

% Determine number of columns
C = textscan(fgetl(fid),'%q','delimiter',delimiter);
numColumns = length(C{1});
frewind(fid);

% Read the file using textscan
C = textscan(fid,repmat('%q',1,numColumns),'delimiter',delimiter,'CollectOutput',true);
C = C{1};

% Close the text file
fclose(fid);

% Create structure of data if requested
if nargout > 1
    if size(C,1) < 2
        S = [];
    else
        num = str2double(C);
        S = struct;
        for i = 1:numColumns
            % Determine header
            header = genvarname(C{1,i});
            assert(~isfield(S,header),'repeated header: %s',C{1,i});
            % Determine if the data is numeric
            if isequal(isnan(num(2:end,i)),cellfun('isempty',C(2:end,i)))
                % Data is numeric or empty
                S.(header) = num(2:end,i);
            else
                S.(header) = C(2:end,i);
            end
        end
    end
end

end
