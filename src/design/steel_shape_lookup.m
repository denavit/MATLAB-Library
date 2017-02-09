function data = steel_shape_lookup(label,type)

if nargin < 2
    % Determine Type
    first = strtok(label,'1234567890');
    if strcmpi(first,'W')
        type = 'Wide_Flange';
    elseif strcmpi(first,'HSS')
        numx = length(strfind(label,'x'));
        if numx == 1
            type = 'Round_HSS';
        elseif numx == 2
            type = 'Rectangular_HSS';
        else
            error('Could not determine section type for shape: %s',label)
        end
    else
        error('Could not determine section type for shape: %s',label)
    end
end

% Lookup section
switch type
    case 'Wide_Flange'
        load 'ShapeData_Wide_Flange.mat';
        labels = {ShapeData_Wide_Flange(:).label};
        ind = strcmpi(label,labels)==1;
        if sum(ind) ~= 1
            error('Could not find section: %s',label)
        end
        data = ShapeData_Wide_Flange(ind);
    case 'Rectangular_HSS'
        load 'ShapeData_Rectangular_HSS.mat';
        labels = {ShapeData_Rectangular_HSS(:).label};
        ind = strcmpi(label,labels)==1;
        if sum(ind) ~= 1
            error('Could not find section: %s',label)
        end
        data = ShapeData_Rectangular_HSS(ind);
    case 'Round_HSS'
        load 'ShapeData_Round_HSS.mat';
        labels = {ShapeData_Round_HSS(:).label};
        ind = strcmpi(label,labels)==1;
        if sum(ind) ~= 1
            error('Could not find section: %s',label)
        end
        data = ShapeData_Round_HSS(ind);
    otherwise
        error('Unknown type: %s',type)
end
        
end

