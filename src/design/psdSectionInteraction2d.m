function [P,M] = psdSectionInteraction2d(obj,axis,quadrant,points)
% Section Interaction Determined from the Plastic Stress
% Distribution Method

[points,remain] = strtok(points,'-');
if strcmpi(points,'psd')    
    [points,remain] = strtok(remain,'-');
end

% Continuous 
if strcmpi(points,'continuous')
    psd = obj.plasticStressDistributionObject;
    switch lower(axis)
        case {'x','z','major','strong'}
            [P,M,~] = psd.interactionSweep(0,100);
        case {'y','minor','weak'}
            [P,~,M] = psd.interactionSweep(pi/2,100);
        otherwise
            error('Bad axis');
    end
    return
end

% Compression Positive
P_cp = [];
M_cp = [];
if ismember('a',lower(points))
    assert(ismethod(obj,'pointA'),'obj does not have a method for point A');
    [Pi,Mi] = obj.pointA(axis);
    P_cp = vertcat(P_cp,Pi);
    M_cp = vertcat(M_cp,Mi);
end
if ismember('e',lower(points))
    assert(ismethod(obj,'pointE'),'obj does not have a method for point E');
    [Pi,Mi] = obj.pointE(axis);
    P_cp = vertcat(P_cp,Pi);
    M_cp = vertcat(M_cp,Mi);
end
if ismember('c',lower(points))
    assert(ismethod(obj,'pointC'),'obj does not have a method for point C');
    [Pi,Mi] = obj.pointC(axis);
    P_cp = vertcat(P_cp,Pi);
    M_cp = vertcat(M_cp,Mi);
end
if ismember('d',lower(points))
    assert(ismethod(obj,'pointD'),'obj does not have a method for point D');
    [Pi,Mi] = obj.pointD(axis);
    P_cp = vertcat(P_cp,Pi);
    M_cp = vertcat(M_cp,Mi);
end
if ismember('b',lower(points))
    assert(ismethod(obj,'pointB'),'obj does not have a method for point B');
    [Pi,Mi] = obj.pointB(axis);
    P_cp = vertcat(P_cp,Pi);
    M_cp = vertcat(M_cp,Mi);
end

switch lower(quadrant)
    case 'all'
        Pt = obj.Pnt;
        P = vertcat(P_cp,Pt ,flipud( P_cp));
        M = vertcat(M_cp,0.0,flipud(-M_cp));
    case {'cp','comppos','compressionpositive'}
        P = P_cp;
        M = M_cp;
    otherwise
        error('Unknown quadrant: %s',quadrant);
end

end