function [P,M] = AISC_H1_interaction_diagram(Pt,Pc,M,quadrant)
    if nargin < 4
        quadrant = 'all';
    end
    switch lower(quadrant)
        case 'all'
            P = vertcat(0,[0.2 1.0 0.2]'*Pc,0,[0.2 1.0 0.2]'*Pt,0);
            M = [-1.0 -0.9 0.0 0.9 1.0 0.9 0.0 -0.9 -1.0]'*M;
        case {'cp','comppos','compressionpositive'}                  
            P = [1.0 0.2 0.0]*Pc;
            M = [0.0 0.9 1.0]*M;
        otherwise
            error('Unknown quadrant');
    end
end