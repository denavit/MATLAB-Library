function [P,M] = trial_interaction_diagram(obj,axis,quadrant,type,factored)

if length(type) >= 7
    if strcmpi(type(1:6),'trial-')
        type = type(7:end);
    end
end
if length(type) >= 15
    if strcmpi(type(1:14),'factoredtrial-')
        type = type(15:end);
    end
end


switch lower(type)
    case 'acdb'

        switch lower(quadrant)
            case {'cp','comppos','compressionpositive'}
                [P,M] = obj.sectionInteraction2d(axis,'psd-acdb',quadrant);

                Pa = P(1);
                Pc = P(2);
                Pd = P(3);
                Mc = M(2);
                Md = M(3);
                Mb = M(4);

                Pnco = obj.Pnco;
                x = obj.stabilityReduction(axis,Pnco);

                if factored 
                    P = [0.75*x*Pa 0.75*x*Pc           0.75*x*Pd      0];
                    M = [        0   0.75*Mc 0.75*(Mb+x*(Md-Mb)) 0.9*Mb];
                else
                    P = [x*Pa x*Pc x*Pd 0];
                    M = [0 Mc (1-x)*Mb+x*Md Mb];
                end
            otherwise
                error('Unknown quadrant');
        end        
       
    otherwise
        error('Unknown type: %s',type);
end

end