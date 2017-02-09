function tau = AISC_tau(pOverPy,tauType)


switch lower(tauType)
    case 'none'
        tau = 1.0;
    case 'aisc'
        if isnan(pOverPy) || isinf(pOverPy)
            tau = NaN;
            return
        end

        assert(pOverPy <= 0,'AISC_tau only valid for compressive loads, pOverPy = %g',pOverPy);
        pOverPy = -pOverPy; % Compressive loads positive in the equations below
        
        if pOverPy <= 0.5
            tau = 1.0;
        elseif pOverPy <= 1.0
            tau = 4*(pOverPy)*(1-(pOverPy));
        else 
            tau = 0.0;
        end
    case 'maleck'
        if isnan(pOverPy) || isinf(pOverPy)
            tau = NaN;
            return
        end

        assert(pOverPy <= 0,'AISC_tau only valid for compressive loads, pOverPy = %g',pOverPy);
        pOverPy = -pOverPy; % Compressive loads positive in the equations below
        
        if pOverPy <= 0.39
            tau = 1.0;
        elseif pOverPy <= 1.0
            tau = -2.724*(pOverPy)*log(pOverPy);
        else
            tau = 0.0;
        end
    case 'composite'
        tau = 0.8;
    otherwise
        error('Unknown tauType')
end

end

