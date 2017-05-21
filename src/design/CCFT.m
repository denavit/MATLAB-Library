classdef CCFT < structural_shape
    %CCFT
    %
    % Notes:
    % - Formulas for plastic strength of the section are taken from
    % Geschwindner (2010) unless otherwise noted.
    %
    % REFERENCES:
    % - American Institute of Steel Construction (AISC) (2010). Specification
    % for Structural Steel Buildings, ANSI/AISC 360-10, American Institute
    % of Steel Construction, Chicago, Illinois.
    % - Roberto T Leon and Jerome F Hajjar, "Limit state response of composite
    % columns and beam-columns part II: Application of design provisions for
    % the 2005 AISC specification," Engineering Journal 45, no. 1 (2008): 21-46.
    % - Louis F. Geschwindner, "DISSCUSSION: Limit state response of composite
    % columns and beam-columns part II: Application of design provisions for
    % the 2005 AISC specification," Engineering Journal 47, no. 2 (2010): 131-140


    properties
        % Geometric Properties
        D        % Outer diameter of the steel tube
        t        % Thickness of the steel tube
        % Material Properties
        Fy       % Yield stress of the steel tube
        Fu = []; % Ultimate stress of the steel tube
        fc       % Compressive strength of the concrete core
        Es       % Modulus of Elasticity of Steel
        Ec       % Modulus of Elasticity of Concrete
        % Information
        shapeName = ''; % Name of the steel shape
        % Design Options
        C2 = 0.95;
        neglectLocalBuckling = false;
        option_EI = 'AISC2016';
    end

    methods
        %% Constructor
        function obj=CCFT(D,t,Fy,fc,units)
            obj.D  = D;
            obj.t  = t;
            obj.Fy = Fy;
            obj.fc = fc;
            obj.units = units;

            switch obj.units
                case 'US'
                    obj.Es = 29000;
                    obj.Ec = 57*sqrt(obj.fc*1000);
                case 'SI'
                    obj.Es = unitConvert('pressure',29000,'ksi','MPa');
                    fc_psi = unitConvert('pressure',obj.fc,'MPa','psi');
                    Ec_ksi = 57*sqrt(fc_psi);
                    obj.Ec = unitConvert('pressure',Ec_ksi,'ksi','MPa');
                otherwise
                    warning('Design:Input','Unknown unit system, setting Es = 0 and Ec = 0')
                    obj.Es = 0;
                    obj.Ec = 0;
            end

        end

        %% Geometric Properties
        function area = As(obj)
            area = (pi/4)*(obj.D)^2 - obj.Ac;
        end
        function area = Ac(obj)
            area = (pi/4)*(obj.D-2*obj.t)^2;
        end
        function area = Ag(obj)
            area = (pi/4)*(obj.D)^2;
        end
        function inertia = Is(obj,axis)
            inertia = (pi/64)*(obj.D)^4 - obj.Ic;
        end
        function inertia = Ic(obj,axis)
            inertia = (pi/64)*(obj.D-2*obj.t)^4;
        end
        function inertia = Ig(obj,axis)
            inertia = (pi/64)*(obj.D)^4;
        end
        function d = depth(obj,axis)
            d = obj.D;
        end
        function p = interfacePerimeter(obj)
            p = pi*(obj.D-2*obj.t);
        end

        %% Plastic Stress Distribution Anchor Points
        function [Pa,Ma] = pointA(obj,axis)
            Pa = -(obj.As*obj.Fy + obj.C2*obj.fc*obj.Ac);
            Ma = 0;
        end
        function [Pb,Mb] = pointB(obj,axis)
            h = obj.D-2*obj.t;
            rm = (obj.D-obj.t)/2;
            Kc = obj.fc*h^2;
            Ks = obj.Fy*rm*obj.t;
            theta = (0.0260*Kc-2*Ks)/(0.0848*Kc)+...
                sqrt((0.0260*Kc+2*Ks)^2 + 0.857*Kc*Ks)/(0.0848*Kc);
            % @todo - check if theta depends on C2
            ZsB = (obj.D^3-h^3) / 6 * sin (theta/2);
            ZcB = (h*sin(theta/2))^3 / 6;
            Mb = ZsB*obj.Fy + 0.5*ZcB*obj.C2*obj.fc;
            Pb = 0;
        end
        function [Pc,Mc] = pointC(obj,axis)
            [~,Mb] = obj.pointB;
            Pc = -(obj.C2*obj.fc*obj.Ac);
            Mc = Mb;
        end
        function [Pd,Md] = pointD(obj,axis)
            h = obj.D-2*obj.t;
            Zc = h^3/6;
            Zs = obj.D^3/6 - Zc;
            Pd = -(obj.C2*obj.fc*obj.Ac/2);
            Md = Zs*obj.Fy + 0.5*Zc*obj.C2*obj.fc;
        end
        function [Pe,Me] = pointE(obj,axis,type)
            if nargin < 3
                type = 'proposed';
            end
            [Pa,~] = obj.pointA;
            h = obj.D-2*obj.t;

            % From Point B
            rm = (obj.D-obj.t)/2;
            Kc = obj.fc*h^2;
            Ks = obj.Fy*rm*obj.t;
            theta = (0.0260*Kc-2*Ks)/(0.0848*Kc)+...
                sqrt((0.0260*Kc+2*Ks)^2 + 0.857*Kc*Ks)/(0.0848*Kc);
            % @todo - check if theta depends on C2
            hn = min((h/2)*sin((pi-theta)/2),h/2);

            he = hn/2 + h/4;
            theta2 = pi - 2*asin(2*he/h);
            ZsE = ((obj.D^3-h^3)/6)*sin(theta2/2);
            ZcE = (h^3/6)*sin(theta2/2)^3;

            switch type
                case 'original'
                    ZsE = ((obj.D^3-h^3)/6)*sin(theta2/2)^(4/3);
                    
                    %X = theta2/(theta2-sin(theta2)) + (2*pi-theta2)/((2*pi-theta2)-sin(2*pi-theta2));
                    %ZsE = ((obj.D^3-h^3)/6)*sin(theta2/2)^3*X;
                    
                    Pe = -(-Pa - 0.5*(obj.Fy*(obj.D^2-h^2)+0.5*obj.C2*obj.fc*h^2)*(theta2-sin(theta2)));
                    Me = obj.Fy*ZsE + 0.5*obj.C2*obj.fc*ZcE;
                case 'current'
                    Pe = -(-Pa - 0.25*(obj.Fy*(obj.D^2-h^2)+0.5*obj.C2*obj.fc*h^2)*(theta2-sin(theta2)));
                    Me = obj.Fy*ZsE + 0.5*obj.C2*obj.fc*ZcE;
                case 'proposed'
                    Pe = -(-Pa - 0.25*obj.Fy*(obj.D^2-h^2)*theta2 - 0.125*obj.C2*obj.fc*h^2*(theta2-sin(theta2)));
                    Me = obj.Fy*ZsE + 0.5*obj.C2*obj.fc*ZcE;
                case 'exact'
                    theta2s = pi - 2*asin(2*he/obj.D);
                    ZsE = (obj.D^3/6)*sin(theta2s/2)^3 - ZcE;

                    AcT = 0.125*h^2*(theta2-sin(theta2));
                    AsT = 0.125*obj.D^2*(theta2s-sin(theta2s)) - AcT;
                    Pe = -(-Pa - 2*obj.Fy*AsT - obj.C2*obj.fc*AcT);
                    Me = obj.Fy*ZsE + 0.5*obj.C2*obj.fc*ZcE;                    
                otherwise
                    error('Unknown type: %s',type);
            end
        end   
        
        %% Section Slenderness
        function lambda = lambda(obj,axis,component)
            if strcmpi(component,'tube')
                lambda = obj.D/obj.t;
            else
                error('Unknown axis: %s and/or component: %s',axis,component)
            end
        end
        function [slenderness,lambda,lambdaP,lambdaR,lambdaM] = ...
                slendernessInCompression(obj)
            lambda = obj.D/obj.t;
            lambdaP = 0.15*obj.Es/obj.Fy;
            lambdaR = 0.19*obj.Es/obj.Fy;
            lambdaM = 0.31*obj.Es/obj.Fy;
            if ( lambda <= lambdaP )
                slenderness = 'compact';
            elseif ( lambda <= lambdaR )
                slenderness = 'noncompact';
            elseif ( lambda <= lambdaM )
                slenderness = 'slender';
            else
                slenderness = 'notpermitted';
            end
        end
        function [slenderness,lambda,lambdaP,lambdaR,lambdaM] = ...
                slendernessInFlexure(obj)
            lambda = obj.D/obj.t;
            lambdaP = 0.09*obj.Es/obj.Fy;
            lambdaR = 0.31*obj.Es/obj.Fy;
            lambdaM = 0.31*obj.Es/obj.Fy;
            % Since lambdaR = lambdaM, no section should be classified as
            % slender
            if ( lambda <= lambdaP )
                slenderness = 'compact';
            elseif ( lambda <= lambdaR )
                slenderness = 'noncompact';
            else
                slenderness = 'notpermitted';
            end
        end
        function tf = isCompact(obj)
            tf = strcmp(obj.slendernessInCompression,'compact') && ...
                strcmp(obj.slendernessInFlexure,'compact');
        end

        %% Design Strengths
        function pnco = Pnco(obj)
            % Stub Column (L=0) Strength
            if obj.neglectLocalBuckling
                [Pa,~] = obj.pointA;
                pnco = -Pa;
            else
                [slenderness,lambda,lambdaP,lambdaR] = ...
                    obj.slendernessInCompression;
                [Pa,~] = obj.pointA;
                Pp = -Pa;
                switch slenderness
                    case 'compact'
                        pnco = Pp;
                    case 'noncompact'
                        Py = obj.Fy*obj.As + 0.7*obj.fc*obj.Ac;
                        pnco = Pp - (Pp-Py)*(lambda-lambdaP)^2/(lambdaR-lambdaP)^2;
                    case 'slender'
                        Fcr = 0.72*obj.Fy / ((obj.D/obj.t)*(obj.Fy/obj.Es))^0.2;
                        pnco = Fcr*obj.As + 0.7*obj.fc*obj.Ac;
                    case 'notpermitted'
                        pnco = 0;
                end
            end
        end
        function C3 = C3(obj)
            switch obj.option_EI
                case 'AISC2010'
                    C3 = min(0.9,0.6+2*(obj.As/(obj.Ac+obj.As)));
                case 'AISC2016'
                    C3 = min(0.9,0.45+3*obj.As/obj.Ag);
                otherwise
                    error('Unknown option_EI: %s',obj.option_EI);
            end            
        end
        function eieff = EIeff(obj,axis)
            eieff = obj.Es*obj.Is + obj.C3*obj.Ec*obj.Ic;
        end
        function x = stabilityReduction(obj,axis,Po)
            % Stability Reduction
            if strcmpi(axis,'min')
                x = min([obj.stabilityReduction('strong',Po)
                    obj.stabilityReduction('weak',Po)]);
                return
            end
            Pe = pi^2*obj.EIeff(axis)/(obj.K(axis)*obj.L(axis))^2;
            x = AISC_column_curve(Po/Pe);
        end
        function pn = Pnc(obj,axis)
            % Compressive Strength
            Pnco = obj.Pnco;
            pn = Pnco*obj.stabilityReduction(axis,Pnco);
        end
        function pnt = Pnt(obj)
            % Tension Strength
            pnt = obj.As*obj.Fy;
        end
        function mno = Mno(obj,axis)
            % L=0 Moment Strength
            mno = obj.Mn(axis);
        end
        function mn = Mn(obj,axis)
            % Flexural Strength

            if obj.neglectLocalBuckling
                [~,Mp] = obj.pointB(axis); % Platic Moment (Point B)
                mn = Mp;
            else
                [slenderness,lambda,lambdaP,lambdaR] = ...
                    obj.slendernessInFlexure;
                [~,Mp] = obj.pointB;
                switch slenderness
                    case 'compact'
                        mn = Mp;
                    case 'noncompact'
                        My = yieldMoment(obj);
                        mn = Mp - (Mp-My)*((lambda-lambdaP)/(lambdaR-lambdaP));
                    case 'notpermitted'
                        mn = 0;
                        return;
                end
            end
        end
        function vn = Vn(obj,Lv)
            % Shear Strength

            % Based on the steel section alone
            if (nargin < 2)
                % Lv not given, assume Equation (G6-2a) does not control
                Fcr = 0.78*obj.Es/(obj.D/obj.t)^1.5;
            else
                Fcr1 = 1.60*obj.Es/(sqrt(Lv/obj.D)*(obj.D/obj.t)^1.25);
                Fcr2 = 0.78*obj.Es/(obj.D/obj.t)^1.5;
                Fcr = max([Fcr1 Fcr2]);
            end
            if (Fcr > 0.6*obj.Fy)
                Fcr = 0.6*obj.Fy;
            end
            vn = Fcr*obj.As/2;
        end

        %% Design Checks
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)

            % Resistance Factors
            phi_Pc = 0.75;
            phi_M  = 0.90;
            phi_Pt = 0.90;
            phi_V  = 0.90;

            % Compressive Load / Moment Interaction
            if obj.isCompact || obj.neglectLocalBuckling
                % Use plastic stress distribution method
                [P_pointCs,M_pointCs] = obj.pointC('strong');
                [~,M_pointCw] = obj.pointC('weak');
                x = obj.stabilityReduction('min',obj.Pnco);
                Pa = phi_Pc*x*obj.Pnco;
                Pc = phi_Pc*x*(-P_pointCs);
                Mcs = phi_M*M_pointCs;
                Mcw = phi_M*M_pointCw;
                ratio_PMc = aisc_PSD_interaction_check(P,Ms,Mw,Pa,Pc,Mcs,Mcw);
            else
                % Use Chapter H equations
                Mcs = phi_M*obj.Mn('strong');
                Mcw = phi_M*obj.Mn('weak');
                Pc = phi_Pc*obj.Pnc('min');
                ratio_PMc = aisc_H11_interaction_check(P,Ms,Mw,Pc,Mcs,Mcw);
            end

            % Tensile Load / Moment Interaction
            Pc = phi_Pt*obj.Pnt;
            ratio_PMt = aisc_H12_interaction_check(P,Ms,Mw,Pc,Mcs,Mcw);

            % Check shear
            if isempty(Vs) && isempty(Vw)
                ratio_V = 0;
            else
                if isempty(Vs)
                    V = Vw;
                elseif isempty(Vw)
                    V = Vs;
                else
                    V = sqrt(Vs.^2+Vw.^2);
                end
                ratio_V = max(abs(V))/(phi_V*obj.Vn);
            end

            % No check on torsion
            if ~isempty(T)
                error('Torsion strength check not impelemented');
            end

            ratio = max([ratio_PMc ratio_PMt ratio_V]);
        end
        function pass_tf = proportioningCheck(obj,checkType,varargin)
            switch lower(checkType)
                case 'steelratio'
                    pass_tf = 0.01 <= obj.As/(obj.As+obj.Ac);
                case 'compact'
                    pass_tf = obj.isCompact;
                case 'moderatelyductilesection'
                    pass_tf = obj.D/obj.t < 0.15*obj.Es/obj.Fy;
                    % Also shear strength based on steel section alone                    
                case 'highlyductilesection'
                    pass_tf = obj.D/obj.t < 0.076*obj.Es/obj.Fy;
                    % Also shear strength based on steel section alone                    
                case 'slendernessratio_klr'
                    limit = varargin{1};
                    KLr_strong = obj.K('strong')*obj.L('strong')/sqrt(obj.Is/obj.As);
                    KLr_weak = obj.K('weak')*obj.L('weak')/sqrt(obj.Is/obj.As);
                    pass_tf = max([KLr_strong KLr_weak]) < limit;
                otherwise
                    error('Unknown checkType');
            end
        end

        %% Interaction Strength
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case 'aisc'
                    if obj.isCompact || obj.neglectLocalBuckling
                        [P,M] = obj.sectionInteraction2d(axis,'psd-acbt',quadrant);
                    else
                        [P,M] = AISC_H1_interaction_diagram(...
                            obj.Pnt,-obj.Pnco,obj.Mn(axis),quadrant);
                    end                
                case 'h1.1'
                    [P,M] = AISC_H1_interaction_diagram(...
                        obj.Pnt,-obj.Pnco,obj.Mno(axis),quadrant);
                case 'psd'
                    [P,M] = psdSectionInteraction2d(obj,axis,quadrant,type);
                case 'aci'
                    scACI = obj.strainCompatibilityAciObject;
                    switch lower(axis)
                        case 'strong'
                            [P,M,~] = scACI.interactionSweep(0,50);
                        case 'weak'
                            [P,~,M] = scACI.interactionSweep(pi/2,50);
                        otherwise
                            error('Bad axis');
                    end 
                case 'factoredaci'
                    scACI = obj.strainCompatibilityAciObject;
                    switch lower(axis)
                        case 'strong'
                            [P,M,~,et] = scACI.interactionSweep(0,50);
                        case 'weak'
                            [P,~,M,et] = scACI.interactionSweep(pi/2,50);
                        otherwise
                            error('Bad axis');
                    end 
                    phi = 0.75+(et-0.002)*50;
                    phi = max(0.75,phi);
                    phi = min(0.9,phi);
                    P = phi.*P;
                    M = phi.*M;
                case 'varma'
                    assert(any(strcmpi(quadrant,{'cp','comppos','compressionpositive'})),...
                        'Unknown quadrant');
                    [cp,cm] = obj.VarmaParameters;
                    P = -[1 cp 0]'*obj.Pnco;
                    M =  [0 cm 1]'*obj.Mno(axis);
                case  'factoredvarma'
                    assert(any(strcmpi(quadrant,{'cp','comppos','compressionpositive'})),...
                        'Unknown quadrant');
                    [cp,cm] = obj.VarmaParameters;
                    P = -[1 cp 0]'*0.75*obj.Pnco;
                    M =  [0 cm 1]'*0.9*obj.Mno(axis);
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case 'aisc'
                    if obj.isCompact || obj.neglectLocalBuckling
                        [P,M] = obj.beamColumnInteraction2d(axis,'psdsimple-acbt',quadrant);
                    else
                        [P,M] = AISC_H1_interaction_diagram(...
                            obj.Pnt,-obj.Pnc(axis),obj.Mn(axis),quadrant);
                    end
                case 'factoredaisc'
                    [P,M] = obj.beamColumnInteraction2d(axis,'AISC',quadrant);
                    P = 0.75*P;
                    M =  0.9*M;                        
                case 'h1.1'
                    [P,M] = AISC_H1_interaction_diagram(...
                        obj.Pnt,-obj.Pnc(axis),obj.Mn(axis),quadrant);
                case 'factoredh1.1'
                    [P,M] = AISC_H1_interaction_diagram(...
                        0.9*obj.Pnt,-0.75*obj.Pnc(axis),0.9*obj.Mn(axis),quadrant);
                case 'psdsimple'
                    [P,M] = psdSectionInteraction2d(obj,axis,quadrant,type(11:end));
                    ind = find(P<0);
                    P(ind) = P(ind)*obj.stabilityReduction(axis,obj.Pnco);
                case 'aci'
                    [P,M] = obj.sectionInteraction2d(axis,'ACI',quadrant);
                case 'factoredaci'
                    [P,M] = obj.sectionInteraction2d(axis,'FactoredACI',quadrant);                    
                case 'varma'
                    assert(any(strcmpi(quadrant,{'cp','comppos','compressionpositive'})),...
                        'Unknown quadrant');
                    [cp,cm] = obj.VarmaParameters;
                    P = -[1 cp 0]'*obj.Pnc(axis);
                    M =  [0 cm 1]'*obj.Mn(axis);
                case  'factoredvarma'
                    assert(any(strcmpi(quadrant,{'cp','comppos','compressionpositive'})),...
                        'Unknown quadrant');
                    [cp,cm] = obj.VarmaParameters;
                    P = -[1 cp 0]'*0.75*obj.Pnc(axis);
                    M =  [0 cm 1]'*0.9*obj.Mn(axis);                    
                case 'proposed'
                    [P,M] = proposedBeamColumnInteraction2d(obj,axis,quadrant,type(10:end),false);
                case 'factoredproposed'
                    [P,M] = proposedBeamColumnInteraction2d(obj,axis,quadrant,type(18:end),true);
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [cp,cm] = VarmaParameters(obj)
            csr = (obj.As*obj.Fy)/(obj.Ac*obj.fc);
            cp = 0.27*csr^-0.4;
            if csr >= 0.5
                cm = max(1.10*csr.^-0.08,1.00);
            else
                cm = min(0.95*csr.^-0.32,1.67);
            end
        end        

        %% Export and Information Functions
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(...
                obj,axis,type)
            switch lower(type)
                case 'columnstrength'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    E = obj.Es;
                    I = obj.EIeff/E;
                    A = EAeff/E;
                case 'gross'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is + obj.Ec*obj.Ic;
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                case 'steelonly'
                    E = obj.Es;
                    I = obj.Is(axis);
                    A = obj.As;                      
                case 'proposed'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is + obj.C3*obj.Ec*obj.Ic;
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                case 'aci'
                    E = obj.Ec;
                    I = 0.7*obj.Ig;
                    A = obj.Ag;
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [E,A,Iz,Iy,GJ] = sectionPropertiesForElasticAnalysis3d(...
                obj,type)
            error('Not yet implemented');
            % @todo
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');

            yExtreme = obj.D/2;
            yExtremeConcrete = obj.D/2 - obj.t;

            switch lower(type)
                case 'maxcompressive'
                    strain = min(axialStrain+yExtreme*curvature,...
                        axialStrain-yExtreme*curvature);
                case 'maxtensile'
                    strain = max(axialStrain+yExtreme*curvature,...
                        axialStrain-yExtreme*curvature);
                case 'maxabsolute'
                    strain = max(abs(axialStrain+yExtreme*curvature),...
                        abs(axialStrain-yExtreme*curvature));
                case 'maxconcretecompressive'
                    strain = min(axialStrain+yExtremeConcrete*curvature,...
                        axialStrain-yExtremeConcrete*curvature);
                case 'maxconcretetensile'
                    strain = max(axialStrain+yExtremeConcrete*curvature,...
                        axialStrain-yExtremeConcrete*curvature);
                otherwise
                    error('Unknown type');
            end
        end
        function x = getSectionData(obj,type,axis)
            switch lower(type)
                case 'steelarea'
                    x = obj.As;
                case 'concretearea'
                    x = obj.Ac;
                case 'grossarea'
                    x = obj.Ag;
                case 'steelratio'
                    x = obj.As/obj.Ag;
                case 'steelstrength'
                    x = obj.Fy;
                case 'concretestrength'
                    x = obj.fc;
                case 'grossconcreteflexuralrigidity'
                    x = obj.Ec*obj.Ic;
                case 'grosssteelflexuralrigidity'
                    x = obj.Es*obj.Is;
                case 'ecig'
                    x = obj.Ec*obj.Ig;
                case {'acfc','grossconcretestrength'}
                    x = obj.Ac*obj.fc;
                case {'asfy','grosssteelstrength'}
                    x = obj.As*obj.Fy;
                case 'grosssectioncompressionstrength'
                    x = obj.As*obj.Fy + obj.fc*obj.Ac;
                case 'aci-pmax'
                    x = 0.85*(0.85*obj.fc*obj.Ac + obj.As*obj.Fy);
                otherwise
                    x = NaN;
            end
        end
        function description = sectionDescription(obj,flag)
            if nargin < 2
                flag = 1;
            end
            switch flag
                case 1
                    description = sprintf('CCFT D = %g, t = %g, Fy = %g, fc = %g',...
                        obj.D,obj.t,obj.Fy,obj.fc);
                otherwise
                    error('Unknown flag');
            end
        end
        function fs = fiberSectionObject(obj,idSteel,idConc)
            fs = fiberSection;
            fs.addCCFT(idSteel,idConc,obj.D,obj.t);
        end
        function psd = plasticStressDistributionObject(obj)
            idSteel = 1;
            idConc  = 2;
            fs = obj.fiberSectionObject(idSteel,idConc);
            psd = plastic_stress_distribution(fs);
            psd.addMaterial(idSteel,obj.Fy,-obj.Fy);
            psd.addMaterial(idConc,0.0,-obj.C2*obj.fc);
        end
        function ess = elasticSectionStiffnessObject(obj)
            idSteel = 1;
            idConc  = 2;
            fs = obj.fiberSectionObject(idSteel,idConc);
            ess = elastic_section_stiffness(fs);
            ess.addMaterial(idSteel,obj.Es);
            ess.addMaterial(idConc ,obj.Ec,'CompressionOnly');
        end
        function scACI = strainCompatibilityAciObject(obj)
            idSteel = 1;
            idConc  = 2;
            fs = obj.fiberSectionObject(idSteel,idConc);
            scACI = ACI_strain_compatibility(fs);
            scACI.addConcreteBoundary(0,0,obj.D/2-obj.t);
            scACI.addSteelBoundary(0,0,obj.D/2);
            scACI.maxCompressiveStrength = -0.85*(obj.As*obj.Fy + 0.85*obj.Ac*obj.fc);
            scACI.addMaterial('steel',idSteel,obj.Fy,obj.Es);
            scACI.addMaterial('concrete',idConc,obj.fc,obj.units);
        end         
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            angles = linspace(0,2*pi,100);
            xo = (obj.D/2)*cos(angles);
            yo = (obj.D/2)*sin(angles);
            xi = (obj.D/2-obj.t)*cos(angles);
            yi = (obj.D/2-obj.t)*sin(angles);
            fill(xo,yo,obj.color_steelFill,'LineStyle','none')
            fill(xi,yi,obj.color_concreteFill,'LineStyle','none')
            plot(xo,yo,'k-','LineWidth',lineWidth);
            plot(xi,yi,'k-','LineWidth',lineWidth);
            axis equal
        end
    end

    methods (Static)
        function type = memberType()
            type = 'ccft';
        end
        function tf = hasConcrete()
            tf = true;
        end
        function tf = hasReinforcement()
            tf = false;
        end
        function t = t_given_rho(rho,D)
            t = (1-sqrt(1-rho))*D/2;
        end
    end

end




function my = yieldMoment(obj)

% Generate Fiber Data
fp_s = fiberPatch_circ(1,0,0,obj.D/2-obj.t,obj.D/2);
[~,As,ys] = fp_s.fiberData2d('strong',obj.D/400);
fp_c = fiberPatch_circ(1,0,0,0,obj.D/2-obj.t);
[~,Ac,yc] = fp_c.fiberData2d('strong',obj.D/400);

% Find Neutral Axis and Determine Moment
nf = length(Ac);
P = zeros(nf/2,1);
M = zeros(nf/2,1);
for i = 1:floor(nf/2)
    % Assume Neutal Axis
    yNA = yc(i);
    % Assign Steel Stresses
    Fs = zeros(size(As));
    fibs = (ys>yNA);
    Fs(fibs) = (ys(fibs)-yNA)/(obj.D/2-yNA) * -obj.Fy;
    fibs = (ys<yNA) & (ys>2*yNA-obj.D/2);
    Fs(fibs) = (yNA-ys(fibs))/(yNA-(2*yNA-obj.D/2)) * obj.Fy;
    fibs = (ys<2*yNA-obj.D/2);
    Fs(fibs) = ones(sum(fibs),1)*obj.Fy;
    % Assign Concrete Stresses
    Fc = zeros(size(Ac));
    fibs = (yc>yNA);
    Fc(fibs) = (yc(fibs)-yNA)/((obj.D/2-obj.t)-yNA) * -0.7*obj.fc;
    % Compute stress resultants
    P(i) = sum(Fs.*As)+sum(Fc.*Ac);
    M(i) = sum(Fs.*As.*-ys)+sum(Fc.*Ac.*-yc);
end

% Determine the best match
[forceErr,ind] = min(abs(P));
ay = obj.D/2-yc(ind);
my = M(ind) - forceErr*(obj.D/2-ay);

end

