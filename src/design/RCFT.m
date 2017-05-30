classdef RCFT < structural_shape
    %RCFT
    %
    % Notes:
    % - Formulas for plastic strength of the section are taken from
    % Geschwindner (2010) unless otherwise noted.
    % - The internal radius of the fillet of steel tube is assumed to be
    % equal to the thickness of the steel tube
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
    % the 2005 AISC specification," Engineering Journal 47, no. 2 (2010): 131-140.

    properties
        % Geometric Properties
        H           % Height of the steel tube
        B           % Width of the steel tube
        t           % Thickness of the steel tube
        ri          % Internal radius of the corner of the steel tube
                    %   ri = 0 indicates a box section
        % Material Properties
        Fy          % Yield stress of the steel tube
        Fu = [];    % Ultimate stress of the steel tube
        fc          % Compressive strength of the concrete core
        Es          % Modulus of Elasticity of Steel
        Ec          % Modulus of Elasticity of Concrete
        vs = 0.3;   % Poisson's ratio of Steel
        vc = 0.2;   % Poisson's ratio of Concrete
        % Doubler Plate Properties
        t_doubler_strong  = [];
        Fy_doubler_strong = [];
        t_doubler_weak    = [];
        Fy_doubler_weak   = [];
        % Information
        shapeName = ''; % Name of the steel shape
        % Design Options
        C2 = 0.85;
        neglectLocalBuckling = false;
        useDefinedInsideCornerRadiusForSlenderness = false;
        option_EI = 'AISC2016';
    end

    methods
        %% Constructor
        function obj=RCFT(H,B,t,Fy,fc,units)
           	obj.H  = H;
            obj.B  = B;
            obj.t  = t;
            obj.Fy = Fy;
            obj.fc = fc;
            obj.units = units;
            
            % Default value of internal radius
            obj.ri = obj.t;

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
        function ro = ro(obj)
            if obj.ri == 0
                ro = 0;
            else
                ro = obj.ri+obj.t;
            end
        end
        function hc = Hc(obj)
            hc = obj.H - 2*obj.t;
        end
        function bc = Bc(obj)
            bc = obj.B - 2*obj.t;
        end
        function area = As(obj)
            area = obj.Ag - obj.Ac;
        end
        function area = Ac(obj)
            area = obj.Bc*obj.Hc - (4-pi)*obj.ri^2;
        end
        function area = Ag(obj)
            area = obj.B*obj.H - (4-pi)*obj.ro^2;
        end
        function i = Is(obj,axis)
            shp = Rectangular_Tube_Shape(obj.H,obj.B,obj.t,obj.ro);
            i = shp.I(axis);
        end
        function i = Ic(obj,axis)
            shp = Rectangle_Shape(obj.Hc,obj.Bc,obj.ri);
            i = shp.I(axis);
        end
        function i = Ig(obj,axis)
            shp = Rectangle_Shape(obj.H,obj.B,obj.ro);
            i = shp.I(axis);
        end
        function [d,b] = depth(obj,axis)
            switch lower(axis)
                case 'strong'
                    d = obj.H;
                    b = obj.B;
                case 'weak'
                    d = obj.B;
                    b = obj.H;
                otherwise
                    error('Bad axis');
            end
            if nargout < 2
                clear b;
            end
        end
        function [t_doubler,Fy_doubler] = doublerPlate(obj,axis)
            switch lower(axis)
                case 'strong'
                    t_doubler  = obj.t_doubler_strong;
                    Fy_doubler = obj.Fy_doubler_strong;
                case 'weak'
                    t_doubler  = obj.t_doubler_weak;
                    Fy_doubler = obj.Fy_doubler_weak;
                otherwise
                    error('Bad axis');
            end            
        end
        function p = interfacePerimeter(obj)
            p = 4*( (obj.H/2-obj.t-obj.ri) + (obj.B/2-obj.t-obj.ri) + (pi/2)*obj.ri );
        end
        function gs = Gs(obj)
            gs = obj.Es/2/(1+obj.vs);
        end
        function gc = Gc(obj)
            gc = obj.Ec/2/(1+obj.vc);
        end

        %% Plastic Stress Distribution Anchor Points
        function [Pa,Ma] = pointA(obj,axis)
            % Point A does not depend on the axis
            Pa = -(obj.As*obj.Fy + obj.C2*obj.fc*obj.Ac);
            Ma = 0;
        end
        function [Pb,Mb] = pointB(obj,axis)
            switch lower(axis)
                case 'strong'
                    [~,Md] = obj.pointD('strong');
                    h1 = obj.Bc;
                    h2 = obj.Hc;
                case 'weak'
                    [~,Md] = obj.pointD('weak');
                    h1 = obj.Hc;
                    h2 = obj.Bc;
                otherwise
                    error('Bad axis');
            end
            hn = min((obj.C2*obj.fc*obj.Ac)/(2*(obj.C2*obj.fc*h1+4*obj.t*obj.Fy)),h2/2);
            Zsn = 2*obj.t*hn^2;
            Zcn = h1*hn^2;
            Mb = Md - Zsn*obj.Fy - 0.5*Zcn*obj.C2*obj.fc;
            Pb = 0;
        end
        function [Pc,Mc] = pointC(obj,axis)
            switch lower(axis)
                case 'strong'
                    [~,Mb] = obj.pointB('strong');
                case 'weak'
                    [~,Mb] = obj.pointB('weak');
                otherwise
                    error('Bad axis');
            end
            Pc = -(obj.C2*obj.fc*obj.Ac);
            Mc = Mb;
        end
        function [Pd,Md] = pointD(obj,axis)
            concrete_shape = Rectangle_Shape(obj.Bc,obj.Hc,obj.ri);
            Zc = concrete_shape.Z(axis);
            steel_shape = Rectangular_Tube_Shape(obj.B,obj.H,obj.t,obj.ro);
            Zs = steel_shape.Z(axis);
            Pd = -(obj.C2*obj.fc*obj.Ac/2);
            Md = Zs*obj.Fy + 0.5*Zc*obj.C2*obj.fc;
        end
        function [Pe,Me] = pointE(obj,axis)
            switch lower(axis)
                case 'strong'
                    h1 = obj.Bc;
                    h2 = obj.Hc;
                    [~,Md] = obj.pointD('strong');
                case 'weak'
                    h1 = obj.Hc;
                    h2 = obj.Bc;
                    [~,Md] = obj.pointD('weak');
                otherwise
                    error('Bad axis');
            end
            % From Point B
            hn = min((obj.C2*obj.fc*obj.Ac)/(2*(obj.C2*obj.fc*h1+4*obj.t*obj.Fy)),h2/2);

            he = hn/2 + obj.H/4;
            ZsE = 2*obj.t*he^2;
            ZcE = h1*he^2;
            Pe = -(0.5*obj.C2*obj.fc*obj.Ac + obj.C2*obj.fc*h1*he + 4*obj.Fy*obj.t*he);
            Me = Md - obj.Fy*ZsE - 0.5*obj.C2*obj.fc*ZcE;
        end


        %% Section Slenderness
        function lambda = lambda(obj,axis,component)
            [d,b] = obj.depth(axis);
            if strcmpi(component,'flange')
                if obj.ri == 0
                    lambda = (b-2*obj.t)/obj.t;
                elseif obj.useDefinedInsideCornerRadiusForSlenderness
                    lambda = (b-2*(obj.t+obj.ri))/obj.t;
                else
                    lambda = (b-3*obj.t)/obj.t;
                end
            elseif strcmpi(component,'web')
                if obj.ri == 0
                    lambda = (d-2*obj.t)/obj.t;
                elseif obj.useDefinedInsideCornerRadiusForSlenderness
                    lambda = (d-2*(obj.t-obj.ri))/obj.t;
                else
                    lambda = (d-3*obj.t)/obj.t;
                end
            else
                error('Unknown axis: %s and/or component: %s',axis,component)
            end
        end
        function [slenderness,lambda,lambdaP,lambdaR,lambdaM] = ...
                slendernessInCompression(obj)
            lambda = max([...
                obj.lambda('strong','web') ...
                obj.lambda('strong','flange')]);
            lambdaP = 2.26*sqrt(obj.Es/obj.Fy);
            lambdaR = 3.00*sqrt(obj.Es/obj.Fy);
            lambdaM = 5.00*sqrt(obj.Es/obj.Fy);
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
                slendernessInFlexure(obj,axis)
            
            % Flange Slenderness Limits
            lambda_flg  = obj.lambda(axis,'flange');
            lambdaP_flg = 2.26*sqrt(obj.Es/obj.Fy);
            lambdaR_flg = 3.00*sqrt(obj.Es/obj.Fy);
            lambdaM_flg = 5.00*sqrt(obj.Es/obj.Fy);
            % Web Slenderness Limits
            lambda_web  = obj.lambda(axis,'web');
            lambdaP_web = 3.00*sqrt(obj.Es/obj.Fy);
            lambdaR_web = 5.70*sqrt(obj.Es/obj.Fy);
            lambdaM_web = 5.70*sqrt(obj.Es/obj.Fy); % i.e., no slender

            if ( lambda_web <= lambdaP_web && lambda_flg <= lambdaP_flg )
                slenderness = 'compact';
            elseif ( lambda_web <= lambdaR_web && lambda_flg <= lambdaR_flg )
                slenderness = 'noncompact';
            elseif ( lambda_web <= lambdaM_web && lambda_flg <= lambdaM_flg )
                slenderness = 'slender';
            else
                slenderness = 'notpermitted';
            end

            % @todo - use web values if they control
            lambda  = lambda_flg;
            lambdaP = lambdaP_flg;
            lambdaR = lambdaR_flg;
            lambdaM = lambdaM_flg; 
        end
        function tf = isCompact(obj)
            tf = strcmp(obj.slendernessInCompression,'compact') && ...
                strcmp(obj.slendernessInFlexure('strong'),'compact') && ...
                strcmp(obj.slendernessInFlexure('weak'),'compact');
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
                        Fcr = 9*obj.Es/lambda^2;
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
            eieff = obj.Es*obj.Is(axis) + obj.C3*obj.Ec*obj.Ic(axis);
        end
        function x = stabilityReduction(obj,axis,Po)
            % Stability Reduction
            if strcmpi(axis,'min')
                x = min([obj.stabilityReduction('strong',Po) ...
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
        function pn = Pnc_steel(obj,axis)
            % Compressive strength of the steel section only
            steelSection = RectangularHSS(obj.H,obj.B,obj.t,obj.Fy,obj.units);
            steelSection.ri = obj.ri;
            steelSection.Lx = obj.Lx;
            steelSection.Ly = obj.Ly;
            steelSection.Kx = obj.Kx;
            steelSection.Ky = obj.Ky;
            steelSection.Fu = obj.Fu;
            steelSection.Es = obj.Es;
            %steelSection.vs = obj.vs;
            steelSection.neglectLocalBuckling = obj.neglectLocalBuckling;
            pn = steelSection.Pnc(axis);
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
                switch lower(axis)
                    case 'strong'
                        d = obj.H;
                        b = obj.B;
                    case 'weak'
                        d = obj.B;
                        b = obj.H;
                    otherwise
                        error('Bad axis');
                end
                [slenderness,lambda,lambdaP,lambdaR] = ...
                    obj.slendernessInFlexure(axis);
                switch slenderness
                    case 'compact'
                        [~,Mp] = obj.pointB(axis); % Plastic Moment (Point B)
                        mn = Mp;
                    case 'noncompact'
                        [~,Mp] = obj.pointB(axis); % Plastic Moment (Point B)
                        bi = b-2*obj.t;
                        ay = (2*obj.Fy*d*obj.t + 0.35*obj.fc*bi*obj.t)/...
                            (4*obj.t*obj.Fy + 0.35*obj.fc*bi);
                        My = (0.35*obj.fc*(ay-obj.t)*bi)*((d-obj.t-(ay-obj.t)/3)-(d/2)) + ...
                            (bi*obj.t*obj.Fy)*((d-0.5*obj.t)-(d/2)) + ...
                            (ay*2*obj.t*0.5*obj.Fy)*((d-ay/3)-(d/2)) + ...
                            -(ay*2*obj.t*0.5*obj.Fy)*((d-ay-(2/3)*ay)-(d/2)) + ...
                            -((d-2*ay)*2*obj.t*obj.Fy)*(((d-2*ay)/2)-(d/2)) + ...
                            -(bi*obj.t*obj.Fy)*((0.5*obj.t)-(d/2));
                        mn = Mp - (Mp-My)*(lambda-lambdaP)/(lambdaR-lambdaP);
                        % @todo - which lambda web or flange
                    case 'slender'
                        Fcr = 9*obj.Es/lambda^2;
                        bi = b-2*obj.t;
                        acr = (obj.Fy*d*obj.t + (0.35*obj.fc+obj.Fy-Fcr)*bi*obj.t)/...
                            (obj.t*(Fcr+obj.Fy) + 0.35*obj.fc*bi);
                        mn = (0.35*obj.fc*(acr-obj.t)*bi)*((d-obj.t-(acr-obj.t)/3)-(d/2)) + ...
                            (bi*obj.t*Fcr)*((d-0.5*obj.t)-(d/2)) + ...
                            (acr*2*obj.t*0.5*Fcr)*((d-acr/3)-(d/2)) + ...
                            -((d-acr)*2*obj.t*0.5*obj.Fy)*(((d-acr)/3)-(d/2)) + ...
                            -(bi*obj.t*obj.Fy)*((0.5*obj.t)-(d/2));
                    case 'notpermitted'
                        mn = 0;
                end
            end
        end
        function vn = Vn(obj,axis)
            % Shear Strength

            % Based on the steel section alone
            lambda = obj.lambda(axis,'web');
            h = lambda*obj.t;
            Aw = 2*h*obj.t;
            kv = 5;
            if ( lambda <= 1.10*sqrt(kv*obj.Es/obj.Fy) )
                Cv = 1.0;
            elseif ( lambda <= 1.37*sqrt(kv*obj.Es/obj.Fy) )
                Cv = 1.10*sqrt(kv*obj.Es/obj.Fy)/lambda;
            else
                Cv = 1.51*obj.Es*kv/(lambda^2*obj.Fy);
            end
            vn = 0.6*obj.Fy*Aw*Cv;
        end
        function vpz = Vpz(obj,axis)
            switch obj.units
                % fcv[psi] = 28*sqrt(fc[psi])
                case 'US'
                    fcv = 0.885*sqrt(obj.fc); 
                case 'SI'
                    fcv = 2.325*sqrt(obj.fc);
                otherwise
                    error('Unknown unit system')
            end
            Vc  = fcv*obj.Ac;
            H = obj.depth(axis);
            Vs  = (0.6*obj.Fy)*H*2*obj.t;
            [t_doubler,Fy_doubler] = obj.doublerPlate(axis);
            if isempty(t_doubler) 
                Vdoubler = 0;
            else
                assert(~isempty(Fy_doubler),'Panel zone thickness defined but not yield strength');
                Vdoubler = (0.6*Fy_doubler)*(H-2*obj.ro)*t_doubler;
            end
            vpz = Vc+Vs+Vdoubler;
        end
        function kepz = Kepz(obj,axis,h,S)   
            [H,B] = obj.depth(axis);
            tanAlpha = h/H;
            Acp = 5/6*(S/tanAlpha)*(B-2*obj.t);
            kec = obj.Gc*Acp;
            kes = obj.Gs*H*2*obj.t;
            [t_doubler,~] = obj.doublerPlate(axis);
            if isempty(t_doubler)
                kedoubler = 0;
            else
                kedoubler = obj.Gs*(H-2*obj.ro)*t_doubler;
            end 
            kepz = kec+kes+kedoubler;
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

            % Check strong axis shear
            if isempty(Vs)
                ratio_Vs = 0;
            else
                ratio_Vs = max(abs(Vs))/(phi_V*obj.Vn('strong'));
            end

            % Check weak axis shear
            if isempty(Vw)
                ratio_Vw = 0;
            else
                ratio_Vw = max(abs(Vw))/(phi_V*obj.Vn('weak'));
            end

            % No check on torsion
            if ~isempty(T)
                error('Torsion strength check not impelemented');
            end

            ratio = max([ratio_PMc ratio_PMt ratio_Vs ratio_Vw]);
        end
        function pass_tf = proportioningCheck(obj,checkType,varargin)
            switch lower(checkType)
                case 'steelratio'
                    pass_tf = 0.01 <= obj.As/(obj.As+obj.Ac);
                case 'compact'
                    pass_tf = obj.isCompact;
                case 'highlyductilesection'
                    lambda = max([...
                        obj.lambda('strong','web') ...
                        obj.lambda('strong','flange')]);
                    pass_tf = lambda < 1.4*sqrt(obj.Es/obj.Fy);
                    % Also shear strength based on steel section alone
                case 'slendernessratio_klr'
                    limit = varargin{1};
                    KLr_strong = obj.K('strong')*obj.L('strong')/sqrt(obj.Is('strong')/obj.As);
                    KLr_weak = obj.K('weak')*obj.L('weak')/sqrt(obj.Is('weak')/obj.As);
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
                    phi = 0.65+(et-0.002)*250/3;
                    phi = max(0.65,phi);
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
            cp = 0.17*csr^-0.4;
            if csr >= 0.5
                cm = max(1.06*csr.^-0.11,1.00);
            else
                cm = min(0.90*csr.^-0.36,1.67);
            end
        end
        
        %% Export and Information Functions
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(...
                obj,axis,type)
            switch lower(type)
                case 'columnstrength'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    E = obj.Es;
                    I = obj.EIeff(axis)/E;
                    A = EAeff/E;
                case 'gross'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is(axis) + obj.Ec*obj.Ic(axis);
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                case 'steelonly'
                    E = obj.Es;
                    I = obj.Is(axis);
                    A = obj.As;                    
                case 'proposed'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is(axis) + obj.C3*obj.Ec*obj.Ic(axis);
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                case 'aci'
                    E = obj.Ec;
                    I = 0.7*obj.Ig(axis);
                    A = obj.Ag;
                case 'servicedefl'
                    EAeff = obj.Es*obj.As + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is(axis) + 0.6*obj.Ec*obj.Ic(axis);
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                otherwise
                    error('Unknown type');
            end
        end
        function [E,A,Iz,Iy,GJ] = sectionPropertiesForElasticAnalysis3d(...
                obj,type)
            error('Not yet implemented');
            % @todo
        end
        function lp = Lp(obj,axis,Li)
            lp = 0.17*Li;
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');
            switch lower(axis)
                case 'strong'
                    yExtreme = obj.H/2;
                    yExtremeConcrete = obj.Hc/2;
                case 'weak'
                    yExtreme = obj.B/2;
                    yExtremeConcrete = obj.Bc/2;
                otherwise
                    error('Bad axis');
            end
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
        function strain = longitudinalStrain3d(obj,axialStrain,curvatureY,curvatureZ,type)
            assert(isequal(size(axialStrain),size(curvatureY),size(curvatureZ)),...
                'axialStrain and curvature should be the same size');
            
            z = obj.B/2 - obj.ro;
            y = obj.H/2;
            strain_s1a = axialStrain + z*curvatureY + y*curvatureZ;
            strain_s2a = axialStrain + z*curvatureY - y*curvatureZ;
            strain_s3a = axialStrain - z*curvatureY - y*curvatureZ;
            strain_s4a = axialStrain - z*curvatureY + y*curvatureZ;
            
            z = obj.B/2 - (1-1/sqrt(2))*obj.ro;
            y = obj.H/2 - (1-1/sqrt(2))*obj.ro;
            strain_s1b = axialStrain + z*curvatureY + y*curvatureZ;
            strain_s2b = axialStrain + z*curvatureY - y*curvatureZ;
            strain_s3b = axialStrain - z*curvatureY - y*curvatureZ;
            strain_s4b = axialStrain - z*curvatureY + y*curvatureZ;
            
            z = obj.B/2;
            y = obj.H/2 - obj.ro;
            strain_s1c = axialStrain + z*curvatureY + y*curvatureZ;
            strain_s2c = axialStrain + z*curvatureY - y*curvatureZ;
            strain_s3c = axialStrain - z*curvatureY - y*curvatureZ;
            strain_s4c = axialStrain - z*curvatureY + y*curvatureZ;            
            
            z = obj.Bc/2 - obj.ri;
            y = obj.Hc/2;
            strain_c1a = axialStrain + z*curvatureY + y*curvatureZ;
            strain_c2a = axialStrain + z*curvatureY - y*curvatureZ;
            strain_c3a = axialStrain - z*curvatureY - y*curvatureZ;
            strain_c4a = axialStrain - z*curvatureY + y*curvatureZ;
            
            z = obj.Bc/2 - (1-1/sqrt(2))*obj.ri;
            y = obj.Hc/2 - (1-1/sqrt(2))*obj.ri;
            strain_c1b = axialStrain + z*curvatureY + y*curvatureZ;
            strain_c2b = axialStrain + z*curvatureY - y*curvatureZ;
            strain_c3b = axialStrain - z*curvatureY - y*curvatureZ;
            strain_c4b = axialStrain - z*curvatureY + y*curvatureZ;
            
            z = obj.Bc/2;
            y = obj.Hc/2 - obj.ri;
            strain_c1c = axialStrain + z*curvatureY + y*curvatureZ;
            strain_c2c = axialStrain + z*curvatureY - y*curvatureZ;
            strain_c3c = axialStrain - z*curvatureY - y*curvatureZ;
            strain_c4c = axialStrain - z*curvatureY + y*curvatureZ;           
            
            switch lower(type)
                case 'maxcompressive'
                    strain = min([...
                        strain_s1a(:) strain_s2a(:) strain_s3a(:) strain_s4a(:) ...
                        strain_s1b(:) strain_s2b(:) strain_s3b(:) strain_s4b(:) ...
                        strain_s1c(:) strain_s2c(:) strain_s3c(:) strain_s4c(:) ...
                        ],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxtensile'
                    strain = max([...
                        strain_s1a(:) strain_s2a(:) strain_s3a(:) strain_s4a(:) ...
                        strain_s1b(:) strain_s2b(:) strain_s3b(:) strain_s4b(:) ...
                        strain_s1c(:) strain_s2c(:) strain_s3c(:) strain_s4c(:) ...
                        ],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxabsolute'
                    strain = max(abs([...
                        strain_s1a(:) strain_s2a(:) strain_s3a(:) strain_s4a(:) ...
                        strain_s1b(:) strain_s2b(:) strain_s3b(:) strain_s4b(:) ...
                        strain_s1c(:) strain_s2c(:) strain_s3c(:) strain_s4c(:) ...
                        ]),[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxconcretecompressive'
                    strain = min([...
                        strain_c1a(:) strain_c2a(:) strain_c3a(:) strain_c4a(:) ...
                        strain_c1b(:) strain_c2b(:) strain_c3b(:) strain_c4b(:) ...
                        strain_c1c(:) strain_c2c(:) strain_c3c(:) strain_c4c(:) ...
                        ],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxconcretetensile'
                    strain = max([...
                        strain_c1a(:) strain_c2a(:) strain_c3a(:) strain_c4a(:) ...
                        strain_c1b(:) strain_c2b(:) strain_c3b(:) strain_c4b(:) ...
                        strain_c1c(:) strain_c2c(:) strain_c3c(:) strain_c4c(:) ...
                        ],[],2);
                    strain = reshape(strain,size(axialStrain));
                otherwise
                    error('Unknown type: %s',type);
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
                case 'concretestrength'
                    x = obj.fc;
                case 'steelstrength'
                    x = obj.Fy;
                case 'grossconcreteflexuralrigidity'
                    x = obj.Ec*obj.Ic(axis);
                case 'grosssteelflexuralrigidity'
                    x = obj.Es*obj.Is(axis);
                case 'ecig'
                    x = obj.Ec*obj.Ig(axis);
                case 'grosssectioncompressionstrength'
                    x = obj.As*obj.Fy + obj.fc*obj.Ac;
                case {'asfy','grosssteelstrength'}
                    x = obj.As*obj.Fy;
                case {'acfc','grossconcretestrength'}
                    x = obj.fc*obj.Ac;
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
                    if ~isempty(obj.shapeName)
                        description = regexprep(lower(obj.shapeName),'hss', 'RCFT');
                    else
                        description = sprintf('RCFT%gx%gx%.3f',...
                            obj.H,obj.B,obj.t);
                    end
                otherwise
                    error('Unknown flag');
            end
        end
        function fs = fiberSectionObject(obj,idSteel,idConc)
            fs = fiberSection;
            fs.addRCFT(idSteel,idConc,obj.H,obj.B,obj.t,obj.ri);
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
            if obj.ri == 0
                scACI.addConcreteBoundary(-obj.B/2+obj.t,-obj.H/2+obj.t,0);
                scACI.addConcreteBoundary(-obj.B/2+obj.t, obj.H/2-obj.t,0);
                scACI.addConcreteBoundary( obj.B/2-obj.t, obj.H/2-obj.t,0);
                scACI.addConcreteBoundary( obj.B/2-obj.t,-obj.H/2+obj.t,0);
                scACI.addSteelBoundary(-obj.B/2,-obj.H/2,0);
                scACI.addSteelBoundary(-obj.B/2, obj.H/2,0);
                scACI.addSteelBoundary( obj.B/2, obj.H/2,0);
                scACI.addSteelBoundary( obj.B/2,-obj.H/2,0);
            else
                scACI.addConcreteBoundary(-obj.B/2+obj.ro,-obj.H/2+obj.ro,obj.ri);
                scACI.addConcreteBoundary(-obj.B/2+obj.ro, obj.H/2-obj.ro,obj.ri);
                scACI.addConcreteBoundary( obj.B/2-obj.ro, obj.H/2-obj.ro,obj.ri);
                scACI.addConcreteBoundary( obj.B/2-obj.ro,-obj.H/2+obj.ro,obj.ri);                       
                scACI.addSteelBoundary(-obj.B/2+obj.ro,-obj.H/2+obj.ro,obj.ro);
                scACI.addSteelBoundary(-obj.B/2+obj.ro, obj.H/2-obj.ro,obj.ro);
                scACI.addSteelBoundary( obj.B/2-obj.ro, obj.H/2-obj.ro,obj.ro);
                scACI.addSteelBoundary( obj.B/2-obj.ro,-obj.H/2+obj.ro,obj.ro);                       
            end
            scACI.maxCompressiveStrength = -0.85*(obj.As*obj.Fy + 0.85*obj.Ac*obj.fc);
            scACI.addMaterial('steel',idSteel,obj.Fy,obj.Es);
            scACI.addMaterial('concrete',idConc,obj.fc,obj.units);
        end         
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            if obj.ri == 0
                xo = obj.B/2  * [1 -1 -1 1 1];
                yo = obj.H/2  * [1 1 -1 -1 1];
                xi = obj.Bc/2 * [1 -1 -1 1 1];
                yi = obj.Hc/2 * [1 1 -1 -1 1];
            else
                angles = linspace(0,pi/2,25);
                d1 = obj.H/2-obj.ro;
                b1 = obj.B/2-obj.ro;
                xo = [ ( b1+obj.ro*cos(angles)) ...
                    (-b1+obj.ro*cos(angles+pi/2)) ...
                    (-b1+obj.ro*cos(angles+pi)) ...
                    ( b1+obj.ro*cos(angles+1.5*pi)) b1+obj.ro];
                yo = [ ( d1+obj.ro*sin(angles)) ...
                    ( d1+obj.ro*sin(angles+pi/2)) ...
                    (-d1+obj.ro*sin(angles+pi)) ...
                    (-d1+obj.ro*sin(angles+1.5*pi)) d1 ];
                xi = [ ( b1+obj.ri*cos(angles)) ...
                    (-b1+obj.ri*cos(angles+pi/2)) ...
                    (-b1+obj.ri*cos(angles+pi)) ...
                    ( b1+obj.ri*cos(angles+1.5*pi)) b1+obj.ri];
                yi = [ ( d1+obj.ri*sin(angles)) ...
                    ( d1+obj.ri*sin(angles+pi/2)) ...
                    (-d1+obj.ri*sin(angles+pi)) ...
                    (-d1+obj.ri*sin(angles+1.5*pi)) d1 ];
            end
            fill(xo,yo,obj.color_steelFill,'LineStyle','none')
            fill(xi,yi,obj.color_concreteFill,'LineStyle','none')
            plot(xo,yo,'k-','LineWidth',lineWidth);
            plot(xi,yi,'k-','LineWidth',lineWidth);
            axis equal
        end
    end

    methods (Static)
        function type = memberType()
            type = 'rcft';
        end
        function tf = hasSteel()
            tf = true;
        end
        function tf = hasConcrete()
            tf = true;
        end
        function tf = hasReinforcement()
            tf = false;
        end
        function t = t_given_rho(rho,H,B,ri)
            if nargin < 4
                % ri = t;
                t0 = H*B*rho/(2*H+2*B);

                options = struct;
                options.Display = 'off';
                [t,~,exitflag] = fsolve(@(x)t_given_rho_err(H,B,x)-rho,t0,options);
                if exitflag <= 0
                    warning('fsolve could not find solution');
                    t = NaN;
                end
            else
                error('Not yet implemented')
            end
        end        
    end

end


function r = t_given_rho_err(H,B,t)
    Ag = H*B-(4-pi)*4*t^2;
    Ac = (H-2*t)*(B-2*t)-(4-pi)*t^2;
    As = Ag-Ac;
    r = As/Ag;
end