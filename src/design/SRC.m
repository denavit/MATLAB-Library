classdef SRC < structural_shape
    %SRC Summary of this class goes here
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
        % Gross Geometric Properties
        H               % Depth of the column
        B               % Width of the column
        % Steel Shape Geometric Properties
        d               % Depth of the steel section
        tw              % Thickness of the web
        bf              % Width of the flange
        tf              % Thickness of the flange
        % Steel Shape Material Properties
        Fy              % Yield stress of the steel shape
        Fu = [];        % Ultimate stress of the steel shape
        Es              % Modulus of Elasticity of Steel
        vs = 0.3        % Poisson's ratio of Steel
        % Reinforcement Geometric Properties
        db              % Diameter of the longitudinal reinforcing bars
        dbTies          % Diameter of the lateral reinforcing bars
        s               % Center-to-center spacing of the lateral reinforcing ties
        cover           % Clear cover to reinforcing steel
        % Reinforcement Material Properties
        Fylr            % Yield stress of the longitudinal reinforcing steel
        Fulr = [];      % Ultimate stress of the longitudinal reinforcing steel
        Eslr            % Modulus of Elasticity of the longitudinal reinforcing steel
        Fytr            % Yield stress of the transverse reinforcing steel
        % Concrete Material Properties
        fc              % Compressive strength of the concrete
        Ec              % Modulus of Elasticity of Concrete
        vc = 0.2        % Poisson's ratio of Concrete
        % Doubler Plate Properties
        t_doubler = [];
        Fy_doubler = [];
        % Information
        shapeName = ''; % Name of the steel shape
        % Clearance at the centerline of the column in each direction which
        % is to be kept clear of vertical bars so that a beam may frame to
        % the embedded rolled shape
        beamClearanceZ = 0;
        beamClearanceY = 0;
        % Design Options
        option_EI = 'AISC2016';
    end

    properties (Dependent)
        rebarConfig     % Configuration of the longitudinal reinforcing bars
    end

    properties (Access = private)
        nbZ
        nbY
        privateRebarConfig
    end

    methods
        %% Constructor
        function obj=SRC(d,tw,bf,tf,Fy,H,B,fc,db,config,Fylr,dbTies,s,Fytr,cover,units)
            obj.d           = d;
            obj.tw          = tw;
            obj.bf          = bf;
            obj.tf          = tf;
            obj.Fy          = Fy;
            obj.H           = H;
            obj.B           = B;
            obj.fc          = fc;
            obj.db          = db;
            obj.rebarConfig = config;
            obj.Fylr        = Fylr;
            obj.dbTies      = dbTies;
            obj.s           = s;
            obj.Fytr        = Fytr;
            obj.cover       = cover;
            obj.units       = units;

            switch obj.units
                case 'US'
                    obj.Es      = 29000;
                    obj.Eslr    = 29000;
                    obj.Ec      = 57*sqrt(obj.fc*1000);
                case 'SI'
                    obj.Es      = unitConvert('pressure',29000,'ksi','MPa');
                    obj.Eslr    = unitConvert('pressure',29000,'ksi','MPa');
                    fc_psi      = unitConvert('pressure',obj.fc,'MPa','psi');
                    Ec_ksi      = 57*sqrt(fc_psi);
                    obj.Ec      = unitConvert('pressure',Ec_ksi,'ksi','MPa');
                otherwise
                    warning('Design:Input','Unknown unit system, setting Es = 0 and Ec = 0')
                    obj.Es      = 0;
                    obj.Eslr    = 0;
                    obj.Ec      = 0;
            end
        end

        %% Set and Get Functions
        function config = get.rebarConfig(obj)
            config = obj.privateRebarConfig;
        end
        function obj = set.rebarConfig(obj,config)
            if strcmpi(config,'none')
                obj.privateRebarConfig = 'none';
                obj.nbZ = 0;
                obj.nbY = 0;
                return;
            end
            temp = sscanf(config,'%i%c-%i%c');
            assert(isequal(size(temp),[4 1]),...
                'invalid rebar configuration: %s',config);
            assert(temp(2)=='x' || temp(2)=='z',...
                'invalid rebar configuration: %s',config);
            assert(temp(4)=='y',...
                'invalid rebar configuration: %s',config);

            nbZt = temp(1);
            nbYt = temp(3);
            assert(nbZt > 0 && nbYt > 0 ,...
                'rebar should be specified in positive numbers');
            assert(mod(nbZt,2)==0 && mod(nbYt,2)==0 ,...
                'rebar should be specified in even numbers');

            obj.privateRebarConfig = config;
            obj.nbZ = nbZt/2;
            obj.nbY = nbYt/2;
        end
        function obj = set.db(obj,size)
            if ischar(size)
                obj.db = rebarSize(size);
            else
                obj.db = size;
            end
        end
        function obj = set.dbTies(obj,size)
            if ischar(size)
                obj.dbTies = rebarSize(size);
            else
                obj.dbTies = size;
            end
        end

        %% Geometric Properties
        function area = As(obj)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf);
            area = shape.A;
        end
        function area = Asr(obj)
            [~,~,A] = obj.rebarData;
            area = sum(A);
        end
        function area = Ac(obj)
            area = obj.H*obj.B - obj.As - obj.Asr;
        end
        function area = Ag(obj)
            area = obj.H*obj.B;
        end
        function inertia = Is(obj,axis)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf);
            inertia = shape.I(axis);
        end
        function inertia = Isr(obj,axis)
            switch lower(axis)
                case 'strong'
                    [~,y,A] = obj.rebarData;
                    inertia = sum(y.^2.*A);
                case 'weak'
                    [z,~,A] = obj.rebarData;
                    inertia = sum(z.^2.*A);
                otherwise
                    error('Bad axis');
            end
        end
        function inertia = Ic(obj,axis)
            switch lower(axis)
                case 'strong'
                    inertia = (1/12)*obj.B*obj.H^3 - ...
                        obj.Is('strong') - obj.Isr('strong');
                case 'weak'
                    inertia = (1/12)*obj.H*obj.B^3 - ...
                        obj.Is('weak') - obj.Isr('weak');
                otherwise
                    error('Bad axis');
            end
        end
        function inertia = Ig(obj,axis)
            switch lower(axis)
                case 'strong'
                    inertia = (1/12)*obj.B*obj.H^3;
                case 'weak'
                    inertia = (1/12)*obj.H*obj.B^3;
                otherwise
                    error('Bad axis');
            end            
        end
        function d = depth(obj,axis)
            switch lower(axis)
                case 'strong'
                    d = obj.H;
                case 'weak'
                    d = obj.B;
                otherwise
                    error('Bad axis');
            end
        end
        function [z,y,A] = rebarData(obj)
            if strcmpi(obj.rebarConfig,'none')
                z = zeros(0,1);
                y = zeros(0,1);
                A = zeros(0,1);
                return;
            end
            cornerZloc = obj.B/2 - obj.cover - obj.dbTies - obj.db/2;
            cornerYloc = obj.H/2 - obj.cover - obj.dbTies - obj.db/2;
            spacing = obj.minClearSpacing + obj.db;
            z1 = zeros((obj.nbZ + obj.nbY - 1),1);
            y1 = zeros((obj.nbZ + obj.nbY - 1),1);
            z1(1) = cornerZloc;
            y1(1) = cornerYloc;
            z1(2:(obj.nbZ)) = (cornerZloc-spacing):-spacing:(cornerZloc-(obj.nbZ-1)*spacing);
            y1(2:(obj.nbZ)) = cornerYloc;
            z1((obj.nbZ+1):end) = cornerZloc;
            y1((obj.nbZ+1):end) = (cornerYloc-spacing):-spacing:(cornerYloc-(obj.nbY-1)*spacing);

            z = vertcat(z1,-z1,-z1,z1);
            y = vertcat(y1,y1,-y1,-y1);

            numBars = 4*(obj.nbZ+obj.nbY-1);
            A = (pi/4)*obj.db^2*ones(numBars,1);
        end
        function spacing = minClearSpacing(obj)
            switch obj.units
                case 'US'
                    minSpacing = 1.5;
                case 'SI'
                    minSpacing = convertUnits.length(1.5,'in','mm');
                otherwise
                    error('Unknown unit system');
            end
            spacing = max([minSpacing 1.5*obj.db]);
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
            Pa = -(obj.As*obj.Fy + obj.Asr*obj.Fylr + 0.85*obj.fc*obj.Ac);
            Ma = 0;
        end
        function [Pb,Mb] = pointB(obj,axis)
            Pb = 0;
            Mb = compute_Mb(obj,axis);
        end
        function [Pc,Mc] = pointC(obj,axis)
            Pc = -(0.85*obj.fc*obj.Ac);
            Mc = compute_Mb(obj,axis);
        end
        function [Pd,Md] = pointD(obj,axis)
            PNA = 0;
            [P1,M1] = PSD_steel(obj,axis,PNA);
            [P2,M2] = PSD_conc(obj,axis,PNA);
            [P3,M3] = PSD_reinf(obj,axis,PNA);
            Pd = P1+P2+P3;
            Md = M1+M2+M3;
        end
        function [Pe,Me] = pointE(obj,axis)
            switch lower(axis)
                case 'strong'
                    PNA = -obj.d/2;
                case 'weak'
                    PNA = -obj.bf/2;
                otherwise
                    error('Unknown axis: %s',axis);
            end
            [P1,M1] = PSD_steel(obj,axis,PNA);
            [P2,M2] = PSD_conc(obj,axis,PNA);
            [P3,M3] = PSD_reinf(obj,axis,PNA);
            Pe = P1+P2+P3;
            Me = M1+M2+M3;
        end

        %% Design Strengths
        function pnco = Pnco(obj)
            % Stub Column (L=0) Strength
            [Pa,~] = obj.pointA;
            pnco = -Pa;
        end
        function C1 = C1(obj)
            switch obj.option_EI
                case 'AISC2010'
                    C1 = min(0.3,0.1+2*(obj.As/(obj.Ac+obj.As)));
                case 'AISC2016'
                    C1 = min(0.7,0.25+3*(obj.As+obj.Asr)/obj.Ag);
                otherwise
                    error('Unknown option_EI: %s',obj.option_EI);
            end
        end
        function eieff = EIeff(obj,axis)
            switch obj.option_EI
                case 'AISC2010'
                    eieff = obj.Es*obj.Is(axis) + 0.5*obj.Eslr*obj.Isr(axis) + ...
                        obj.C1*obj.Ec*obj.Ic(axis);
                case 'AISC2016'
                    eieff = obj.Es*obj.Is(axis) + obj.Eslr*obj.Isr(axis) + ...
                        obj.C1*obj.Ec*obj.Ic(axis);
                otherwise
                    error('Unknown option_EI: %s',obj.option_EI);
            end
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
            k = obj.tf; % neglect the fillet 
            steelSection = wfaisc2010(obj.d,obj.tw,obj.bf,obj.tf,k,obj.Fy,obj.units,obj.L,obj.Kstrong,obj.Kweak);
            steelSection.Fu = obj.Fu;
            steelSection.Es = obj.Es;
            steelSection.vs = obj.vs;
            pn = steelSection.Pnc(axis);
        end        
        function pnt = Pnt(obj)
            % Tension Strength
            pnt = obj.As*obj.Fy + obj.Asr*obj.Fylr;
        end
        function mno = Mno(obj,axis)
            % L=0 Moment Strength
            mno = obj.Mn(axis);
        end
        function mn = Mn(obj,axis)
            % Flexural Strength
            mn = compute_Mb(obj,axis);
        end
        function vn = Vn(obj,axis)
            % Shear Strength

            % Based on the steel section alone
            switch lower(axis)
                case 'strong'
                    axis = 'strong';
                    Aw = obj.d*obj.tw;
                    h = obj.d - 2*obj.tf;
                    lambda = h/obj.tw;
                    if ( lambda < 260 )
                        kv = 5;
                    else
                        warning('design:limit','web is too slender and should be stiffened')
                        vn = 0;
                        return
                    end
                case 'weak'
                    axis = 'weak';
                    Aw = 2*obj.bf*obj.tf;
                    lambda = obj.bf/2/obj.tf;
                    kv = 1.2;
                otherwise
                    error('Bad axis');
            end

            if ( strcmpi(axis,'strong') && lambda <= 2.24*sqrt(obj.Es/obj.Fy) )
                Cv = 1;
            else
                if ( lambda <= 1.10*sqrt(kv*obj.Es/obj.Fy) )
                    Cv = 1.0;
                elseif ( lambda <= 1.37*sqrt(kv*obj.Es/obj.Fy) )
                    Cv = 1.10*sqrt(kv*obj.Es/obj.Fy)/lambda;
                else
                    Cv = 1.51*obj.Es*kv/lambda^2/obj.Fy;
                end
            end
            vn = 0.6*obj.Fy*Aw*Cv;
        end
        function vpz = Vpz(obj,axis)
            switch lower(axis)
                case 'strong'

                    switch obj.units
                        % fcv[psi] = 15*sqrt(fc[psi])
                        case 'US'
                            fcv = 0.474*sqrt(obj.fc);
                        case 'SI'
                            fcv = 1.246*sqrt(obj.fc);
                        otherwise
                            error('Unknown unit system')
                    end
                    Vc  = fcv*obj.Ac;
                    Vs  = (0.6*obj.Fy)*obj.d*obj.tw;
                    if isempty(obj.t_doubler) || isempty(obj.Fy_doubler)
                        Vd = 0;
                    else
                        d_doubler = obj.d - 2*obj.tf;
                        Vd = (0.6*obj.Fy_doubler)*d_doubler*obj.t_doubler;
                    end
                    vpz = Vc+Vs+Vd;

                case 'weak'
                    error('Vpz not implemented for weak axis');
                otherwise
                    error('Bad axis');
            end
        end
        function kepz = Kepz(obj,axis,h,S)
            switch lower(axis)
                case 'strong'

                    tanAlpha = h/obj.H;
                    Acp = 5/6*(S/tanAlpha)*obj.B;
                    kec = obj.Gc*Acp;
                    kes = obj.Gs*obj.d*obj.tw;
                    if isempty(obj.t_doubler) || isempty(obj.Fy_doubler)
                        ked = 0;
                    else
                        d_doubler = obj.d - 2*obj.tf;
                        ked = obj.Gs*d_doubler*obj.t_doubler;
                    end
                    kepz = kec+kes+ked;

                case 'weak'
                    error('Kepz not implemented for weak axis');
                otherwise
                    error('Bad axis');
            end
        end
        function Kbs = bondSlip(obj,axis)          
            Kbs = (obj.Es*obj.Is(axis) + obj.Eslr*obj.Isr(axis) + 0.5*obj.Ec*obj.Ic(axis))/(4*obj.depth(axis));
        end

        %% Design Checks
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)

            % Resistance Factors
            phi_Pc = 0.75;
            phi_M  = 0.90;
            phi_Pt = 0.90;
            if ( ((obj.d-2*obj.tf)/obj.tw) <= 2.24*sqrt(obj.Es/obj.Fy) )
                phi_Vstrong  = 1.00;
            else
                phi_Vstrong  = 0.90;
            end
            phi_Vweak  = 0.90;

            % Compressive Load / Moment Interaction
            [P_pointCs,M_pointCs] = obj.pointC('strong');
            [~,M_pointCw] = obj.pointC('weak');
            x = obj.stabilityReduction('min',obj.Pnco);
            Pa = phi_Pc*x*obj.Pnco;
            Pc = phi_Pc*x*(-P_pointCs);
            Mcs = phi_M*M_pointCs;
            Mcw = phi_M*M_pointCw;
            ratio_PMc = aisc2010.interactionCheck_PSD(P,Ms,Mw,Pa,Pc,Mcs,Mcw);

            % Tensile Load / Moment Interaction
            Pc = phi_Pt*obj.Pnt;
            ratio_PMt = aisc2010.interactionCheck_H12(P,Ms,Mw,Pc,Mcs,Mcw);

            % Check strong axis shear
            if isempty(Vs)
                ratio_Vs = 0;
            else
                ratio_Vs = max(abs(Vs))/(phi_Vstrong*obj.Vn('strong'));
            end

            % Check weak axis shear
            if isempty(Vw)
                ratio_Vw = 0;
            else
                ratio_Vw = max(abs(Vw))/(phi_Vweak*obj.Vn('weak'));
            end

            % No check on torsion
            if ~isempty(T)
                error('Torsion strength check not impelemented');
            end

            ratio = max([ratio_PMc ratio_PMt ratio_Vs ratio_Vw]);
        end
        function pass_tf = proportioningCheck(obj,checkType,varargin)
            switch lower(checkType)
                case 'reinforcementspacing'
                    spacing = obj.minClearSpacing;
                    cornerZ = obj.B/2 - obj.cover - obj.dbTies - obj.db/2;
                    cornerY = obj.H/2 - obj.cover - obj.dbTies - obj.db/2;
                    lastBarZ = cornerZ - (obj.nbZ-1)*(spacing+obj.db);
                    lastBarY = cornerY - (obj.nbY-1)*(spacing+obj.db);
                    % Spacing to the steel shape (Y direction)
                    if obj.bf/2 >= lastBarZ
                        pass_tf1 = cornerY - obj.db/2 - obj.d/2 >= spacing;
                    else
                        pass_tf1 = (hypot(lastBarZ-obj.bf/2,cornerY-obj.d/2) >= spacing + obj.db/2) && ...
                            (obj.H/2 - obj.cover - obj.dbTies - obj.d/2 >= spacing);
                    end
                    % Spacing to the steel shape (Z direction)
                    if obj.d/2 >= lastBarY
                        pass_tf2 = cornerZ - obj.db/2 - obj.bf/2 >= spacing;
                    else
                        pass_tf2 = (hypot(cornerZ-obj.bf/2,lastBarY-obj.d/2) >= spacing + obj.db/2) && ...
                            (obj.B/2 - obj.cover - obj.dbTies - obj.bf/2 >= spacing);
                    end
                    % Beam clearance (Y direction)
                    pass_tf3 = lastBarY-obj.db/2-obj.beamClearanceY/2 >= spacing;
                    % Beam clearance (Z direction)
                    pass_tf4 = lastBarZ-obj.db/2-obj.beamClearanceZ/2 >= spacing;
                    % All checks must pass
                    pass_tf = pass_tf1 && pass_tf2 && pass_tf3 && pass_tf4;
                case 'steelratio'
                    pass_tf = 0.01 <= obj.As/(obj.H*obj.B);
                case 'reinforcementratio'
                    pass_tf = 0.004 <= obj.Asr/(obj.H*obj.B);
                case 'slendernessratio_klr'
                    limit = varargin{1};
                    KLr_strong = obj.Kstrong*obj.L/sqrt(obj.Is('strong')/obj.As);
                    KLr_weak = obj.Kweak*obj.L/sqrt(obj.Is('weak')/obj.As);
                    pass_tf = max([KLr_strong KLr_weak]) < limit;
                case 'moderatelyductilesection'
                    switch obj.units
                        case 'US'
                            maxSpacingTransverseReinf4 = 12;
                        case 'SI'
                            maxSpacingTransverseReinf4 = 300;
                        otherwise
                            error('Unknown unit system');
                    end

                    % AISC 341-10, D1.4b(1)(1)
                    maxSpacingTransverseReinf = min([
                        0.5*min([obj.H obj.B])
                        8*obj.db
                        24*obj.dbTies
                        maxSpacingTransverseReinf4]);
                    pass_tf = obj.s <= maxSpacingTransverseReinf;

                    % AISC 341-10, D1.4b(1)(2) & D1.4b(1)(3)
                    % Spacing assumed to be constant long entire length

                    % AISC 341-10, D1.4b(1)(4)
                    % Any splices and end bearing details are assumed to be
                    % satisfied

                    % AISC 341-10, D1.4b(1)(5)
                    % No welded wire fabric included
                case 'highlyductilesection'
                    switch obj.units
                        case 'US'
                            maxSpacingTransverseReinf_2ii = 6;
                        case 'SI'
                            maxSpacingTransverseReinf_2ii = 150;
                        otherwise
                            error('Unknown unit system');
                    end

                    % AISC 341-10, D1.4b(2)
                    pass_tf_md = obj.proportioningCheck('ModeratelyDuctileSection');

                    % AISC 341-10, D1.4b(2)(1) -> ACI 318-08, 21.6.3
                    pass_tf_1 = (obj.Asr >= 0.01*obj.H*obj.B) && (obj.Asr <= 0.06*obj.H*obj.B);

                    % AISC 341-10, D1.4b(2)(2)(i)
                    hcc = max([obj.H obj.B]) - 2*obj.cover;
                    Ash_min = 0.09*hcc*obj.s*(1-obj.Fy*obj.As/obj.Pnco)*(obj.fc/obj.Fytr);
                    pass_tf_2i = 2*(pi/4)*obj.dbTies^2 >= Ash_min;

                    % AISC 341-10, D1.4b(2)(2)(ii)
                    maxSpacingTransverseReinf = min([
                        6*obj.db
                        maxSpacingTransverseReinf_2ii]);
                    pass_tf_2ii = obj.s <= maxSpacingTransverseReinf;

                    % AISC 341-10, D1.4b(2)(3) - D1.4b(2)(6)
                    % Not checked

                    pass_tf = pass_tf_md && pass_tf_1 && pass_tf_2i && pass_tf_2ii;
                otherwise
                    error('Unknown checkType: %s',checkType);
            end
        end

        %% Interaction Strength
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case 'aisc'
                    [P,M] = obj.sectionInteraction2d(axis,type,'psd-acbt');
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
                    phi = ACI_phi('other',et,obj.Fy/obj.Es);
                    P = phi.*P;
                    M = phi.*M;                     
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case 'aisc'
                    [P,M] = obj.beamColumnInteraction2d(axis,'psdsimple-acbt',quadrant);
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
                case 'trial'
                    [P,M] = trial_interaction_diagram(obj,axis,quadrant,type(7:end),false);
                case 'factoredtrial'
                    [P,M] = trial_interaction_diagram(obj,axis,quadrant,type(15:end),true);
                otherwise
                    error('Unknown type: %s',type);
            end
        end

        %% Export and Information Functions
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(obj,axis,type)
            switch lower(type)
                case 'columnstrength'
                    EAeff = obj.Es*obj.As + obj.Eslr*obj.Asr + obj.Ec*obj.Ac;
                    E = obj.Es;
                    I = obj.EIeff(axis)/E;
                    A = EAeff/E;
                case 'gross'
                    EAeff = obj.Es*obj.As + obj.Eslr*obj.Asr + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is(axis) + obj.Eslr*obj.Isr(axis) + ...
                        obj.Ec*obj.Ic(axis);
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                case 'steelonly'
                    E = obj.Es;
                    I = obj.Is(axis);
                    A = obj.As;
                case 'proposed'
                    EAeff = obj.Es*obj.As + obj.Eslr*obj.Asr + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is(axis) + obj.Eslr*obj.Isr(axis) + ...
                        obj.C1*obj.Ec*obj.Ic(axis);
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;
                case 'aci'
                    E = obj.Ec;
                    I = 0.7*obj.Ig(axis);
                    A = obj.Ag;
                case 'servicedefl'
                    EAeff = obj.Es*obj.As + obj.Eslr*obj.Asr + obj.Ec*obj.Ac;
                    EIeff = obj.Es*obj.Is(axis) + obj.Eslr*obj.Isr(axis) + ...
                        0.4*obj.Ec*obj.Ic(axis);
                    E = obj.Es;
                    I = EIeff/E;
                    A = EAeff/E;                    
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [E,A,Iz,Iy,GJ] = sectionPropertiesForElasticAnalysis3d(obj,type)
            switch lower(type)
                case 'columnstrength'
                    E   = obj.Es;
                    A   = (obj.Es*obj.As + obj.Eslr*obj.Asr + obj.Ec*obj.Ac)/E;
                    Iz  = obj.EIeff('z')/E;
                    Iy  = obj.EIeff('y')/E;
                    GJ  = min(obj.Gs*obj.Js,obj.Gc*obj.Jc);
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function lp = Lp(obj,axis,Li)
            switch lower(axis)
                case 'strong'
                    lp = 0.12*Li;
                case 'weak'
                    lp = 0.20*Li;
                otherwise
                    error('Bad axis');
            end
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');
            switch lower(axis)
                case 'strong'
                    yExtreme = obj.H/2;
                    yExtremeConcrete = obj.H/2;
                case 'weak'
                    yExtreme = obj.B/2;
                    yExtremeConcrete = obj.B/2;
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
            
            z =  obj.B/2;
            y =  obj.H/2;
            strain_c1 = axialStrain + z*curvatureY + y*curvatureZ;
            strain_c2 = axialStrain + z*curvatureY - y*curvatureZ;
            strain_c3 = axialStrain - z*curvatureY - y*curvatureZ;
            strain_c4 = axialStrain - z*curvatureY + y*curvatureZ;
            
            switch lower(type)
                case 'maxcompressive'
                    strain = min([strain_c1(:) strain_c2(:) strain_c3(:) strain_c4(:)],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxtensile'
                    strain = max([strain_c1(:) strain_c2(:) strain_c3(:) strain_c4(:)],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxabsolute'
                    strain = max(abs([strain_c1(:) strain_c2(:) strain_c3(:) strain_c4(:)]),[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxconcretecompressive'
                    strain = min([strain_c1(:) strain_c2(:) strain_c3(:) strain_c4(:)],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxconcretetensile'
                    strain = max([strain_c1(:) strain_c2(:) strain_c3(:) strain_c4(:)],[],2);
                    strain = reshape(strain,size(axialStrain));
                otherwise
                    error('Unknown type: %s',type);
            end
        end               
        function x = getSectionData(obj,type,axis)
            switch lower(type)
                case 'steelarea'
                    x = obj.As;
                case 'reinforcingarea'
                    x = obj.Asr;
                case 'concretearea'
                    x = obj.Ac;
                case 'grossarea'
                    x = obj.Ag;
                case 'steelratio'
                    x = obj.As/obj.Ag;
                case 'reinforcingratio'
                    x = obj.Asr/obj.Ag;
                case 'steelstrength'
                    x = obj.Fy;
                case 'reinforcingstrength'
                    x = obj.Fylr;
                case 'concretestrength'
                    x = obj.fc;
                case 'grossconcreteflexuralrigidity'
                    x = obj.Ec*obj.Ic(axis);
                case 'grosssteelflexuralrigidity'
                    x = obj.Es*obj.Is(axis);
                case 'grossreinforcingflexuralrigidity'
                    x = obj.Eslr*obj.Isr(axis);
                case 'ecig'
                    x = obj.Ec*obj.Ig(axis);
                case 'grosssectioncompressionstrength'
                    x = obj.As*obj.Fy + obj.Asr*obj.Fylr + obj.fc*obj.Ac;
                case {'asfy','grosssteelstrength'}
                    x = obj.As*obj.Fy;
                case {'asrfysr','grossreinforcingstrength'}
                    x = obj.Asr*obj.Fylr;
                case {'acfc','grossconcretestrength'}
                    x = obj.Ac*obj.fc;
                case 'aci-pmax'
                    x = 0.85*(0.85*obj.fc*obj.Ac + obj.As*obj.Fy + obj.Asr*obj.Fylr);
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
                    description = sprintf('SRC %gx%g %s',...
                        obj.H,obj.B,obj.shapeName);
                otherwise
                    error('Unknown flag');
            end
        end
        function fs = fiberSectionObject(obj,idSteel,idConc,idReinf)
            [zReinf,yReinf,~] = obj.rebarData;
            fs = fiberSection;
            fs.addSRC(idSteel,idConc,idReinf,...
                obj.H,obj.B,obj.d,obj.tw,obj.bf,obj.tf,...
                zReinf,yReinf,obj.db);
        end
        function psd = plasticStressDistributionObject(obj)
            idSteel = 1;
            idConc  = 2;
            idReinf = 3;
            fs = obj.fiberSectionObject(idSteel,idConc,idReinf);
            psd = plastic_stress_distribution(fs);
            psd.addMaterial(idSteel,obj.Fy,-obj.Fy);
            psd.addMaterial(idConc,0.0,-0.85*obj.fc);
            psd.addMaterial(idReinf,obj.Fylr,-obj.Fylr);
        end
        function ess = elasticSectionStiffnessObject(obj)
            idSteel = 1;
            idConc  = 2;
            idReinf = 3;
            fs = obj.fiberSectionObject(idSteel,idConc,idReinf);
            ess = elastic_section_stiffness(fs);
            ess.addMaterial(idSteel,obj.Es);
            ess.addMaterial(idConc ,obj.Ec,'CompressionOnly');
            ess.addMaterial(idReinf,obj.Eslr);
        end
        function scACI = strainCompatibilityAciObject(obj)
            idSteel = 1;
            idConc  = 2;
            idReinf = 3;
            fs = obj.fiberSectionObject(idSteel,idConc,idReinf);
            scACI = ACI_strain_compatibility(fs);
            scACI.addConcreteBoundary(-obj.B/2,-obj.H/2,0);
            scACI.addConcreteBoundary(-obj.B/2, obj.H/2,0);
            scACI.addConcreteBoundary( obj.B/2, obj.H/2,0);
            scACI.addConcreteBoundary( obj.B/2,-obj.H/2,0);
            scACI.addSteelBoundary(-obj.bf/2,-obj.d/2,0);
            scACI.addSteelBoundary(-obj.bf/2, obj.d/2,0);
            scACI.addSteelBoundary( obj.bf/2, obj.d/2,0);
            scACI.addSteelBoundary( obj.bf/2,-obj.d/2,0);
            [z,y,A] = obj.rebarData;
            for i = 1:length(z)
                scACI.addSteelBoundary(z(i),y(i),sqrt(A(i)/pi)); % Assuming circular reinforcement
            end
            scACI.maxCompressiveStrength = -0.85*(obj.As*obj.Fy + obj.Asr*obj.Fylr + 0.85*obj.Ac*obj.fc);
            scACI.addMaterial('steel',idSteel,obj.Fy,obj.Es);
            scACI.addMaterial('concrete',idConc,obj.fc,obj.units);
            scACI.addMaterial('steel',idReinf,obj.Fylr,obj.Eslr);
        end         
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            % Gross Section
            x = [-0.5 0.5 0.5 -0.5 -0.5]*obj.B;
            y = [0.5 0.5 -0.5 -0.5 0.5]*obj.H;
            fill(x,y,obj.color_concreteFill,'LineStyle','none')
            plot(x,y,'k-','LineWidth',lineWidth);
            % Steel Shape
            x = [-0.5*obj.bf 0.5*obj.bf 0.5*obj.bf 0.5*obj.tw 0.5*obj.tw ...
                0.5*obj.bf 0.5*obj.bf -0.5*obj.bf -0.5*obj.bf -0.5*obj.tw ...
                -0.5*obj.tw -0.5*obj.bf -0.5*obj.bf];
            y = [0.5*obj.d 0.5*obj.d 0.5*obj.d-obj.tf 0.5*obj.d-obj.tf ...
                -0.5*obj.d+obj.tf -0.5*obj.d+obj.tf -0.5*obj.d -0.5*obj.d ...
                -0.5*obj.d+obj.tf -0.5*obj.d+obj.tf 0.5*obj.d-obj.tf ...
                0.5*obj.d-obj.tf 0.5*obj.d];
            fill(x,y,obj.color_steelFill,'LineStyle','none')
            plot(x,y,'k-','LineWidth',lineWidth);
            % Transverse Rebar
            angles = linspace(0,pi/2,25);
            riBend = 2*obj.dbTies;
            roBend = riBend + obj.dbTies;
            d1 = obj.H/2 - obj.cover - obj.dbTies - riBend;
            b1 = obj.B/2 - obj.cover - obj.dbTies - riBend;
            xo = [ ( b1+roBend*cos(angles)) ...
                   (-b1+roBend*cos(angles+pi/2)) ...
                   (-b1+roBend*cos(angles+pi)) ...
                   ( b1+roBend*cos(angles+1.5*pi)) b1+roBend ];
            yo = [ ( d1+roBend*sin(angles)) ...
                   ( d1+roBend*sin(angles+pi/2)) ...
                   (-d1+roBend*sin(angles+pi)) ...
                   (-d1+roBend*sin(angles+1.5*pi)) d1 ];
            xi = [ ( b1+riBend*cos(angles)) ...
                   (-b1+riBend*cos(angles+pi/2)) ...
                   (-b1+riBend*cos(angles+pi)) ...
                   ( b1+riBend*cos(angles+1.5*pi)) b1+riBend ];
            yi = [ ( d1+riBend*sin(angles)) ...
                   ( d1+riBend*sin(angles+pi/2)) ...
                   (-d1+riBend*sin(angles+pi)) ...
                   (-d1+riBend*sin(angles+1.5*pi)) d1 ];
            fill([xi fliplr(xo) xi(1)],[yi fliplr(yo) yi(1)],...
                obj.color_steelFill,'LineStyle','none')
            plot(xi,yi,'k-','LineWidth',lineWidth);
            plot(xo,yo,'k-','LineWidth',lineWidth);
            % Longitudinal Rebar
            angles = linspace(0,2*pi,25);
            circX = (obj.db/2*cos(angles));
            circY = (obj.db/2*sin(angles));
            [rebarZ,rebarY,~] = obj.rebarData;
            for i = 1:length(rebarZ)
                x =  rebarZ(i) + circX;
                y =  rebarY(i) + circY;
                fill(x,y,obj.color_steelFill,'LineStyle','none')
                plot(x,y,'k-','LineWidth',lineWidth);
            end
            axis equal
        end
        
        function tf = hasReinforcement(obj)
            tf = ~strcmpi(obj.rebarConfig,'none');
        end
    end

    methods (Static)
        function type = memberType()
            type = 'src';
        end
        function tf = hasConcrete()
            tf = true;
        end
    end

end


function Mb = compute_Mb(obj,axis)

switch lower(axis)
    case 'strong'
        [~,ybar,~] = obj.rebarData;
        PNAb = vertcat(ybar(ybar>0),0,obj.H/2,obj.d/2,obj.d/2-obj.tf);
    case 'weak'
        [ybar,~,~] = obj.rebarData;
        PNAb = vertcat(ybar(ybar>0),0,obj.B/2,obj.bf/2,obj.tw/2);
    otherwise
        error('Bad axis');
end

PNAb = sort(unique(PNAb));

for i = 2:length(PNAb)
    PNAi = (PNAb(i)+PNAb(i-1))/2;
    [P1,~,dPdPNA1] = PSD_steel(obj,axis,PNAi);
    [P2,~,dPdPNA2] = PSD_conc(obj,axis,PNAi);
    [P3,~] = PSD_reinf(obj,axis,PNAi);
    P = P1+P2+P3;
    dPdPNA = dPdPNA1+dPdPNA2;
    PNAc = PNAi - P/dPdPNA;
    if PNAb(i-1) <= PNAc && PNAc <= PNAb(i)
        [~,M1] = PSD_steel(obj,axis,PNAc);
        [~,M2] = PSD_conc(obj,axis,PNAc);
        [~,M3] = PSD_reinf(obj,axis,PNAc);
        Mb = M1+M2+M3;
        return
    end
end
ySmall = 1e-6;
for i = 2:length(PNAb)
    [P11,~] = PSD_steel(obj,axis,PNAb(i)-ySmall);
    [P21,~] = PSD_conc(obj,axis,PNAb(i)-ySmall);
    [P31,~] = PSD_reinf(obj,axis,PNAb(i)-ySmall);
    P1 = P11+P21+P31;
    [P12,~] = PSD_steel(obj,axis,PNAb(i)+ySmall);
    [P22,~] = PSD_conc(obj,axis,PNAb(i)+ySmall);
    [P32,~] = PSD_reinf(obj,axis,PNAb(i)+ySmall);
    P2 = P12+P22+P32;
    if P1 <= 0 && P2 >= 0
        [~,M1] = PSD_steel(obj,axis,PNAb(i));
        [~,M2] = PSD_conc(obj,axis,PNAb(i));
        [~,M3] = PSD_reinf(obj,axis,PNAb(i));
        Mb = M1+M2+M3;
        return
    end
end
error('Mb not found')
end

function [P,M,dPdPNA] = PSD_steel(obj,axis,PNA)

switch lower(axis)
    case 'strong'
        if PNA >= obj.d/2
            P = obj.As*obj.Fy;
            M = 0;
            dPdPNA = 0;
        elseif PNA >= obj.d/2 - obj.tf
            Ft = ((obj.d/2-PNA)*obj.bf)*obj.Fy;
            P = obj.As*obj.Fy - 2*Ft;
            M = 2*Ft*(obj.d/2-(obj.d/2-PNA)/2);
            dPdPNA = 2*obj.bf*obj.Fy;
        elseif PNA >= -obj.d/2 + obj.tf
            P = 2*PNA*obj.tw*obj.Fy;
            temp = obj.d/2 - obj.tf - PNA;
            M = 2*(obj.bf*obj.tf)*obj.Fy*(obj.d/2-obj.tf/2) + ...
                2*(obj.tw*temp)*obj.Fy*(obj.d/2-obj.tf-temp/2);
            dPdPNA = 2*obj.tw*obj.Fy;
        elseif PNA >= -obj.d/2
            Fc = ((obj.d/2+PNA)*obj.bf)*obj.Fy;
            P = -obj.As*obj.Fy + 2*Fc;
            M = 2*Fc*(obj.d/2-(obj.d/2+PNA)/2);
            dPdPNA = 2*obj.bf*obj.Fy;
        else
            P = -obj.As*obj.Fy;
            M = 0;
            dPdPNA = 0;
        end
    case 'weak'
        if PNA >= obj.bf/2
            P = obj.As*obj.Fy;
            M = 0;
            dPdPNA = 0;
        elseif PNA >= obj.tw/2
            Ft = ((obj.bf/2-PNA)*2*obj.tf)*obj.Fy;
            P = obj.As*obj.Fy - 2*Ft;
            M = 2*Ft*(obj.bf/2-(obj.bf/2-PNA)/2);
            dPdPNA = 2*2*obj.tf*obj.Fy;
        elseif PNA >= -obj.tw/2
            P = 2*PNA*obj.d*obj.Fy;
            temp = obj.tw/2 - PNA;
            M = 2*(2*obj.tf*(obj.bf/2-obj.tw/2))*obj.Fy*(obj.bf/4+obj.tw/4) + ...
                2*(obj.d*temp)*obj.Fy*(obj.tw/2-temp/2);
            dPdPNA = 2*obj.d*obj.Fy;
        elseif PNA >= -obj.bf/2
            Fc = ((obj.bf/2+PNA)*2*obj.tf)*obj.Fy;
            P = -obj.As*obj.Fy + 2*Fc;
            M = 2*Fc*(obj.bf/2-(obj.bf/2+PNA)/2);
            dPdPNA = 2*2*obj.tf*obj.Fy;
        else
            P = -obj.As*obj.Fy;
            M = 0;
            dPdPNA = 0;
        end
    otherwise
        error('Bad axis');
end
if nargout < 3
    clear dPdPNA
end

end

function [P,M,dPdPNA] = PSD_conc(obj,axis,PNA)

C2 = 0.85;
switch lower(axis)
    case 'strong'
        [~,ybar,A] = obj.rebarData;
        if PNA >= obj.H/2
            P = 0;
            M = 0;
            dPdPNA = 0;
        elseif PNA >= obj.d/2
            temp = obj.H/2-PNA;
            P = obj.B*temp*-C2*obj.fc;
            M = -P*(obj.H/2 - temp/2);
            dPdPNA = obj.B*C2*obj.fc;
        elseif PNA >= obj.d/2 - obj.tf

            P1 = obj.B*(obj.H/2-obj.d/2)*-C2*obj.fc;
            M1 = -P1*(obj.H/4 + obj.d/4);

            temp = obj.d/2 - PNA;
            P2 = (obj.B-obj.bf)*temp*-C2*obj.fc;
            M2 = -P2*(obj.d/2 - temp/2);

            P = P1+P2;
            M = M1+M2;
            dPdPNA = (obj.B-obj.bf)*C2*obj.fc;
        elseif PNA >= -obj.d/2 + obj.tf
            P1 = obj.B*(obj.H/2-obj.d/2)*-C2*obj.fc;
            M1 = -P1*(obj.H/4 + obj.d/4);
            P2 = (obj.B-obj.bf)*obj.tf*-C2*obj.fc;
            M2 = -P2*(obj.d/2 - obj.tf/2);

            temp = obj.d/2 -obj.tf - PNA;
            P3 = (obj.B-obj.tw)*temp*-C2*obj.fc;
            M3 = -P3*(obj.d/2 - obj.tf - temp/2);

            P = P1+P2+P3;
            M = M1+M2+M3;
            dPdPNA = (obj.B-obj.tw)*C2*obj.fc;

        elseif PNA >= -obj.d/2

            P1 = obj.B*(obj.H/2-obj.d/2)*-C2*obj.fc;
            M1 = -P1*(obj.H/4 + obj.d/4);
            P2 = (obj.B-obj.bf)*obj.tf*-C2*obj.fc;
            M2 = -P2*(obj.d/2 - obj.tf/2);
            P3 = (obj.B-obj.tw)*(obj.d-2*obj.tf)*-C2*obj.fc;

            temp = -obj.d/2 + obj.tf - PNA;
            P4 = (obj.B-obj.bf)*temp*-C2*obj.fc;
            M4 = -P4*( -obj.d/2 + obj.tf - temp/2);

            P = P1+P2+P3+P4;
            M = M1+M2+M4;
            dPdPNA = (obj.B-obj.bf)*C2*obj.fc;

        elseif PNA >= -obj.H/2

            P1 = obj.B*(obj.H/2-obj.d/2)*-C2*obj.fc;
            M1 = -P1*(obj.H/4 + obj.d/4);
            P2 = (obj.B-obj.bf)*obj.tf*-C2*obj.fc;
            P3 = (obj.B-obj.tw)*(obj.d-2*obj.tf)*-C2*obj.fc;
            P4 = (obj.B-obj.bf)*obj.tf*-C2*obj.fc;

            temp = -obj.d/2 - PNA;
            P5 = obj.B*temp*-C2*obj.fc;
            M5 = -P5*( -obj.d/2 - temp/2);

            P = P1+P2+P3+P4+P5;
            M = M1+M5;
            dPdPNA = obj.B*C2*obj.fc;

        else
            P = -0.85*obj.Ac*obj.fc;
            M = 0;
            dPdPNA = 0;
        end
    case 'weak'
        [ybar,~,A] = obj.rebarData;
        if PNA >= obj.B/2
            P = 0;
            M = 0;
            dPdPNA = 0;
        elseif PNA >= obj.bf/2
            temp = obj.B/2-PNA;
            P = obj.H*temp*-C2*obj.fc;
            M = -P*(obj.B/2 - temp/2);
            dPdPNA = obj.H*C2*obj.fc;
        elseif PNA >= obj.tw/2

            P1 = obj.H*(obj.B/2-obj.bf/2)*-C2*obj.fc;
            M1 = -P1*(obj.B/4 + obj.bf/4);

            temp = obj.bf/2 - PNA;
            P2 = (obj.H-2*obj.tf)*temp*-C2*obj.fc;
            M2 = -P2*(obj.bf/2 - temp/2);

            P = P1+P2;
            M = M1+M2;
            dPdPNA = (obj.H-2*obj.tf)*C2*obj.fc;
        elseif PNA >= -obj.tw/2
            P1 = obj.H*(obj.B/2-obj.bf/2)*-C2*obj.fc;
            M1 = -P1*(obj.B/4 + obj.bf/4);
            P2 = (obj.H-2*obj.tf)*(obj.bf/2-obj.tw/2)*-C2*obj.fc;
            M2 = -P2*(obj.bf/4 + obj.tw/4);

            temp = obj.tw/2 - PNA;
            P3 = (obj.H-obj.d)*temp*-C2*obj.fc;
            M3 = -P3*(obj.tw/2 - temp/2);

            P = P1+P2+P3;
            M = M1+M2+M3;
            dPdPNA = (obj.H-obj.d)*C2*obj.fc;

        elseif PNA >= -obj.bf/2

            P1 = obj.H*(obj.B/2-obj.bf/2)*-C2*obj.fc;
            M1 = -P1*(obj.B/4 + obj.bf/4);
            P2 = (obj.H-2*obj.tf)*(obj.bf/2-obj.tw/2)*-C2*obj.fc;
            M2 = -P2*(obj.bf/4 + obj.tw/4);
            P3 = (obj.H-obj.d)*obj.tw*-C2*obj.fc;

            temp = -obj.tw/2 - PNA;
            P4 = (obj.H-2*obj.tf)*temp*-C2*obj.fc;
            M4 = -P4*( -obj.tw/2 - temp/2);

            P = P1+P2+P3+P4;
            M = M1+M2+M4;
            dPdPNA = (obj.H-2*obj.tf)*C2*obj.fc;

        elseif PNA >= -obj.B/2

            P1 = obj.H*(obj.B/2-obj.bf/2)*-C2*obj.fc;
            M1 = -P1*(obj.B/4 + obj.bf/4);
            P2 = (obj.H-2*obj.tf)*(obj.bf/2-obj.tw/2)*-C2*obj.fc;
            P3 = (obj.H-obj.d)*obj.tw*-C2*obj.fc;
            P4 = (obj.H-2*obj.tf)*(obj.bf/2-obj.tw/2)*-C2*obj.fc;

            temp = -obj.bf/2 - PNA;
            P5 = obj.H*temp*-C2*obj.fc;
            M5 = -P5*( -obj.bf/2 - temp/2);

            P = P1+P2+P3+P4+P5;
            M = M1+M5;
            dPdPNA = obj.H*C2*obj.fc;

        else
            P = -0.85*obj.Ac*obj.fc;
            M = 0;
            dPdPNA = 0;
        end
    otherwise
        error('Bad axis');
end

F = zeros(size(A));
ind = ybar > PNA;
F(ind) = C2*obj.fc*A(ind);

P = P + sum(F);
M = M + sum(-F.*ybar);

if nargout < 3
    clear dPdPNA
end

end

function [P,M] = PSD_reinf(obj,axis,PNA)

switch lower(axis)
    case 'strong'
        [~,ybar,A] = obj.rebarData;
    case 'weak'
        [ybar,~,A] = obj.rebarData;
    otherwise
        error('Bad axis');
end

F = zeros(size(A));
ind = ybar < PNA;
F(ind) =  obj.Fylr*A(ind);
ind = ybar > PNA;
F(ind) = -obj.Fylr*A(ind);

P = sum(F);
M = sum(-F.*ybar);

end

function db = rebarSize(str)
    num = sscanf(str,'#%i');
    assert(numel(num)==1,'invalid bar size: %s',str)
    rebarData = [
        3	0.375
        4	0.5
        5	0.625
        6	0.75
        7	0.875
        8	1
        9	1.128
        10	1.27
        11	1.41
        14	1.693
        18	2.257];
    ind = find(rebarData(:,1)==num);
    assert(~isempty(ind),'invalid bar size: %i',num)
    db = rebarData(ind,2);
end
