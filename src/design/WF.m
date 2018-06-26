classdef WF < structural_shape
    % WF 
    %
    % Notes:
    %
    % REFERENCES:
    % - American Institute of Steel Construction (AISC) (2010). Specification 
    % for Structural Steel Buildings, ANSI/AISC 360-10, American Institute 
    % of Steel Construction, Chicago, Illinois.
    % 
    
    properties
        % Geometric Properties
        d           % Depth of the steel section 
        tw          % Thickness of the web 
        bf          % Width of the flange
        tf          % Thickness of the flange
        k           % Distance from foot of fillet to outside edge of flange
        % Material Properties
        Fy          % Steel yield stress
        Fu = [];    % Steel Ultimate Stress
        eu = [];    % Strain at Ultimate Stress
        Es          % Modulus of Elasticity of Steel
        G           % Shear Modulus of Elasticity of Steel
        % Expected Strength Parameters
        Ry = [];
        Rt = [];
        % Information
        shapeName = ''; % Name of the steel shape
        % Design Options
        Lb = [];
        Cb = [];
        neglectLocalBuckling = false;
        neglectLateralTorsionalBuckling = false;
    end
    
    methods
        %% Constructor
        function obj=WF(d,tw,bf,tf,k,Fy,units)
            obj.d  = d;
            obj.tw = tw;
            obj.bf = bf;
            obj.tf = tf;
            obj.k  = k;
            obj.Fy = Fy;
            obj.units = units;

            switch obj.units
                case 'US'
                    obj.Es = 29000;
                    obj.G  = 11200;
                case 'SI'
                    obj.Es = convertUnits.pressure(29000,'ksi','MPa');
                    obj.G  = convertUnits.pressure(11200,'ksi','MPa');
                otherwise
                    warning('WF:badInput','Unknown unit system, setting Es = 0')
                    obj.Es = 0;
                    obj.G  = 0;
            end
            
        end     
        
        %% Geometric Properties
        function area = A(obj)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf,obj.k);
            area = shape.A;
        end
        function j = J(obj)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf,obj.k);
            j = shape.J;
        end        
        function inertia = I(obj,axis)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf,obj.k);
            inertia = shape.I(axis);       
        end   
        function z = Z(obj,axis)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf,obj.k);
            z = shape.Z(axis);              
        end
        function s = S(obj,axis)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf,obj.k);
            s = shape.S(axis);
        end
        function radius = r(obj,axis)
            shape = I_Shape(obj.d,obj.tw,obj.bf,obj.tf,obj.k);
            radius = shape.r(axis);
        end
        function h = hw(obj)
            h = obj.d - 2*obj.k;
        end
        function d = depth(obj,axis)
            switch lower(axis)
                case 'strong'
                    d = obj.d;
                case 'weak'
                    d = obj.bf;
                otherwise
                    error('Bad axis'); 
            end            
        end
        
        %% Design Strengths
        function lambda = lambda(obj,axis,component)
            if strcmpi(component,'flange')
                lambda = (obj.bf/2)/obj.tf;
            elseif strcmpi(component,'web')
                lambda = obj.hw/obj.tw;
            else
                error('Unknown axis: %s and/or component: %s',axis,component)
            end
        end
        function pnco = Pnco(obj)
            % Stub Column (L=0) Strength
            
            if obj.neglectLocalBuckling
                pnco = obj.Fy*obj.A;
                return
            end
            
            lambdaf = obj.lambda('strong','flange');
            lambdaRf = 0.56*sqrt(obj.Es/obj.Fy);
            if ( lambdaf <= lambdaRf )
                flangeSlenderness = 'noncompact';
            else 
                flangeSlenderness = 'slender';
            end
            
            lambdaw = obj.lambda('strong','web');
            lambdaRw = 1.49*sqrt(obj.Es/obj.Fy);
            if ( lambdaw <= lambdaRw )
                webSlenderness = 'noncompact';
            else 
                webSlenderness = 'slender';
            end
            
            if strcmp(flangeSlenderness,'noncompact') && strcmp(webSlenderness,'noncompact')    
                pnco = obj.Fy*obj.A;               
            else
                if strcmp(flangeSlenderness,'slender')
                    kc = 4/sqrt(lambdaw);
                    if (kc > 0.76)
                        if( lambdaf <= 1.03*sqrt(obj.E/obj.Fy) )
                            Qs = 1.415 - 0.74*lambdaf*sqrt(obj.Fy/obj.E);
                        else
                            Qs = 0.69*obj.E/(obj.Fy*lambdaf^2);
                        end
                    else
                        if kc < 0.35; kc = 0.35; end
                        if( lambdaf <= 0.64*sqrt(obj.E*kc/obj.Fy) )
                            Qs = 1;
                        elseif( lambdaf <= 0.64*sqrt(obj.E*kc/obj.Fy) )
                            Qs = 1.415 - 0.65*lambdaf*sqrt(obj.Fy/obj.E/kc);
                        else
                            Qs = 0.90*obj.E*kc/(obj.Fy*lambdaf^2);
                        end                        
                    end
                else 
                    Qs = 1;
                end
                if strcmpi(webSlenderness,'slender')
                    Fe = min([...
                        pi^2*obj.Es/(obj.K('strong')*obj.L('strong')/obj.r('strong'))^2;
                        pi^2*obj.Es/(obj.K('weak')*obj.L('weak')/obj.r('weak'))^2]);
                    f = AISC_column_curve(obj.Fy/Fe)*obj.Fy;
                    hwe = min([1.92*obj.tw*sqrt(obj.Es/f)*(1-0.34*sqrt(obj.Es/f)/lambdaw) obj.hw]);
                    if hwe < 0
                        hwe = 0;
                    end
                    Qa = 1 - obj.tw*(obj.hw-hwe)/obj.A;
                else 
                    Qa = 1;
                end  
                Q = Qa*Qs;
                pnco = Q*obj.Fy*obj.A;
            end            
            
        end
        function pnc = Pnc(obj,axis)
            % Compressive Strength
            if strcmpi(axis,'min')
                pnc = min([obj.Pnc('strong') obj.Pnc('weak')]);
                return;
            elseif strcmpi(axis,'max')
                pnc = max([obj.Pnc('strong') obj.Pnc('weak')]);
                return;
            end
            % @todo - add torsional buckling (E4)
            Fe = pi^2*obj.Es/(obj.K(axis)*obj.L(axis)/obj.r(axis))^2;
            Pnco = obj.Pnco;
            pnc = Pnco*AISC_column_curve((Pnco/obj.A)/Fe);
        end
        function pnt = Pnt(obj)
            % Tensile Strength
            pnt = obj.A*obj.Fy;
        end
        function mno = Mno(obj,axis)
            % L=0 Moment Strength
            mno = obj.Mn(axis,1,0);
        end        
        function mn = Mn(obj,axis,Cb,Lb)
            % Flexural Strength 
            if ( nargin < 3 )
                if isempty(obj.Cb)
                    Cb = 1;
                else
                    Cb = obj.Cb;
                end
            end
            if ( nargin < 4 )
                if isempty(obj.Lb)
                    Lb = obj.L;
                else
                    Lb = obj.Lb;
                end
            end
            
            if obj.neglectLateralTorsionalBuckling && obj.neglectLocalBuckling
                mn = obj.Fy*obj.Z(axis);
                return
            end       
            
            lambdaf  = obj.lambda(axis,'flange');
            lambdaPf = 0.38*sqrt(obj.Es/obj.Fy);
            lambdaRf = 1.00*sqrt(obj.Es/obj.Fy);
            if ( lambdaf <= lambdaPf )
                flangeSlenderness = 'compact';
            elseif ( lambdaf <= lambdaRf )
                flangeSlenderness = 'noncompact';
            else 
                flangeSlenderness = 'slender';
            end
                  
            lambdaw  = obj.lambda(axis,'web');
            lambdaPw = 3.76*sqrt(obj.Es/obj.Fy);
            lambdaRw = 5.70*sqrt(obj.Es/obj.Fy);
            if ( lambdaw <= lambdaPw )
                webSlenderness = 'compact';
            elseif ( lambdaw <= lambdaRw )
                webSlenderness = 'noncompact';
            else 
                webSlenderness = 'slender';
            end            
           
            if obj.neglectLocalBuckling
                flangeSlenderness = 'compact';
                webSlenderness = 'compact';
            end
            if obj.neglectLateralTorsionalBuckling
                Lb = 0;
            end
            
            if strcmpi(axis,'strong')
                % Strong axis bending
                mp = obj.Fy*obj.Z(axis);
                
                % Lateral Torsional Buckling
                % Doubly symmetric I-shape: c = 1
                ho = obj.d-obj.tf;
                Cw = obj.I('weak')*ho^2/4;
                rts = sqrt(sqrt(obj.I('weak')*Cw)/obj.S('strong'));
                Lp = 1.76*obj.r('weak')*sqrt(obj.Es/obj.Fy);
                Lr = 1.95*rts*(obj.Es/0.7/obj.Fy)* ...
                    sqrt(obj.J/obj.S(axis)/ho)* ...
                    sqrt(1+sqrt(1+6.76*(0.7*obj.Fy/obj.Es * obj.S(axis)*ho/obj.J )));
                if ( Lb <= Lp )
                    mltb = mp;
                elseif ( Lb <= Lr )
                    mltb = Cb*( mp - (mp-0.7*obj.Fy*obj.S(axis))*((Lb-Lp)/(Lr-Lp)) );
                    mltb = min([mltb mp]);
                else
                    Fcr = Cb*pi^2*obj.Es / (Lb/rts)^2 * ...
                        sqrt(1+0.078*(obj.J/obj.S(axis)/ho)*(Lb/rts)^2);
                    mltb = Fcr*obj.S(axis);
                    mltb = min([mltb mp]);
                end     
                
                if strcmpi(webSlenderness,'compact')
                    if strcmpi(flangeSlenderness,'compact')                     
                        mn = min([mp mltb]);
                    else
                        if strcmpi(flangeSlenderness,'noncompact')
                            mcflb = mp - (mp-0.7*obj.Fy*obj.S(axis))*...
                                ((lambdaf-lambdaPf)/(lambdaRf-lambdaPf));
                        else
                            kc = 4/sqrt(lambdaw);
                            if ( kc < 0.35 ); kc = 0.35; end
                            if ( kc > 0.76 ); kc = 0.76; end
                            mcflb = 0.9*obj.Es*kc*obj.S(axis)/lambdaf^2;
                        end
                        mn = min([mp mltb mcflb]);  
                    end
                elseif strcmpi(webSlenderness,'noncompact')
                    warning('design:notImplemented',...
                        'Mn for wide flanges shapes with noncompact webs not yet implemented')
                    mn = 0;
                else
                    warning('design:notImplemented',...
                        'Mn for wide flanges shapes with slender webs not yet implemented')
                    mn = 0;
                end
            elseif strcmpi(axis,'weak') 
                % Weak axis bending
                mp = min([obj.Fy*obj.Z(axis) 1.6*obj.Fy*obj.S(axis)]);
                if strcmpi(flangeSlenderness,'compact')
                    mflb = mp;
                elseif strcmpi(flangeSlenderness,'noncompact')
                    mflb = mp - (mp - 0.7*obj.Fy*obj.S(axis))*...
                        ((lambdaf - lambdaPf)/(lambdaRf - lambdaPf));
                else 
                    Fcr = 0.69*obj.Es/lambdaf^2;
                    mflb = Fcr*obj.S(axis);
                end
                mn = min([mp mflb]);
            else
                error('Unknown axis');           
            end
            
        end
        function vn = Vn(obj,axis) 
            % Shear Strength
            switch lower(axis)
                case 'strong'
                    axis = 'strong';
                    Aw = obj.d*obj.tw;
                    lambda = obj.lambda(axis,'web');
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
                    lambda = obj.lambda(axis,'flange');
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
        function pnt = Pnt_expected(obj)
            pnt = obj.Ry*obj.Fy*obj.A;
        end
        function pnc = Pnc_expected(obj,axis)
            if strcmpi(axis,'min')
                pnc = min([obj.Pnc_expected('strong') obj.Pnc_expected('weak')]);
                return
            end     
            assert(obj.lambda(axis,'flange') < 0.56*sqrt(obj.Es/obj.Fy),...
                'Pnc_expected not yet implemented for sections with slender flanges');
            assert(obj.lambda(axis,'web') < 1.49*sqrt(obj.Es/obj.Fy),...
                'Pnc_expected not yet implemented for sections with slender webs'); 
            r = sqrt(obj.I(axis)/obj.A);
            Fe = pi^2*obj.Es/(obj.K(axis)*obj.L/r)^2;
            Fcre = obj.Ry*obj.Fy*AISC_column_curve(obj.Ry*obj.Fy/Fe);
            pnc = min([obj.Ry*obj.Fy*obj.A 1.14*Fcre*obj.A]);
        end
        function mn = Mn_expected(obj,axis)
            assert(obj.lambda(axis,'flange') < 0.38*sqrt(obj.Es/obj.Fy),...
                'Pnc_expected not yet implemented for sections with non-compact or slender flanges');
            assert(obj.lambda(axis,'web') < 3.78*sqrt(obj.Es/obj.Fy),...
                'Pnc_expected not yet implemented for sections with non-compact or slender webs');
            % @todo - check unbraced length?
            mn = obj.Ry*obj.Fy*obj.Z(axis);
        end
        
        %% Design Checks
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)          
            
            % Resistance factors
            phi_Pc = 0.90;
            phi_M  = 0.90;
            phi_Pt = 0.90;
            if ( ((obj.d-2*obj.tf)/obj.tw) <= 2.24*sqrt(obj.Es/obj.Fy) )
                phi_Vstrong  = 1.00;
            else
                phi_Vstrong  = 0.90;
            end
            phi_Vweak  = 0.90;
            
            % Moment Strengths
            if isempty(obj.Cb)
                if isempty(Ms)
                    Cb = 1;
                else
                    Cb = aisc2010.Cb(xi,Ms);
                end
            else
                Cb = obj.Cb;
            end
            if isempty(obj.Lb)
                Lb = obj.L; % If Lb is not defined, assume Lb = L
            else
                Lb = obj.Lb;
            end
            Mcs = phi_M*obj.Mn('strong',Lb,Cb); 
            Mcw = phi_M*obj.Mn('weak');
            
            % Compressive Load / Moment Interaction
            Pc = phi_Pc*obj.Pnc('min');
            ratio_PMc = aisc_H11_interaction_check(P,Ms,Mw,Pc,Mcs,Mcw);
            
            % Tensile Load / Moment Interaction
            Pc = phi_Pt*obj.Pnt;
            ratio_PMt = aisc_H12_interaction_check(P,Ms,Mw,Pc,Mcs,Mcw);
                 
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
                case 'highlyductilesection'
                    Pu = varargin{1};
                    
                    % Web
                    lambdaw = obj.lambda('strong','web');
                    if strcmpi(Pu,'brace')
                        lambdaHDw = 1.49*sqrt(obj.Es/obj.Fy);
                    else
                        phi_Pc = 0.90;
                        Ca = Pu/(phi_Pc*obj.Pnco);
                        if Ca <= 0.125
                            lambdaHDw = 2.45*sqrt(obj.Es/obj.Fy)*(1-0.93*Ca);
                        else
                            lambdaHDw = max([...
                                0.77*sqrt(obj.Es/obj.Fy)*(2.93-Ca)
                                1.49*sqrt(obj.Es/obj.Fy)]);
                        end
                    end
                    tf_web = lambdaw <= lambdaHDw;
                    % Flange
                    lambdaf = obj.lambda('strong','flange');
                    lambdaHDf = 0.30*sqrt(obj.Es/obj.Fy);
                    tf_flange = lambdaf <= lambdaHDf;
                    
                    pass_tf = tf_web && tf_flange;
                case 'moderatelyductilesection'
                    Pu = varargin{1};
                    
                    % Web
                    lambdaw = obj.lambda('strong','web');
                    if strcmpi(Pu,'brace')
                        lambdaHDw = 1.49*sqrt(obj.Es/obj.Fy);
                    else
                        phi_Pc = 0.90;
                        Ca = Pu/(phi_Pc*obj.Pnco);
                        if Ca <= 0.125
                            lambdaHDw = 3.76*sqrt(obj.Es/obj.Fy)*(1-2.75*Ca);
                        else
                            lambdaHDw = max([...
                                1.12*sqrt(obj.Es/obj.Fy)*(2.33-Ca)
                                1.49*sqrt(obj.Es/obj.Fy)]);
                        end
                    end
                    tf_web = lambdaw <= lambdaHDw;
                    % Flange
                    lambdaf = obj.lambda('strong','flange');
                    lambdaHDf = 0.38*sqrt(obj.Es/obj.Fy);
                    tf_flange = lambdaf <= lambdaHDf;
                    
                    pass_tf = tf_web && tf_flange;  
                case 'moderatelyductilebracespacing'
                    pass_tf = obj.L('weak') < 0.17*obj.r('weak')*obj.Es/obj.Fy;
                case 'highlyductilebracespacing'
                    pass_tf = obj.L('weak') < 0.086*obj.r('weak')*obj.Es/obj.Fy;
                case 'slendernessratio_klr'
                    limit = varargin{1};                    
                    KLr_strong = obj.K('strong')*obj.L('strong')/sqrt(obj.I('strong')/obj.A);
                    KLr_weak = obj.K('weak')*obj.L('weak')/sqrt(obj.I('weak')/obj.A);
                    pass_tf = max([KLr_strong KLr_weak]) < limit;
                case 'requiredmomentofinertia'
                    axis        = varargin{1};
                    lowerLimit  = varargin{2};
                    pass_tf = lowerLimit < obj.I(axis);
                otherwise
                    error('Unknown checkType');
            end
        end
        
        %% Interaction Strength
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case {'aisc','h1.1'}
                    [P,M] = AISC_H1_interaction_diagram(...
                        obj.Pnt,-obj.Pnco,obj.Mno(axis),quadrant);
                otherwise
                    error('Unknown type: %s',type);
            end            
        end
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case {'aisc','h1.1'}
                    [P,M] = AISC_H1_interaction_diagram(...
                        obj.Pnt,-obj.Pnc(axis),obj.Mn(axis,1,0),quadrant);
                case 'factoredaisc'
                    [P,M] = obj.beamColumnInteraction2d(axis,'AISC',quadrant);
                    P = 0.9*P;
                    M = 0.9*M;
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        
        %% Export and Information Functions
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(...
                obj,axis,type)
            switch lower(type)
                case {'gross','columnstrength'}
                    E = obj.Es;
                    A = obj.A;
                    I = obj.I(axis);
                otherwise
                    error('Unknown type');
            end
        end
        function [E,A,Iz,Iy,GJ] = sectionPropertiesForElasticAnalysis3d(...
                obj,type)
            switch lower(type)
                case {'gross','columnstrength'}
                    E = obj.Es;
                    A = obj.A;
                    Iz = obj.I('strong');
                    Iy = obj.I('weak');
                    GJ = obj.G*obj.J;
                otherwise
                    error('Unknown type');
            end            
        end
        function lp = Lp(obj,axis,Li)
            switch lower(axis)
                case 'strong'
                    assert(~isempty(obj.Fu),'Fu is required to compute Lp');
                    lp = (0.405-0.0033*(obj.hw/obj.tw)-0.0268*(obj.bf/2/obj.tf)+0.184*(obj.Fu/obj.Fy-1))*Li;
                case 'weak'
                    error('Lp not implemented for weak axis bending');
                otherwise
                    error('Bad axis'); 
            end             
        end        
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');
            switch lower(axis)
                case 'strong'
                    yExtreme = obj.d/2;
                case 'weak'
                    yExtreme = obj.bf/2;
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
                otherwise
                    error('Unknown type');
            end
        end
        function strain = longitudinalStrain3d(obj,axialStrain,curvatureY,curvatureZ,type)
            assert(isequal(size(axialStrain),size(curvatureY),size(curvatureZ)),...
                'axialStrain and curvature should be the same size');

            z = obj.bf/2;
            y = obj.d/2;
            strain_s1 = axialStrain + z*curvatureY + y*curvatureZ;
            strain_s2 = axialStrain + z*curvatureY - y*curvatureZ;
            strain_s3 = axialStrain - z*curvatureY - y*curvatureZ;
            strain_s4 = axialStrain - z*curvatureY + y*curvatureZ;

            switch lower(type)
                case 'maxcompressive'
                    strain = min([strain_s1(:) strain_s2(:) strain_s3(:) strain_s4(:)],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxtensile'
                    strain = max([strain_s1(:) strain_s2(:) strain_s3(:) strain_s4(:)],[],2);
                    strain = reshape(strain,size(axialStrain));
                case 'maxabsolute'
                    strain = max(abs([strain_s1(:) strain_s2(:) strain_s3(:) strain_s4(:)]),[],2);
                    strain = reshape(strain,size(axialStrain));
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function x = getSectionData(obj,type,axis)
            switch lower(type)
                case 'steelstrength'
                    x = obj.Fy;
                case 'grosssteelflexuralrigidity'
                    x = obj.Es*obj.I(axis);
                case 'grosssectioncompressionstrength'
                    x = obj.A*obj.Fy;
                case 'shapefactor'
                    x = obj.Z(axis)/obj.S(axis);
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
                        description = obj.shapeName;
                    else
                        description = sprintf('WF d = %g, Fy = %g',...
                            obj.d,obj.Fy);
                    end
                otherwise
                    error('Unknown flag');
            end
        end  
        function fs = fiberSectionObject(obj,id)
            fs = fiberSection;
            fs.addWFShape(id,obj.d,obj.tw,obj.bf,obj.tf,obj.k);
        end
        function psd = plasticStressDistributionObject(obj)
            idSteel = 1;
            fs = obj.fiberSectionObject(id);
            psd = plastic_stress_distribution(fs);       
            psd.addMaterial(idSteel,obj.Fy,-obj.Fy);
        end        
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            angles = linspace(0,pi/2,25);
            r = obj.k-obj.tf;
            xo = obj.tw/2 + r;
            yo = obj.d/2 - obj.k;
            x = [-0.5*obj.bf -0.5*obj.bf 0.5*obj.bf 0.5*obj.bf ...
                xo+r*cos(angles+pi/2) xo+r*cos(angles+pi) ...
                0.5*obj.bf 0.5*obj.bf -0.5*obj.bf -0.5*obj.bf ...
                -xo+r*cos(angles+3*pi/2) -xo+r*cos(angles) ...
                -0.5*obj.bf];
            y = [0.5*obj.d-obj.tf 0.5*obj.d 0.5*obj.d 0.5*obj.d-obj.tf ...
                yo+r*sin(angles+pi/2) -yo+r*sin(angles+pi) ...
                -0.5*obj.d+obj.tf -0.5*obj.d -0.5*obj.d -0.5*obj.d+obj.tf ...
                -yo+r*sin(angles+3*pi/2) yo+r*sin(angles) ...
                0.5*obj.d-obj.tf ];
            fill(x,y,obj.color_steelFill,'LineStyle','none')
            plot(x,y,'k-','LineWidth',lineWidth);
            axis equal
        end        
    end
    
    methods (Static)
        function type = memberType()
            type = 'wf';
        end 
        function tf = hasConcrete()
            tf = false;
        end
    end    
    
end
