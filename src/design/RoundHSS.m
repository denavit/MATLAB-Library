classdef RoundHSS < structural_shape
    % RoundHSS
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
        D               % Outside diameter of the steel tube
        t               % Thickness of the steel tube
        % Material Properties
        Fy              % Yield stress of the steel tube
        Fu = [];        % Ultimate stress of the steel tube
        Es              % Modulus of Elasticity of Steel 
        G               % Shear Modulus of Elasticity of Steel
        % Expected Strength Parameters
        Ry = [];
        Rt = [];
        % Information
        shapeName = ''; % Name of the steel shape
        % Design Options
        neglectLocalBuckling = false;
    end
    
    methods
        %% Constructor
        function obj=RoundHSS(D,t,Fy,units)
            obj.D  = D;
            obj.t  = t;
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
                    warning('Design:Input','Unknown unit system, setting Es = 0')
                    obj.Es = 0;
                    obj.G  = 0;
            end
        end
        
        %% Geometric Properties
        function A = A(obj)
            A = (pi/4)*(obj.D^2-obj.Di^2);
        end
        function I = I(obj,axis)
            I = (pi/64)*(obj.D^4-obj.Di^4);
        end   
        function Z = Z(obj,axis)
            Z = (1/6)*(obj.D^3-obj.Di^3);
        end
        function S = S(obj,axis)
            S = (pi/32)*(obj.D^4-obj.Di^4)/obj.D;
        end
        function j = J(obj)
            j = (pi/32)*(obj.D^4-obj.Di^4);
        end
        function d = depth(obj,axis)
            d = obj.D;
        end
        function d = Di(obj)
            d = obj.D-2*obj.t;
        end
        
        %% Design Strengths  
        function pnco = Pnco(obj)
            % Stub Column (L=0) Strength
            if obj.neglectLocalBuckling
                pnco = obj.Fy*obj.A;
                return
            end
            
            tubeSlenderness = obj.D/obj.t;

            if tubeSlenderness < 0.11*obj.Es/obj.Fy
                % Without slender elements
                pnco = obj.Fy*obj.A;
            elseif tubeSlenderness < 0.45*obj.Es/obj.Fy
                % With slender elements
                Q = (0.038*obj.Es)/(obj.Fy*(obj.D/obj.t)) + 2/3;
                pnco = Q*obj.Fy*obj.A;  
            else
                pnco = 0;
            end            
        end
        function pn = Pnc(obj,axis)
            % Compressive Strength 
            if strcmpi(axis,'min')
                pn = min([obj.Pnc('strong') obj.Pnc('weak')]);
                return;
            elseif strcmpi(axis,'max')
                pn = max([obj.Pnc('strong') obj.Pnc('weak')]);
                return;
            end
            r = sqrt(obj.I(axis)/obj.A);
            Fe = pi^2*obj.Es/(obj.K(axis)*obj.L/r)^2;
            Pnco = obj.Pnco;
            pn = Pnco*AISC_column_curve((Pnco/obj.A)/Fe);          
        end
        function pnt = Pnt(obj)
            % Tensile Strength                        
            pnt = obj.A*obj.Fy;
        end
        function mno = Mno(obj,axis)
            % L=0 Moment Strength
            mno = obj.Mn(axis);
        end        
        function mn = Mn(obj,axis)
            % Flexural Strength        
            if obj.neglectLocalBuckling
                mn = obj.Z(axis)*obj.Fy;
                return
            end           
            
            mp = obj.Z(axis)*obj.Fy;
            
            tubeSlenderness = obj.D/obj.t;

            if tubeSlenderness < 0.07*obj.Es/obj.Fy
                % Compact Tube
                mn = mp;
            elseif tubeSlenderness < 0.31*obj.Es/obj.Fy
                % Noncompact Tube
                mnLB = ((0.021*obj.Es)/(obj.D/obj.t) + obj.Fy)*obj.S(axis);
                mn = min([mp mnLB]);
            else
                % Slender Tube
                Fcr = (0.33*obj.Es)/(obj.D/obj.t);
                mnLB = Fcr*obj.S(axis); 
                mn = min([mp mnLB]);
            end
        end
        function vn = Vn(obj,axis,Lv) 
            % Shear Strength           
            if nargin < 3
                Lv = Inf;
            end 
            Fcr1 = (1.60*obj.Es)/(sqrt(Lv/obj.D)*(obj.D/obj.t)^5/4);
            Fcr2 = (0.78*obj.Es)/(obj.D/obj.t)^3/2;
            Fcr = min([max([Fcr1 Fcr2]) 0.6*obj.Fy]);
            vn = Fcr*obj.A/2;
        end
        function pnt = Pnt_expected(obj)
            pnt = obj.Ry*obj.Fy*obj.A;
        end
        function pnc = Pnc_expected(obj,axis)
            if strcmpi(axis,'min')
                pnc = min([obj.Pnc_expected('strong') obj.Pnc_expected('weak')]);
                return
            end
            tubeSlenderness = obj.D/obj.t;
            assert(tubeSlenderness < 0.11*obj.Es/obj.Fy,...
                'Pnc_expected not yet implemented for tubes with slender elements');
            r = sqrt(obj.I(axis)/obj.A);
            Fe = pi^2*obj.Es/(obj.K(axis)*obj.L/r)^2;
            Fcre = obj.Ry*obj.Fy*AISC_column_curve(obj.Ry*obj.Fy/Fe);
            pnc = min([obj.Ry*obj.Fy*obj.A 1.14*Fcre*obj.A]);
        end
        
        %% Design Checks
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)
            
            % Resistance factors
            phi_Pc = 0.90;
            phi_M  = 0.90;
            phi_Pt = 0.90;
            phi_V  = 0.90;
            
            % Moment Strengths
            Mcs = phi_M*obj.Mn('strong');
            Mcw = phi_M*obj.Mn('weak');
            
            % Compressive Load / Moment Interaction
            Pc = phi_Pc*obj.Pnc('min');
            ratio_PMc = aisc2010.interactionCheck_H11(P,Ms,Mw,Pc,Mcs,Mcw);
            
            % Tensile Load / Moment Interaction
            Pc = phi_Pt*obj.Pnt;
            ratio_PMt = aisc2010.interactionCheck_H12(P,Ms,Mw,Pc,Mcs,Mcw);
            
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
                error('Torsion strength check not implemented');
            end
            
            ratio = max([ratio_PMc ratio_PMt ratio_Vs ratio_Vw]);
        end
        function pass_tf = proportioningCheck(obj,checkType,varargin)
            switch lower(checkType)
                case 'highlyductilesection'                   
                    pass_tf = obj.D/obj.t < 0.038*obj.Es/obj.Fy;
                case 'moderatelyductilesection'                   
                    pass_tf = obj.D/obj.t < 0.044*obj.Es/obj.Fy;
                case 'slendernessratio_klr'
                    limit = varargin{1};
                    KLr_strong = obj.Kstrong*obj.L/sqrt(obj.I('strong')/obj.A);
                    KLr_weak = obj.Kweak*obj.L/sqrt(obj.I('weak')/obj.A);
                    pass_tf = max([KLr_strong KLr_weak]) < limit;
                otherwise
                    error('Unknown checkType');
            end
        end
        
        %% Interaction Strength
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            switch lower(type)
                case {'aisc','h1.1'}
                    [P,M] = AISC_H1_interaction_diagram(...
                        obj.Pnt,-obj.Pnco,obj.Mno(axis),quadrant);
                otherwise
                    error('Unknown type: %s',type);
            end            
        end  
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            switch lower(type)
                case {'aisc','h1.1'}
                    [P,M] = AISC_H1_interaction_diagram(...
                        obj.Pnt,-obj.Pnc(axis),obj.Mn(axis),quadrant);
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
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');
            switch lower(axis)
                case 'strong'
                    yExtreme = obj.D/2;
                case 'weak'
                    yExtreme = obj.D/2;
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
                        description = sprintf('HSS%.3fx%.3f',obj.D,obj.t);
                    end
                    
                otherwise
                    error('Unknown flag');
            end
        end 
        function fs = fiberSectionObject(obj,id)
            fs = fiberSection;
            fs.addPatch('circ',id,0,0,obj.Di/2,obj.D/2);
        end
        function psd = plasticStressDistributionObject(obj)
            idSteel = 1;
            fs = obj.fiberSectionObject(idSteel);
            psd = plastic_stress_distribution(fs);       
            psd.addMaterial(idSteel,obj.Fy,-obj.Fy);
        end
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            angles = linspace(0,2*pi,100);
            xo = (obj.D/2)*cos(angles);
            yo = (obj.D/2)*sin(angles);
            xi = (obj.Di/2)*cos(angles);
            yi = (obj.Di/2)*sin(angles);
            fill([xi fliplr(xo) xi(1)],[yi fliplr(yo) yi(1)],...
                obj.color_steelFill,'LineStyle','none')
            plot(xo,yo,'k-','LineWidth',lineWidth);
            plot(xi,yi,'k-','LineWidth',lineWidth);
            axis equal            
        end        
    end
    
    methods (Static)
        function type = memberType()
            type = 'roundhss';
        end  
        function tf = hasConcrete()
            tf = false;
        end 
    end
    
end
