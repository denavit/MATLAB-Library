classdef RectangularHSS < structural_shape
    % RectangularHSS
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
        H               % Depth of the steel tube
        B               % Width of the steel tube
        t               % Thickness of the steel tube
        ri              % Internal radius of the corner of the steel tube
                        %   ri = 0 indicates a box section        
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
        function obj=RectangularHSS(H,B,t,Fy,units)
            obj.H  = H;
            obj.B  = B;
            obj.t  = t;
            obj.Fy = Fy;
            obj.units = units;
            
            % Default value of internal radius
            obj.ri = obj.t;            
            
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
        function ro = ro(obj)
            if obj.ri == 0
                ro = 0;
            else
                ro = obj.ri+obj.t;
            end
        end
        function area = A(obj)
            shp = Rectangular_Tube_Shape(obj.H,obj.B,obj.t,obj.ro);
            area = shp.A;
        end
        function i = I(obj,axis)
            shp = Rectangular_Tube_Shape(obj.H,obj.B,obj.t,obj.ro);
            i = shp.I(axis);     
        end  
        function s = S(obj,axis)
            shp = Rectangular_Tube_Shape(obj.H,obj.B,obj.t,obj.ro);
            s = shp.S(axis);
        end         
        function z = Z(obj,axis)
            shp = Rectangular_Tube_Shape(obj.H,obj.B,obj.t,obj.ro);
            z = shp.Z(axis); 
        end
        function j = J(obj)
            shp = Rectangular_Tube_Shape(obj.H,obj.B,obj.t,obj.ro);
            j = shp.J(axis); 
        end
        function [d,b] = depth(obj,axis)
            switch lower(axis)
                case {'z','x','strong'}
                    d = obj.H;
                    b = obj.B;
                case {'y','weak'}
                    d = obj.B;
                    b = obj.H;
                otherwise
                    error('Bad axis');
            end
            if nargout < 2
                clear b;
            end
        end
        
        %% Design Strengths  
        function pnco = Pnco(obj)
            % Stub Column (L=0) Strength
            if obj.neglectLocalBuckling
                pnco = obj.Fy*obj.A;
                return
            end
            
            flangeSlenderness = obj.B/obj.t-3;
            webSlenderness = obj.H/obj.t-3;

            if ( max([flangeSlenderness webSlenderness]) < 1.40*sqrt(obj.Es/obj.Fy) )
                % Without slender elements
                pnco = obj.Fy*obj.A;
            else
                % With slender elements
                f = obj.Fy;  % conservative assumption given in user note
                if ( webSlenderness > 1.40*sqrt(obj.Es/obj.Fy) ) 
                    he = 1.92*obj.t*sqrt(obj.Es/f)*(1-0.34/webSlenderness*sqrt(obj.Es/f));
                    he = min([obj.H he]);
                else
                    he = obj.H;
                end
                if ( flangeSlenderness > 1.40*sqrt(obj.Es/obj.Fy) ) 
                    be = 1.92*obj.t*sqrt(obj.Es/f)*(1-0.34/flangeSlenderness*sqrt(obj.Es/f));
                    be = min([obj.B be]);
                else
                    be = obj.B;
                end
                Aeff = obj.A - 2*(obj.B-be)*obj.t - 2*(obj.H-he)*obj.t;
                Q = Aeff/obj.A;
                pnco = Q*obj.Fy*obj.A;                
            end            
        end
        function pn = Pnc(obj,axis)
            % Compressive Strength 
            if strcmpi(axis,'min')
                pn = min([obj.Pnc('strong') obj.Pnc('weak')]);
                return
            end
            r = sqrt(obj.I(axis)/obj.A);
            Fe = pi^2*obj.Es/(obj.K(axis)*obj.L(axis)/r)^2;
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
            switch lower(axis)
                case {'z','x','strong'}
                    b = obj.B - 2*obj.t;
                    d = obj.H - 2*obj.t;
                    Z = obj.Z('strong');
                    S = obj.I('strong')/(d/2);
                case {'y','weak'}
                    b = obj.H - 2*obj.t;
                    d = obj.B - 2*obj.t; 
                    Z = obj.Z('weak');
                    S = obj.I('weak')/(d/2);
                otherwise
                    error('Bad axis'); 
            end 
            
            if obj.neglectLocalBuckling
                mn = Z*obj.Fy;
                return
            end
            
            mp = Z*obj.Fy;
            
            flangeSlenderness = b/obj.t-3;
            if ( flangeSlenderness < 1.12*sqrt(obj.Es/obj.Fy) )
                % Compact Flange
                mnFLB = mp;
            elseif ( flangeSlenderness < 1.40*sqrt(obj.Es/obj.Fy) )
                % Noncompact Flange
                mnFLB = mp - (mp-obj.Fy*S)*(3.57*(b/obj.t)*sqrt(obj.Es/obj.Fy)-4.0);
                mnFLB = min([mp mnFLB]);
            else
                % Slender Flange
                be = 1.92*obj.t*sqrt(obj.Es/obj.Fy)*(1-0.38/(b/obj.t)*sqrt(obj.Es/obj.Fy));
                be = min([b be]);
                Ieff = sectionProperties('Iy','RoundedRectangularTube',be,d,obj.t,2*obj.t);
                Seff = Ieff/(d/2);
                mnFLB = obj.Fy*Seff;
            end
            
            webSlenderness = d/obj.t-3;
            if ( webSlenderness < 2.42*sqrt(obj.Es/obj.Fy) )
                % Compact Web
                mnWLB = mp;
            elseif ( webSlenderness < 5.70*sqrt(obj.Es/obj.Fy) )
                % Noncompact Web
                h = d-3*obj.t;
                mnWLB = mp - (mp-obj.Fy*S)*(0.305*(h/obj.t)*sqrt(obj.Es/obj.Fy)-0.738);
                mnFLB = min([mp mnFLB]);
            else
                % Slender Web
                mnWLB = 0;
            end            
                
            % Minimum of the three limit states
            mn = min([mp mnFLB mnWLB]);
        end
        function vn = Vn(obj,axis) 
            % Shear Strength
            switch lower(axis)
                case {'z','x','strong'}
                    d = obj.H - 2*obj.t;
                case {'y','weak'}
                    d = obj.B - 2*obj.t; 
                otherwise
                    error('Bad axis'); 
            end             
          
            h = d-4*obj.t;
            Aw = 2*h*obj.t;
            kv = 5;
            if ( h/obj.t <= 1.10*sqrt(kv*obj.Es/obj.Fy) )
                Cv = 1.0;
            elseif ( h/obj.t <= 1.37*sqrt(kv*obj.Es/obj.Fy) )
                Cv = 1.10*sqrt(kv*obj.Es/obj.Fy)/(h/obj.t);
            else
                Cv = 1.51*obj.Es*kv/((h/obj.t)^2*obj.Fy);
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
            flangeSlenderness = obj.B/obj.t-3;
            webSlenderness = obj.H/obj.t-3;
            assert(max([flangeSlenderness webSlenderness]) < 1.40*sqrt(obj.Es/obj.Fy),...
                'Pnc_expected not yet implemented for tubes with slender elements');
            r = sqrt(obj.I(axis)/obj.A);
            Fe = pi^2*obj.Es/(obj.K(axis)*obj.L(axis)/r)^2;
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
            ratio_PMc = aisc_H11_interaction_check(P,Ms,Mw,Pc,Mcs,Mcw);
            
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
                case 'highlyductilesection'                   
                    pass_tf = max([obj.B/obj.t-3 obj.H/obj.t-3]) < 0.55*sqrt(obj.Es/obj.Fy);
                case 'moderatelyductilesection'                   
                    pass_tf = max([obj.B/obj.t-3 obj.H/obj.t-3]) < 0.64*sqrt(obj.Es/obj.Fy);
                case 'slendernessratio_klr'
                    limit = varargin{1};
                    KLr_strong = obj.K('strong')*obj.L('strong')/sqrt(obj.I('strong')/obj.A);
                    KLr_weak = obj.K('weak')*obj.L('weak')/sqrt(obj.I('weak')/obj.A);
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
        function [E,A,Iz,Iy,G,J] = sectionPropertiesForElasticAnalysis3d(obj,type)
            switch lower(type)
                case {'gross','columnstrength'}
                    E  = obj.Es;
                    A  = obj.A;
                    Iz = obj.I('strong');
                    Iy = obj.I('weak');
                    G  = obj.G;
                    J  = obj.J;
                otherwise
                    error('Unknown type');
            end
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');
            switch lower(axis)
                case {'z','x','strong'}
                    yExtreme = obj.H/2;
                case {'y','weak'}
                    yExtreme = obj.B/2;
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
                        description = sprintf('HSS%gx%gx%g',...
                            obj.H,obj.B,obj.t);
                    end
                    
                otherwise
                    error('Unknown flag');
            end
        end 
        function fs = fiberSectionObject(obj,id)
            fs = fiberSection;
            fs.addRectHSS(id,obj.H,obj.B,obj.t);
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
            angles = linspace(0,pi/2,25);
            d1 = obj.H/2-2*obj.t;
            b1 = obj.B/2-2*obj.t;
            xo = [ ( b1+2*obj.t*cos(angles)) ...
                  (-b1+2*obj.t*cos(angles+pi/2)) ...
                  (-b1+2*obj.t*cos(angles+pi)) ...
                  ( b1+2*obj.t*cos(angles+1.5*pi)) b1+2*obj.t ];
            yo = [ ( d1+2*obj.t*sin(angles)) ...
                  ( d1+2*obj.t*sin(angles+pi/2)) ...
                  (-d1+2*obj.t*sin(angles+pi)) ...
                  (-d1+2*obj.t*sin(angles+1.5*pi)) d1 ]; 
            xi = [ ( b1+obj.t*cos(angles)) ...
                  (-b1+obj.t*cos(angles+pi/2)) ...
                  (-b1+obj.t*cos(angles+pi)) ...
                  ( b1+obj.t*cos(angles+1.5*pi)) b1+obj.t];
            yi = [ ( d1+obj.t*sin(angles)) ...
                  ( d1+obj.t*sin(angles+pi/2)) ...
                  (-d1+obj.t*sin(angles+pi)) ...
                  (-d1+obj.t*sin(angles+1.5*pi)) d1 ];
            fill([xi fliplr(xo) xi(1)],[yi fliplr(yo) yi(1)],...
                obj.color_steelFill,'LineStyle','none')
            plot(xi,yi,'k-','LineWidth',lineWidth);            
            plot(xo,yo,'k-','LineWidth',lineWidth);            
            axis equal
        end        
    end
    
    methods (Static)
        function type = memberType()
            type = 'recthss';
        end  
        function tf = hasConcrete()
            tf = false;
        end 
    end
    
end
