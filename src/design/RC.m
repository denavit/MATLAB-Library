classdef RC < structural_shape
    
    properties
        fc
        Ec
        
        fy
        Es
        
        conc_cross_section
        reinforcement
        
        transverse_reinf_type = '';

        treat_reinforcement_as_point = true;
    end
    
    methods
        function obj = RC(fc,fy,conc_cross_section,reinforcement,units)
            obj.fc = fc;
            obj.fy = fy;
            obj.conc_cross_section = conc_cross_section;
            if iscell(reinforcement)
                obj.reinforcement = reinforcement;
            else
                obj.reinforcement = {reinforcement};
            end
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
        function d = depth(obj,axis)
            d = obj.conc_cross_section.depth(axis);
        end
        function a = Ag(obj)
            a = obj.conc_cross_section.A;
        end
        function a = Ac(obj)
            a = obj.Ag-obj.Asr;
        end
        function a = Asr(obj)
            a = 0;
            for i = 1:length(obj.reinforcement)
                a = a + obj.reinforcement{i}.num_bars*obj.reinforcement{i}.Ab;
            end
        end
        function i = Ig(obj,axis)
            i = obj.conc_cross_section.I(axis);
        end
        function i = Ic(obj,axis)
            i = obj.Ig(axis)-obj.Isr(axis);
        end
        function i = Isr(obj,axis)
            i = 0;
            for j = 1:length(obj.reinforcement)
                i = i + obj.reinforcement{j}.I(axis);
            end
        end
        function r = rg(obj,axis)
            r = obj.conc_cross_section.r(axis);
        end
        function po = Po(obj)
            po = 0.85*obj.fc*(obj.Ag-obj.Asr)+obj.fy*obj.Asr;
        end
        function pnco = Pnco(obj)
            % See Section 22.4.2 of ACI318-19
            switch lower(obj.transverse_reinf_type)
                case 'ties'
                    pnco = 0.80*obj.Po;
                case {'spiral','spirals'}
                    pnco = 0.85*obj.Po;
                otherwise
                    error('Unknown transverse_reinf_type: %s',obj.transverse_reinf_type);
            end
        end
        function pnc = Pnc(obj,axis)
            pnc = Pnco(obj);
        end
        function pnt = Pnt(obj)
            error('Not yet implemented')
        end
        function mno = Mno(obj,axis)
            error('Not yet implemented')
        end
        function f = phi(obj,et)
            f = ACI_phi(obj.transverse_reinf_type,et,obj.fy/obj.Es);
        end
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case 'aci'
                    scACI = obj.strainCompatibilityAciObject;
                    switch lower(axis)
                        case {'x','z'}
                            [P,M,~] = scACI.interactionSweep(0,50);
                        case 'y'
                            [P,~,M] = scACI.interactionSweep(pi/2,50);
                        otherwise
                            error('Unknown axis: %s',axis);
                    end
                case 'factoredaci'
                    scACI = obj.strainCompatibilityAciObject;
                    switch lower(axis)
                        case {'x','z'}
                            [P,M,~,et] = scACI.interactionSweep(0,50);
                        case 'y'
                            [P,~,M,et] = scACI.interactionSweep(pi/2,50);
                        otherwise
                            error('Unknown axis: %s',axis);
                    end
                    phi = obj.phi(et);
                    P = phi.*P;
                    M = phi.*M;
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case {'aci','factoredaci'}
                    [P,M] = sectionInteraction2d(obj,axis,type,quadrant);
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(obj,axis,type)
            switch lower(type)
                case 'gross'
                    E = obj.Ec;
                    I = (obj.Es*obj.Isr(axis) + obj.Ec*obj.Ic(axis))/E;
                    A = (obj.Es*obj.Asr + obj.Ec*obj.Ac)/E;
                case '0.7ecig'
                    E = obj.Ec;
                    I = 0.7*obj.Ig(axis);
                    A = obj.Ag;
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [E,A,Iz,Iy,G,J] = sectionPropertiesForElasticAnalysis3d(obj,type)
            error('Not yet implemented')
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            assert(isequal(size(axialStrain),size(curvature)),...
                'axialStrain and curvature should be the same size');
            yExtreme = obj.depth(axis)/2;
            switch lower(type)
                case {'maxcompressive','maxconcretecompressive'}
                    strain = min(axialStrain+yExtreme*curvature,...
                        axialStrain-yExtreme*curvature);
                case {'maxtensile','maxconcretetensile'}
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
                case 'reinforcingratio'
                    x = obj.Asr/obj.Ag;
                case 'concretestrength'
                    x = obj.fc;
                otherwise
                    x = NaN;
            end
        end
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)
            error('Not yet implemented')
        end
        function pass_tf = proportioningCheck(obj,checkType,varagrin)
            error('Not yet implemented')
        end
        function description = sectionDescription(obj,flag)
            error('Not yet implemented')
        end
        function fs = fiberSectionObject(obj,idConc,idReinf)
            fs = fiberSection;
            obj.conc_cross_section.add_to_fiber_section(fs,idConc);
            for i = 1:length(obj.reinforcement)
                Ab    = obj.reinforcement{i}.Ab;
                [z,y] = obj.reinforcement{i}.coordinates;
                for j = 1:length(z)
                    if obj.treat_reinforcement_as_point
                        fs.addFiber( idConc,-Ab,z(j),y(j))
                        fs.addFiber(idReinf, Ab,z(j),y(j))
                    else
                        error('Not yet implemented');
                    end
                end
            end
        end        
        function psd = plasticStressDistributionObject(obj)
            error('Not yet implemented')
        end
        function scACI = strainCompatibilityAciObject(obj)
            idConc  = 1;
            idReinf = 2;
            fs = obj.fiberSectionObject(idConc,idReinf);
            scACI = ACI_strain_compatibility(fs);
            
            % Add Concrete Boundaries
            [x,y,r] = obj.conc_cross_section.boundary_points();
            for i = 1:length(x)
                scACI.addConcreteBoundary(x(i),y(i),r(i));
            end
            
            % Add Steel Boundaries
            for i = 1:length(obj.reinforcement)
                [z,y] = obj.reinforcement{i}.coordinates;
                for j = 1:length(z)
                    if obj.treat_reinforcement_as_point
                        scACI.addSteelBoundary(z(j),y(j),0); 
                    else 
                        error('Not yet implemented');
                    end
                end
            end
            
            scACI.maxCompressiveStrength = -0.85*(obj.Asr*obj.fy + 0.85*obj.Ac*obj.fc);
            
            scACI.addMaterial('concrete', idConc,obj.fc,obj.units);
            scACI.addMaterial(   'steel',idReinf,obj.fy,obj.Es);
        end           
        
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            obj.conc_cross_section.plotSection(lineWidth);
            for i = 1:length(obj.reinforcement)
                obj.reinforcement{i}.plotSection(lineWidth);
            end
            axis equal
        end
    end    
    methods (Static)
        function type = memberType()
            type = 'RC';
        end
        function tf = hasConcrete()
            tf = true;
        end
    end    
    
    
end
