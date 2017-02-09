classdef UniaxialMaterialAnalysis < OpenSeesAnalysis
    
    properties
        uniaxial_material_def
        uniaxial_material_tag = 1;
    end
    
    methods
        %% Constructor
        function obj = UniaxialMaterialAnalysis(def,tag)
            obj.uniaxial_material_def = def;
            if nargin > 1
                obj.uniaxial_material_tag = tag;
            end
        end
        
        %% Set and get functions
        function set.uniaxial_material_def(obj,def)
            if ischar(def)
                if exist(def,'file') == 2
                    fid = fopen(def);
                    A = textscan(fid,'%s\n','Whitespace','\n');
                    obj.uniaxial_material_def = A{1};
                else
                    obj.uniaxial_material_def = {def};
                end
            elseif iscell(def) && isvector(def)
                obj.uniaxial_material_def = def;
            else
                error('uniaxial_material_def should be a cell vector');
            end
        end
        
        %% Run Analysis
        function results = runAnalysis(obj,peakPoints,rateType,rateValue)

            if nargin < 3
                rateType  = 'None';
            end
            
            if size(peakPoints,2) == 1
                peakPoints = peakPoints';
            end
            
            switch rateType
                case 'StrainRate'
                    rate = rateValue;
                case 'Steps'
                    rate = sum(abs(peakPoints-[0 peakPoints(1:end-1)]))/rateValue;
                case 'None'
                    rate = max(abs(peakPoints-[0 peakPoints(1:end-1)]));
                otherwise
                    error('Unknown rateType');
            end          
            
            % Filenames
            filename_input          = obj.scratchFile('runMaterialAnalysisInputFile.tcl');
            filename_disp           = obj.scratchFile('runMaterialAnalysisPattern.txt');
            filename_output_force   = obj.scratchFile('runMaterialAnalysisPatternForce.out');
            filename_output_def     = obj.scratchFile('runMaterialAnalysisPatternDisp.out');
            filename_output_stiff   = obj.scratchFile('runMaterialAnalysisPatternStiff.out');
            
            % Generate imposed displacment history
            numbers = fillOutNumbers(peakPoints,rate);
            numbers = vertcat(numbers(1),numbers,numbers(end));
            numSteps = length(numbers)-2;
            save(filename_disp,'numbers','-ASCII');
            
            % Create .tcl file
            fid = fopen(filename_input,'w');
            fprintf(fid,'model BasicBuilder -ndm 1 -ndf 1 \n');
            fprintf(fid,'node 1 0.0 \n');
            fprintf(fid,'node 2 1.0 \n');
            fprintf(fid,'fix 1 1 \n');
            for i = 1:length(obj.uniaxial_material_def)
                fprintf(fid,'%s\n',obj.uniaxial_material_def{i});
            end
            fprintf(fid,'element truss 1 1 2 1.0 %i \n',obj.uniaxial_material_tag);
            fprintf(fid,'pattern Plain 1 "Series -dt 1.0 -filePath {%s} -factor 1.0" { \n',OpenSeesAnalysis.path_for_tcl(filename_disp));
            fprintf(fid,'    sp 2 1 1.0 \n');
            fprintf(fid,'} \n');
            fprintf(fid,'recorder Element -file {%s} -precision 10 -ele 1  force \n',OpenSeesAnalysis.path_for_tcl(filename_output_force));
            fprintf(fid,'recorder Element -file {%s} -precision 10 -ele 1  deformations \n',OpenSeesAnalysis.path_for_tcl(filename_output_def));
            fprintf(fid,'recorder Element -file {%s} -precision 10 -ele 1  stiff \n',OpenSeesAnalysis.path_for_tcl(filename_output_stiff));
            fprintf(fid,'system UmfPack \n');
            fprintf(fid,'constraints Penalty 1.0e12 1.0e12 \n');
            fprintf(fid,'test NormDispIncr 1.0e-8 10 0 \n');
            fprintf(fid,'algorithm Newton \n');
            fprintf(fid,'numberer RCM \n');
            fprintf(fid,'integrator LoadControl 1.0 \n');
            fprintf(fid,'analysis Static \n');
            fprintf(fid,'set ok [analyze %i] \n',numSteps);
            fprintf(fid,'if {$ok != 0} {exit 2} \n');
            fprintf(fid,'exit 1 \n');
            fclose(fid);
            
            % Run OpenSees
            [status,result] = obj.runOpenSees(filename_input);
            results.textOutput  = result;
                        
            switch status
                case 1
                    results.status = 'Analysis Successful';
                case 2
                    results.status = 'Analysis Failed';
                    fprintf('%s\n',result);
                    warning('UniaxialMaterialAnalysis:Analysis','Analysis failed\n');
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis ended in an unknown manner, exit code = %i',status);
            end            
                        
            % Read Results
            results.disp = dlmread(filename_output_def);
            temp = dlmread(filename_output_force);
            results.force = temp(:,2);
            try
                results.stiff = dlmread(filename_output_stiff);
            catch
                results.stiff = [];
            end
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_disp,filename_output_force,...
                    filename_output_def,filename_output_stiff);
            end
        end        
    end
end
