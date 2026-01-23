    clear all; close all;

    % --- 1. SETUP & PATHS ---
    disp(' ');
    disp('%%%%%%%%%%');
    disp('Generate PCA coronary geom and run CFD');
    starttime = datetime('now');
    disp(['Started ' datestr(starttime)]);

    % Execution mode: 1 = Run CFD, 0 = Generate files only
    iexe = 1; 
   
    % ANSYS Paths
    CFX_pre_exe = 'C:\Program Files\ANSYS Inc\ANSYS Student\v252\CFX\bin\cfx5pre';
    CFX_sol_exe = 'C:\Program Files\ANSYS Inc\ANSYS Student\v252\CFX\bin\cfx5solve';
    CFX_post_exe = 'C:\Program Files\ANSYS Inc\ANSYS Student\v252\CFX\bin\cfx5post';

    % Base name of file and directry
    currentdir=pwd;
    casedir='trial';
    basename='trial';
    if exist(casedir)~=7
      mkdir(casedir);
    end

% Save Inflow data
flow_record_file = fullfile(currentdir, 'flow_rates.csv');

% Create the header if the file doesn't exist
if ~exist(flow_record_file, 'file')
    fid_log = fopen(flow_record_file, 'w');
    fprintf(fid_log, 'Filename,FlowRate_kg_s\n');
    fclose(fid_log);
end

    % --- 2. LOAD DATA ---
    
    % Load PCA Geometries
    fprintf('Loading generated geometries...\n');
    if exist('generated_arteries.mat', 'file')
        ws = load('generated_arteries.mat');
        if isfield(ws, 'generated_arteries')
            GeomCells = ws.generated_arteries;
        else
            vars = fieldnames(ws);
            GeomCells = ws.(vars{1});
        end
        nvessel = length(GeomCells);
        fprintf('Found %d vessels to process.\n', nvessel);
    else
        error('File generated_arteries.mat not found.');
    end

    % Load Template Mesh Data
    if exist('Standard_mesh_nodes.dat', 'file') && exist('Standard_mesh_node_table.dat', 'file')
        nodes = load('Standard_mesh_nodes.dat');       
        ntable = load('Standard_mesh_node_table.dat'); 
    else
        error('Standard_mesh template files missing.');
    end

    % Mesh Parameters
    npcl = 401;  % Points along centerline
    nnxc = 801;  % Points in cross section
    r00  = 1.5;  % Template radius
    Vlen = 80;   % Vessel length

    % --- 3. MAIN LOOP ---
    for ivessel = 1:nvessel
        
        fprintf('\n==================================\n');
        fprintf('Working on Model %05d ...\n', ivessel);
        
        % Initialise new mesh array
        nnew = zeros(size(nodes));
        
        % MESH GENERATION ---
        
        % Get Data & Interpolate
        raw_data = GeomCells{ivessel};
        raw_xyz = raw_data(:, 1:3);
        raw_r   = raw_data(:, 4);
        
        s_input = linspace(0, 1, size(raw_xyz, 1)); 
        s_mesh  = linspace(0, 1, npcl);    
        
        xyzc = zeros(npcl, 3);
        xyzc(:,1) = interp1(s_input, raw_xyz(:,1), s_mesh, 'pchip');
        xyzc(:,2) = interp1(s_input, raw_xyz(:,2), s_mesh, 'pchip');
        xyzc(:,3) = interp1(s_input, raw_xyz(:,3), s_mesh, 'pchip');
        rsl       = interp1(s_input, raw_r,       s_mesh, 'linear');
        
        %Morph Mesh
        tv0 = [0, 0, 1]; 
        
        for ip = 1:npcl
            z00 = Vlen * double(ip-1) / double(npcl-1);
            
            % Tangent Vector
            if ip == 1
                tv = xyzc(ip+1,:) - xyzc(ip,:);
            elseif ip == npcl
                tv = xyzc(ip,:) - xyzc(ip-1,:);
            else
                tv = xyzc(ip+1,:) - xyzc(ip-1,:);
            end
            tv = normr(tv);
            
            % Rotation
            dp = dot(tv0, tv);
            if dp > 1, dp = 1; elseif dp < -1, dp = -1; end
            theta = acos(dp);
            rvec = cross(tv0, tv);
            if norm(rvec) < 1e-6, rvec = [1, 0, 0]; else, rvec = normr(rvec); end
            
            % Transform
            ista = nnxc * (ip-1) + 1;
            ilas = nnxc * ip;
            
            slice_pts = nodes(ntable(ista:ilas), 1:3);
            slice_pts(:, 3) = slice_pts(:, 3) - z00; 
            
            % External rotation function
            rot_pts = rot_pts_3dax(slice_pts, nnxc, rvec, theta);
            
            radius_ratio = rsl(ip) / r00;
            nnew(ntable(ista:ilas), 1:3) = rot_pts * radius_ratio + xyzc(ip, :);
        end
        
        % WRITE MESH FILES ---
        caseID = sprintf('%05d', ivessel);
        fname_base = [basename caseID];
        
        clfile   = fullfile(casedir, [fname_base '_cl.dat']);
        nodefile = fullfile(casedir, [fname_base '_nodes.dat']);
        meshfile = fullfile(casedir, [fname_base '.msh']);
        
        writematrix(xyzc, clfile, 'Delimiter', ' ');
        writematrix(nnew, nodefile, 'Delimiter', ' ');
        
        % Assemble .msh
        if ispc
            cmd = sprintf('type Standard_mesh1.msh "%s" Standard_mesh2.msh > "%s"', nodefile, meshfile);
        else
            cmd = sprintf('cat Standard_mesh1.msh "%s" Standard_mesh2.msh > "%s"', nodefile, meshfile);
        end
        [s, m] = system(cmd);
        if s ~= 0, error('Mesh generation failed: %s', m); end
        
        
        % --- PREPARE CFD INPUT ---
        if iexe == 1
            
            % 1. GENERATE CCL FILE 
            
            % Calculate Flow
            Q_val = 0.5 + (1.5 - 0.5) * rand(); 
            inflow = Q_val * 1060 * 1.0e-6; % kg/s

            fid_log = fopen(flow_record_file, 'a');
            if fid_log ~= 1
                fprintf(fid_log, '%s,%e\n', fname_base, inflow);
                fclose(fid_log);
            else
                fprintf(2,'could not save flow rate for %s\n', fname_base);
            end 
            
            % Create String & Temp File
            tmpstr = sprintf('\n      Mass Flow Rate = %e [kg s^-1]\n', inflow);
            tmpfile = 'tmp.txt';
            
            fid_tmp = fopen(tmpfile, 'w');
            fprintf(fid_tmp, '%s', tmpstr);
            fclose(fid_tmp);
            
            % Sandwich Files
            cclfile = fullfile(casedir, [fname_base '.ccl']);
            if exist('Standard_CFD1.ccl','file') && exist('Standard_CFD2.ccl','file')
                if ispc
                    cmd = sprintf('type Standard_CFD1.ccl "%s" Standard_CFD2.ccl > "%s"', tmpfile, cclfile);
                else
                    cmd = sprintf('cat Standard_CFD1.ccl "%s" Standard_CFD2.ccl > "%s"', tmpfile, cclfile);
                end
                [s, m] = system(cmd);
                if s~=0, error('CCL Sandwich failed: %s', m); end
            else
                error('Standard_CFD templates missing.');
            end
            delete(tmpfile);
            
            
            % GENERATE PRE SCRIPT      
            prefile = fullfile(casedir, [fname_base '.pre']);
            deffile = fullfile(currentdir, casedir, [fname_base '.def']);
            resfile = fullfile(casedir, [fname_base '_001.res']);
            
            % Use absolute paths for the content inside the file
            fmeshfile = fullfile(currentdir, meshfile);
            fcclfile  = fullfile(currentdir, cclfile);
            
            % Open File
            fout = fopen(prefile, 'w');
            
            % Header 
           
            fprintf(fout, 'COMMAND FILE:\n  CFX Pre Version = 19.5\nEND\n');
            fprintf(fout,'\n');
            fprintf(fout, '>load mode=new\n>update\n');
            fprintf(fout,'\n');
            
            % Mesh Import
            fprintf(fout, '>gtmImport filename=%s, type=Fluent, units=mm, nameStrategy=Assembly', fmeshfile);
            fprintf(fout, '\n>update\n\n');
            
            % CCL Import
            fprintf(fout, '>importccl filename=%s, mode=replace, autoLoadLibrary=none', fcclfile);
            fprintf(fout, '\n>update\n\n');
            fprintf(fout,'\n');
            
            % Write Def
            fprintf(fout, '>writeCaseFile filename=%s, operation=write def file', deffile);
            fprintf(fout, '\n>update\n\n');
            fprintf(fout,'\n');

            % Footer
            fprintf(fout, '>close, deleteLibrary= On\n>update');
            fprintf(fout,'\n');
            fprintf(fout,'\n>quit\n>update');
            fprintf(fout,'\n');
            fclose(fout);
            
            
            % --- RUN BATCH ---
                
            % Run CFX-Pre
            if ~exist(deffile, 'file')
                fprintf('...Running CFX-Pre...\n');
                % Use absolute path for prefile in the argument
                fprefile_abs = fullfile(currentdir, prefile);
                cmd_pre = sprintf('"%s" -batch "%s"', CFX_pre_exe, fprefile_abs);
                [s, r] = system(cmd_pre);
                
                if ~exist(deffile, 'file')
                     fprintf(2, 'Pre-processing failed. Log:\n%s\n', r);
                     continue; 
                end
            end
            
            % Run CFX-Solver
            if ~exist(resfile, 'file')
                fprintf('...Running CFX-Solver...\n');
                cmd_sol = sprintf('"%s" -batch -def "%s" -par-local -part 4', CFX_sol_exe, deffile);
                system(cmd_sol);
            end
        end
    end 

%%
% =========================================================================
% 4. BATCH POST-PROCESSING LOOP (Run after all solvers finish)
% =========================================================================
fprintf('\n\n%%%%%%%% STARTING PHASE 2: BATCH POST-PROCESSING %%%%%%%%\n');

casedir = 'C:\Users\yawmi\Downloads\ML-Stenosis_Studentship-main\ML-Stenosis_Studentship-main\CFD processing\U-net\Trial';
nvint = 1;
%%% Post processing %%%
npcx=72; % number of points along circumference
for ivessel=1:nvessel

    disp('Post process')

    % Re-read centreline
    clfile = fullfile(casedir, [basename num2str(ivessel, '%05d') '_cl.dat']);
    xyzc=load(clfile);

    % Call postprocessing script 
    CFX_PostProcess(CFX_post_exe,xyzc,currentdir,...
                          [basename num2str(ivessel,'%05d')],'Wall',npcx);
 
    % Delete files (other than 1 per X cases)
    if iexe==1
        if mod(ivessel,nvint)~=0
            if ispc
                fclose('all');
                nodefile=[casedir '\' basename num2str(ivessel,'%05d') '_nodes.dat'];
                meshfile=[casedir '\' basename num2str(ivessel,'%05d') '.msh'];
                prefile= [caseir '\' basename num2str(ivessel,'%05d') '.pre'];
                cclfile= [casedir '\' basename num2str(ivessel,'%05d') '.ccl'];
                meshfile=[casedir '\' basename num2str(ivessel,'%05d') '.msh'];
                outfile= [currentdir '\' basename num2str(ivessel,'%05d') '_001.out'];
                resfile= [currentdir '\' basename num2str(ivessel,'%05d') '_001.res'];
                [status,message]=system(['del ' nodefile]);
                [status,message]=system(['del ' meshfile]);
                [status,message]=system(['del ' prefile]);
                [status,message]=system(['del ' cclfile]);
                [status,message]=system(['del ' meshfile]);
                [status,message]=system(['del ' outfile]);
                [status,message]=system(['del ' resfile]);
            else
                disp('Not a recognised system and cannot delete files');
            end
        end
    end

end


% Final message
disp(' ');
disp('Program completed');
endtime=datetime('now');
disp(datestr(endtime));
disp(['Elasped time: ' num2str(minutes(endtime-starttime)) ' minutes']);
disp('%%%%%%%%%%');
disp(' ');
