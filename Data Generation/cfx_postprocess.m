%%% Post processing function by calling CFD post
function CFX_PostProcess(CFX_pst_exe,cxyz,casedir,casename,wall_part,npcx)

%%% Input arguents
% CFX_pst_exe: executable of CFD post
% cxyz: centreline coordinate
% casedir: case directory
% casename: case name
% wall_part: name of wall part
% npcx: number of points along circumference

%%% Some standard parameters
Cunit = '[mm]'; % Centreline coordinate unit
Export = 'Yes'; % No for generating post-processing file for just visualisation of frames
CFXarg = '-batch '; % Option for CFD post

%%% Derived parameters
nslice=size(cxyz,1);
RESfile=[casedir '\' casename '_001.res'];
datadir=[casedir '\' casename '_P_ESS_XC'];
if exist(datadir,'dir')~=7
    mkdir(datadir);
end

%%% Calculation of tangential vector
tv(1,1:3)=cxyz(2,1:3)-cxyz(1,1:3);
tv(2:nslice-1,1:3)=cxyz(3:nslice,1:3)-cxyz(1:nslice-2,1:3);
tv(nslice,1:3)=cxyz(nslice,1:3)-cxyz(nslice-1,1:3);
tv=normr(tv);

%%% Calculation of inward-curvature vector (need to be carefully
%%%             re-considered for the patient-specific geometries as it may become jittery)
nv(2:nslice,1:3)=tv(2:nslice,1:3)-tv(1:nslice-1,1:3);
nv(1,1:3)=nv(2,1:3);
nv=normr(nv);

%%% Calculation of distance along centreline
ddist=zeros(nslice,1);
ddist(1)=0.0;
ddist(2:nslice)=vecnorm(cxyz(2:nslice,1:3)-cxyz(1:nslice-1,1:3),2,2);
dist=cumsum(ddist);

%%% Generate CFD-post script
ScriptFile=fullfile(casedir,['\PP_',casename,'.cse']);
fscript=fopen(ScriptFile,'w');
fprintf(fscript,'%s\n',['>load filename=',RESfile,', force_reload=\',]);
fprintf(fscript,'true\n');
fprintf(fscript, '\n>update\n');

last_file=[datadir,'/export',num2str(nslice,'%03d'),'.csv'];
for islice=1:nslice
    fprintf(fscript,'PLANE:Plane %d\n',islice);
    fprintf(fscript,'  Option = Point and Normal\n');
    fprintf(fscript,'  Point = %f %s, %f %s, %f %s\n',cxyz(islice,1),Cunit,cxyz(islice,2),Cunit,cxyz(islice,3),Cunit);
    fprintf(fscript,'  Normal = %f, %f, %f\n',tv(islice,1),tv(islice,2),tv(islice,3));
    fprintf(fscript,'  Plane Bound = Circular\n');
    fprintf(fscript,'  Bound Radius = 0.003 [m]\n');
    fprintf(fscript,'  Visibility = Off\n');
    fprintf(fscript,'END\n');
    fprintf(fscript,'\n');
    fprintf(fscript,'POLYLINE:Polyline %d\n',islice);
    fprintf(fscript,'  Option = Boundary Intersection\n');
    fprintf(fscript,'  Boundary List = %s\n',wall_part);
    fprintf(fscript,'  Location = /PLANE:Plane %d\n',islice);
    fprintf(fscript,'  Visibility = Off\n');
    fprintf(fscript,'END\n');
    fprintf(fscript,'\n');
    if(strcmp(Export,'Yes')==1)
        fprintf(fscript,'EXPORT:\n');
        fprintf(fscript,'  ANSYS Export Data = Element Heat Flux\n');
        fprintf(fscript,'  ANSYS File Format = ANSYS\n');
        fprintf(fscript,'  ANSYS Reference Temperature = 0.0 [K]\n');
        fprintf(fscript,'  ANSYS Specify Reference Temperature = Off\n');
        fprintf(fscript,'  Additional Variable List =\n');
        fprintf(fscript,'  BC Profile Type = Inlet Velocity\n');
        fprintf(fscript,'  CSV Type = CSV\n');
        fprintf(fscript,'  Case Name = Case %s\n', [casename,'_001']);
        fprintf(fscript,'  Export Connectivity = Off\n');
        fprintf(fscript,'  Export Coord Frame = Global\n');
        fprintf(fscript,'  Export File = %s\n',[datadir,'/export',num2str(islice,'%03d'),'.csv']);
        fprintf(fscript,'  Export Geometry = On\n');
        fprintf(fscript,'  Export Location Aliases=\n');
        fprintf(fscript,'  Export Node Numbers = Off\n');
        fprintf(fscript,'  Export Null Data = On\n');
        fprintf(fscript,'  Export Type = Generic\n');
        fprintf(fscript,'  Export Units System = Current\n');
        fprintf(fscript,'  Export Variable Type = Current\n');
        fprintf(fscript,'  External Export Data = None\n');
        fprintf(fscript,'  Include File Information = Off\n');
        fprintf(fscript,'  Include Header = Off\n');
        fprintf(fscript,'  Location List = /POLYLINE:Polyline %d\n',islice);
        fprintf(fscript,'  Null Token = null\n');
        fprintf(fscript,'  Overwrite = On\n');
        fprintf(fscript,'  Precision = 8\n');
        fprintf(fscript,'  Separator = \", \"\n');
        fprintf(fscript,'  Spatial Variables = X,Y,Z\n');
        fprintf(fscript,'  Variable List = Pressure, Wall Shear\n');
        fprintf(fscript,'  Vector Brackets = ()\n');
        fprintf(fscript,'  Vector Display = Scalar\n');
        fprintf(fscript,'END\n');
        fprintf(fscript,'>export\n');
        fprintf(fscript,'\n');
    end
end
fclose(fscript);

%%% Launch CFD-post for batch processing
proc=System.Diagnostics.Process;
quotedScriptFile = ['"' ScriptFile '"'];
proc.Start(CFX_pst_exe,[CFXarg quotedScriptFile]);
fprintf('...Reslicing result...\n');
iwrite=0;

maxWaitTime = 120; %Adjust this number
waitedTime = 0;    % Counter initialisation


while ~exist(last_file,'file')
    pause(1);

    % Update timer and check for timeout
    waitedTime = waitedTime + 1;
    if waitedTime > maxWaitTime
        fprintf('\nWARNING: File creation timed out after %d seconds. Skipping to next...\n', maxWaitTime);
        break; % This breaks the while loop to continue the script
    end

    dtmp = dir([datadir,'/export*.csv']);
    progress = length(dtmp)/nslice;
    if progress >0.25 && iwrite<1
        fprintf('.25%%.');
        iwrite=1;
    elseif progress >0.5 && iwrite<2
        fprintf('.50%%.');
        iwrite=2;
    elseif progress >0.75 && iwrite<3
        fprintf('.75%%.');
        iwrite=3;
    else
        fprintf('.');
    end
end
fprintf('.Completed\n');

%%% Summarise data into unfolded map
fout=fopen([casedir '\' casename '_summary.csv'],'w'); % output file
nvar=8; % assuming x,y,z,r,theta,z,p,ess
ista=2; iend=nslice-1; % exclude 1st and last slices
dmap=zeros(iend-ista+1,npcx,nvar); % prepare data container
% figure;
% plot3(cxyz(:,1),cxyz(:,2),cxyz(:,3),'k--');
% hold on;
% quiver3(cxyz(1:20:end,1),cxyz(1:20:end,2),cxyz(1:20:end,3),...
%         tv(1:20:end,1),tv(1:20:end,2),tv(1:20:end,3));
% quiver3(cxyz(1:20:end,1),cxyz(1:20:end,2),cxyz(1:20:end,3),...
%         nv(1:20:end,1),nv(1:20:end,2),nv(1:20:end,3));
for islice=ista:iend
   xcfile=[datadir,'/export',num2str(islice,'%03d'),'.csv'];
   xcdata=load(xcfile);
   xcdata(:,1:3)=1000*xcdata(:,1:3);
%    if mod(islice,20)==0
%        plot3(xcdata(:,1),xcdata(:,2),xcdata(:,3),'r-');
%        axis equal;
%    end
   rv=xcdata(:,1:3)-cxyz(islice,1:3);
   % Calculate angle from reference vector
   if islice==ista % inner curvature reference vector for the 1st slice
       nvref(1:3)=nv(islice,1:3);
   else % minimum rotation vector for the 2nd slice onwards
       nvcor=nvref*tv(islice,1:3)';
       nvref(1:3)=nvref(1:3)-nvcor*tv(islice,1:3);
       nvref=normr(nvref);
   end
   dv=cross(tv(islice,:),nvref(:));
   % Calculate angle from inner curvature reference vector
   prod=rv*nvref(1:3)'./vecnorm(rv,2,2);
   itmp=find(abs(prod)>1);
   prod(itmp)=sign(prod(itmp)); % capping acos() argument to abs(1)
   theta=acos(prod);
   % Calculate direction of theta and correct+
   drot=rv*dv';
   indx=find(drot<0);
   theta(indx)=2*pi-theta(indx);
   % Finding first point
   [val,i0]=min(abs(theta));
   % Reorder the array
   theta_r=zeros(size(theta));
   xcdata_r=zeros(size(xcdata));
   nptmp=size(theta_r,1);
   theta_r(1:nptmp-i0+1)=theta(i0:nptmp);
   xcdata_r(1:nptmp-i0+1,:)=xcdata(i0:nptmp,:);
   if i0~=1
       theta_r(nptmp-i0+2:nptmp)=theta(1:i0-1);
       xcdata_r(nptmp-i0+2:nptmp,:)=xcdata(1:i0-1,:);
   end
   % Eliminate overlap
   dtheta=zeros(size(theta_r));
   dtheta(2:end)=theta_r(2:end)-theta_r(1:end-1);
   dtheta(1)=1e+6;
   indx_nonzero=find(abs(dtheta)>1e-8);
   theta_rsm=theta_r(indx_nonzero);
   xcdata_rsm=xcdata_r(indx_nonzero,:);
   % Interpolate data
   splf=csape(theta_rsm',xcdata_rsm','periodic');
   theta_intp=linspace(0,2*pi,npcx+1)';
   xcdata_intp=fnval(splf,theta_intp)';
%   nindx=size(xcdata_rsm,2);
%   for indx=1:nindx
%       xcdata_intp(:,indx)=interp1(theta_rsm,xcdata_rsm(:,indx),theta_intp,'spline');
%   end
   % Now recalculate the angle theta from the inner-curvature vector
   rv_intp=xcdata_intp(:,1:3)-cxyz(islice,1:3);
   dv_crv=cross(tv(islice,:),nv(islice,:));
   theta_crv=acos(rv_intp*nv(islice,1:3)'./vecnorm(rv_intp,2,2));
   drot_crv=rv_intp*dv_crv';
   indx_crv=find(drot_crv<0);
   theta_crv(indx_crv)=-theta_crv(indx_crv); % this time -pi to +pi
   % Storing the data in global data matrix
   dmap(islice,1:npcx,1:3)=xcdata_intp(1:npcx,1:3); % xyz
   dmap(islice,1:npcx,4)=vecnorm(rv_intp(1:npcx),2,2)'; % r
%   dmap(islice,1:npcx,5)=theta_intp(1:npcx); % theta
   dmap(islice,1:npcx,5)=theta_crv(1:npcx); % theta wrt inner curvature
   dmap(islice,1:npcx,6)=dist(islice)*ones(npcx,1)'; % z
   dmap(islice,1:npcx,7)=xcdata_intp(1:npcx,4); % pressure
   dmap(islice,1:npcx,8)=xcdata_intp(1:npcx,5); % ess
   for ipcx=1:npcx
       fprintf(fout,'%e %e %e %e %e %e %e %d\n',dmap(islice,ipcx,:));
   end
   
%    if islice==nslice-1
%        figure; plot(theta*360/(2*pi));
%        hold on, plot(theta_r*360/(2*pi));
%        plot(theta_rsm*360/(2*pi));
%        figure; plot3(xcdata(:,1),xcdata(:,2),xcdata(:,3),'ko');
%        hold on, plot3(xcdata_intp(:,1),xcdata_intp(:,2),xcdata_intp(:,3),'rx-');
%        figure; plot(theta,xcdata(:,4),'ko');
%        hold on, plot(theta_intp,xcdata_intp(:,4),'rx-');
%    end
end
fclose(fout);
% clear the data folder
if ispc
    close('all');
    [status,message]=system(['rmdir /s /q ' datadir]);
    [status,message]=system(['del ' ScriptFile]);
else
    disp('Not supported system and data folder not deleted');
end

end
