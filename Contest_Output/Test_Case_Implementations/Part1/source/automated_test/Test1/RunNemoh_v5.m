%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     %
% Full single body Nemoh execution from .GDF geometry %
%                                                     %
% by Morten Th�tt Andersen,                           %
% Aalbrog University 2014,                            %
% mta@civil.aau.dk                                    %
%                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Requirements:
%   Nemoh directory containing:
%       mesh.exe
%       preProcessor.exe
%       Solver.exe
%       postProcessor.exe
%
%   Project directory containing:
%       *.gdf (e.g. geometry exported from Rhino)
%       <raw nemoh input and output files will be stored here>
%
%   Output directory:
%       <structured output file will be saved here>
%
%


%% CLEAR
clear; close all; clc

%% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMPLE - Freely setup your desired run parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ID='SingleCyl_6dof_100frq_3beta'; % Indentifier for run
%geometry='singleCylinder'; % Filename of *.gdf (without extension)

%binit92 -changes here 
geometry='Test1/flap-meshed-quad';

path.project='.'; % Path of *.gdf
path.nemoh='F:\Codes\Nemoh\Release'; % Path of Nemoh executables
path.output='.'; % Output folder

h=0;                     % Water depth (0 for deep water) [m]
rho=1000;                % Fluid density [kg/m^3]
frq=[100 0.05 3.0*2*pi]; % Wave frequencies {number | min | max} [rad/s]
beta=[3 0 60];           % Wave directions {number | min | max} [deg]

dofs=[1 1 1 1 1 1]; % On/Off - {surge | sway | heave | roll | pitch | yaw}
CoG=[0 0 -0.037];   % Center of gravity [m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADVANCED USER INPUTS - No changes required to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IRF=[1 0.1 10];  % Calculation of IRF { on/off | timestep | duration} [s]
                 % Note: The IRF output file also provides the added mass at infinity

tasks=[1 1 1];   % On/off - {Setup Nemoh input files | Execute Nemoh | Structure output}
                 % Note: Only the Nemoh execution is time consuming

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAUTION! - Don't make changes below this point unless strictly required!
% Changes may cause instability in results and/or data sorting!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP Nemoh input files
if tasks(1)
    
% READ *.gdf and WRITE temp nemoh mesh *_nemoh
disp('TASK 1.1: Read *.gdf')
cd(path.project)
if exist('Mesh','dir')~=7; mkdir('Mesh'); end
gdfmesh = [geometry,'.gdf'];
nemohmesh = ['Mesh',filesep,geometry,'_nemoh'];

fid = fopen(gdfmesh);
header=fgetl(fid);
temp=fscanf(fid,'%*f %f ULEN GRAV \n %d %*d ISX ISY \n %d \n',[3,1]);
g=temp(1); Xsym=temp(2); N=temp(3);
A = fscanf(fid,'%g %g %g', [3 inf])';
fclose(fid);

fid = fopen(nemohmesh,'w');
fprintf(fid,'%d\r\n',N*4);
fprintf(fid,'%d\r\n',N);
fprintf(fid,'%0.7f %0.7f %0.7f\n',A');
fprintf(fid,'%d %d %d %d\n',[1:4:length(A);2:4:length(A);3:4:length(A);4:4:length(A)]);
fclose(fid);

% WRITE Mesh.cal
disp('TASK 1.2: Write Mesh.cal')
fid=fopen('Mesh.cal','w');
fprintf(fid,[geometry,'_nemoh','\r\n'],1);
fprintf(fid,'%g \r\n0. 0.\r\n',Xsym);
fprintf(fid,'%g %g %g \r\n',CoG);
fprintf(fid,'%g \r\n2\r\n0. \r\n1.\r\n',N);
fclose(fid);

% WRITE Nemoh.cal
disp('TASK 1.3: Write Nemoh.cal')
fid=fopen('Nemoh.cal','w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \r\n');
fprintf(fid,'%g				! RHO 			! kg/m^3 	! Fluid density \r\n',rho);
fprintf(fid,'%g				! G			! m/s^2	! Gravity \r\n',g);
fprintf(fid,'%g                 ! DEPTH			! m		! Water depth \r\n',h);
fprintf(fid,'0.	0.              ! XEFF YEFF		! m		! Wave measurement point\r\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\r\n');
fprintf(fid,'1				! Number of bodies\r\n');
fprintf(fid,'--- Body 1 -----------------------------------------------------------------------------------------------------------------------\r\n');
fprintf(fid,'%s.dat         ! Name of mesh file \r\n',nemohmesh);
fprintf(fid,'%g %g			! Number of points and number of panels 	\r\n',N*4,N);
fprintf(fid,'%g				! Number of degrees of freedom\r\n',length(nonzeros(dofs)));
if dofs(1); fprintf(fid,'1 1. 0. 0. 0. 0. 0.		! Surge\r\n'); end;
if dofs(2); fprintf(fid,'1 0. 1. 0. 0. 0. 0.		! Sway\r\n'); end;
if dofs(3); fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Heave\r\n'); end;
if dofs(4); fprintf(fid,'2 1. 0. 0. %g %g %g		! Roll about a point\r\n',CoG); end;
if dofs(5); fprintf(fid,'2 0. 1. 0. %g %g %g		! Pitch about a point\r\n',CoG); end;
if dofs(6); fprintf(fid,'2 0. 0. 1. %g %g %g		! Yaw about a point\r\n',CoG); end;
fprintf(fid,'%g				! Number of resulting generalised forces\r\n',length(nonzeros(dofs)));
if dofs(1); fprintf(fid,'1 1. 0. 0. 0. 0. 0.		! Force in x direction\r\n'); end;
if dofs(2); fprintf(fid,'1 0. 1. 0. 0. 0. 0.		! Force in y direction\r\n'); end;
if dofs(3); fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Force in z direction\r\n'); end;
if dofs(4); fprintf(fid,'2 1. 0. 0. %g %g %g		! Moment force in x direction about a point\r\n',CoG); end;
if dofs(5); fprintf(fid,'2 0. 1. 0. %g %g %g		! Moment force in y direction about a point\r\n',CoG); end;
if dofs(6); fprintf(fid,'2 0. 0. 1. %g %g %g		! Moment force in z direction about a point\r\n',CoG); end;
fprintf(fid,'0				! Number of lines of additional information \r\n');
fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\r\n');
fprintf(fid,'%g %g %g		! Number of wave frequencies, Min, and Max (rad/s)\r\n',frq);
fprintf(fid,'%g %g %g		! Number of wave directions, Min and Max (degrees)\r\n',beta);
fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\r\n');
fprintf(fid,'%g %g %g		! IRF 				! IRF calculation (0 for no calculation), time step and duration\r\n',IRF);
fprintf(fid,'0				! Show pressure\r\n');
fprintf(fid,'0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\r\n');
fprintf(fid,'0	2	100.	100.	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction');
fclose(fid);

% WRITE input.txt
disp('TASK 1.4: Write input.txt')
fid=fopen('input.txt','w');
fprintf(fid,'--- Calculation parameters ------------------------------------------------------------------------------------\r\n');
fprintf(fid,'0			! Indiq_solver	! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)\r\n');
fprintf(fid,'20			! IRES		! Restart parameter for GMRES\r\n');
fprintf(fid,'5.E-07			! TOL_GMRES	! Stopping criterion for GMRES\r\n');
fprintf(fid,'100			! MAXIT		! Maximum iterations for GMRES\r\n');
fprintf(fid,'1			! Sav_potential	! Save potential for visualization\r\n');
fclose(fid);

% WRITE ID.dat
disp('TASK 1.5: Write ID.dat') % Printed as dummy option
fid=fopen('ID.dat','w');
fprintf(fid,'%g \r\n',1);
fprintf(fid,'%s \r\n','.');
fclose(fid);

clear line N A gdfmesh nemohmesh ULEN Xsym Ysym ans

disp('TASK 1 DONE')
end

%% EXECUTE Nemoh
if tasks(2)

cd(path.project)
% RUN Mesh.exe
% Mesh for solver will be *_nemoh.dat
disp('TASK 2.1: Mesh.exe')
system([path.nemoh,filesep,'Mesh.exe']);

% RUN preProcessor.exe
disp('TASK 2.2: preProcessor.exe')
if exist('Results','dir')~=7; mkdir('Results'); end
system([path.nemoh,filesep,'preProcessor.exe']);

% RUN Solver.exe
disp('TASK 2.3: Solver.exe')
system([path.nemoh,filesep,'Solver.exe']);

% RUN postProcessor.exe
disp('TASK 2.4: postProcessor.exe')
system([path.nemoh,filesep,'postProcessor.exe']);

disp('TASK 2 DONE')
end

%% STRUCTURE output
if tasks(3)

% TRIM workspace
disp('TASK 3.1: Trim workspace')
nemoh.info.ID=ID; clear ID
nemoh.info.geometry=geometry; clear geometry
nemoh.info.path=path; clear path
nemoh.info.dofs=dofs; clear dofs
nemoh.info.CoG=CoG; clear CoG
nemoh.info.h=h; clear h
nemoh.info.rho=rho; clear rho
if exist('g','var'); nemoh.info.g=g; clear g;
else warning('g undefined! Run task 1.1'); end

nemoh.ndof=length(nonzeros(nemoh.info.dofs));
nemoh.w=linspace(frq(2),frq(3),frq(1)); clear frq
nemoh.T=nemoh.w(end:-1:1)/2*pi;
nemoh.beta=linspace(beta(2),beta(3),beta(1)); clear beta

% LOAD
disp('TASK 3.2: Load results')

% KH.dat - Hydrostatic stiffness
filename=[nemoh.info.path.project,filesep,'Mesh',filesep,'KH.dat'];
if exist(filename,'file')==2
    nemoh.Khyd=load(filename);
else warning('KH.dat missing! Run task 2.1'); end
clear filename

% Hydrostatics.dat - CoB, displacement, waterplane area
filename=[nemoh.info.path.project,filesep,'Mesh',filesep,'Hydrostatics.dat'];
if exist(filename,'file')==2
    fid=fopen(filename);
    nemoh.info.CoB=fscanf(fid,' XF = %g - XG = %*g \n YF = %g - YG = %*g \n ZF = %g - ZG = %*g \n')';
    nemoh.info.displacement=fscanf(fid,' Displacement = %f \n');
    nemoh.info.waterplaneArea=fscanf(fid,' Waterplane area = %f \n');
    fclose(fid);
    if length(fieldnames(nemoh.info))==11
        nemoh.info=orderfields(nemoh.info,[1 2 3 4 5 9 6 7 8 10 11]);
    end
else warning('Hydrostatics.dat missing! Run task 2.1'); end
clear filename

% CA.dat - Damping
filename=[nemoh.info.path.project,filesep,'Results',filesep,'CA.dat'];
if exist(filename,'file')==2
    nemoh.damp=zeros(nemoh.ndof,nemoh.ndof,length(nemoh.w));
    fid = fopen(filename);
    Np = fscanf(fid,'Nb de periode :  %d \n');
    w = zeros(Np,1);
    for i=1:Np
        w(i) = fscanf(fid,'%g\n',[1,1]);
        temp = fscanf(fid,'%g',[nemoh.ndof nemoh.ndof])';
        nemoh.damp(:,:,i)=temp;
    end
    fclose(fid);
else warning('CA.dat missing! Run task 2.4'); end
clear filename Np i w temp

% CM.dat - Added mass
filename=[nemoh.info.path.project,filesep,'Results',filesep,'CM.dat'];
if exist(filename,'file')==2
    nemoh.addedmass=zeros(nemoh.ndof,nemoh.ndof,length(nemoh.w));
    fid = fopen(filename);
    Np = fscanf(fid,'Nb de periode :  %d \n');
    w = zeros(Np,1);
    for i=1:Np
        w(i) = fscanf(fid,'%g\n',[1,1]);
        temp = fscanf(fid,'%g',[nemoh.ndof nemoh.ndof])';
        nemoh.addedmass(:,:,i)=temp;
    end
    fclose(fid);
else warning('CM.dat missing! Run task 2.4'); end
clear filename Np i w temp

% IRF.tec - Added mass at inf freq, IRF for radiation
filename=[nemoh.info.path.project,filesep,'Results',filesep,'IRF.tec'];
if exist(filename,'file')==2
    nemoh.infaddedmass=zeros(nemoh.ndof,nemoh.ndof);
    nemoh.IRF=zeros(nemoh.ndof,nemoh.ndof,IRF(3)/IRF(2));
    fid = fopen(filename);
    for i=1:nemoh.ndof+1; header=fgetl(fid); end
    for i=1:nemoh.ndof
        seperator=fgetl(fid);
        temp=fscanf(fid,'%g \n',[nemoh.ndof*2+1,IRF(3)/IRF(2)])';
        nemoh.infaddedmass(i,:)=temp(i,2:2:end);
        nemoh.IRF(:,i,:)=temp(:,3:2:end)';
    end
    nemoh.IRFt=temp(:,1)';
    fclose(fid);
else warning('IRF.tec missing! Run task 2.4'); end
clear filename i header seperator temp IRF

% ExcitationForce.tec - Excitation force
filename=[nemoh.info.path.project,filesep,'Results',filesep,'ExcitationForce.tec'];
if exist(filename,'file')==2
    nemoh.excForce.amp=zeros(length(nemoh.beta),length(nemoh.w),nemoh.ndof);
    nemoh.excForce.phase=zeros(length(nemoh.beta),length(nemoh.w),nemoh.ndof);
    fid = fopen(filename);
    for i=1:nemoh.ndof+1; header=fgetl(fid); end
    for i=1:length(nemoh.beta)
        seperator=fgetl(fid);
        temp=fscanf(fid,'%g \n',[nemoh.ndof*2+1,length(nemoh.w)])';
        nemoh.excForce.amp(i,:,:)=temp(:,2:2:end);
        nemoh.excForce.phase(i,:,:)=temp(:,3:2:end);
    end
    fclose(fid);
else warning('ExcitationForce.tec missing! Run task 2.4'); end
clear filename i header seperator temp

% DiffractionForce.tec - Diffraction force
filename=[nemoh.info.path.project,filesep,'Results',filesep,'DiffractionForce.tec'];
if exist(filename,'file')==2
    nemoh.diffForce.amp=zeros(length(nemoh.beta),length(nemoh.w),nemoh.ndof);
    nemoh.diffForce.phase=zeros(length(nemoh.beta),length(nemoh.w),nemoh.ndof);
    fid = fopen(filename);
    for i=1:nemoh.ndof+1; header=fgetl(fid); end
    for i=1:length(nemoh.beta)
        seperator=fgetl(fid);
        temp=fscanf(fid,'%g \n',[nemoh.ndof*2+1,length(nemoh.w)])';
        nemoh.diffForce.amp(i,:,:)=temp(:,2:2:end);
        nemoh.diffForce.phase(i,:,:)=temp(:,3:2:end);
    end
    fclose(fid);
else warning('DiffractionForce.tec missing! Run task 2.4'); end
clear filename i header seperator temp

% FKForce.tec - Froude-Krylov force
filename=[nemoh.info.path.project,filesep,'Results',filesep,'FKForce.tec'];
if exist(filename,'file')==2
    nemoh.FKForce.amp=zeros(length(nemoh.beta),length(nemoh.w),nemoh.ndof);
    nemoh.FKForce.phase=zeros(length(nemoh.beta),length(nemoh.w),nemoh.ndof);
    fid = fopen(filename);
    for i=1:nemoh.ndof+1; header=fgetl(fid); end
    for i=1:length(nemoh.beta)
        seperator=fgetl(fid);
        temp=fscanf(fid,'%g \n',[nemoh.ndof*2+1,length(nemoh.w)])';
        nemoh.FKForce.amp(i,:,:)=temp(:,2:2:end);
        nemoh.FKForce.phase(i,:,:)=temp(:,3:2:end);
    end
    fclose(fid);
else warning('FKForce.tec.tec missing! Run task 2.2'); end
clear filename i header seperator temp

clear ans fid tasks
cd(nemoh.info.path.output);
save([nemoh.info.ID,'_nemoh']);

disp('TASK 3 DONE')
end

