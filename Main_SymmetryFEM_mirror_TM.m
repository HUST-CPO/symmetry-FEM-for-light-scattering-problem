clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');

%read mesh data
load mirror-unit.mat
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
Mesh.Dof=Mesh.nbrVertex;
%treat mesh data
Mesh=GetEdge(Mesh);
%physic model
physic.c_const=299792458;
physic.lam0=1e-6;%wavelength
physic.k0=2*pi/physic.lam0;
%material
physic.epsilonr=ones(48,1);%Corresponding to triID
physic.epsilonr(1)=3.5*3.5;
physic.epsilonr(3:48)=3.5*3.5;
physic.mur=ones(48,1);%Corresponding to triID



%% Sub-task A
%boundary condition setting
Mesh.out=FindEdge([1:4, 6:18, 21, 22, 25, 26, 29, 30, 33, 34, 37, 38, 41:44, 47, 48, 51, 52],Mesh);
Mesh.inc=FindEdge([5],Mesh);
physic.E0=1/2;

%init matrix
solver.Ai=[];
solver.Aj=[];
solver.Av=[];

%assemble
[solver]=AssemblyOfEqu_mirror(Mesh,physic,solver);
solver.b=complex(zeros(Mesh.Dof,1));
solver=AssemblyOfInc_mirror(Mesh,physic,solver);
solver=AssemblyOfOut_mirror(Mesh,physic,solver);

%solution
solver.A=sparse(solver.Ai,solver.Aj,solver.Av);
tic
solver.x=solver.A\solver.b;
disp('sub-task A: time of solution:')
toc

Ez1=solver.x;

%% Sub-task B
%boundary condition setting
Mesh.out=FindEdge([1:4, 6:18, 21, 22, 25, 26, 29, 30, 33, 34, 37, 38, 41:44, 47, 48, 51, 52],Mesh);
Mesh.inc=FindEdge([5],Mesh);
physic.E0=1/2;
Mesh.PEC=FindIndex([19, 20, 23, 24, 27, 28, 31, 32, 35, 36, 39, 40, 45, 46, 49, 50, 53],Mesh);
Index=1:Mesh.Dof;
Index(Mesh.PEC)=[];

%init matrix
solver.Ai=[];
solver.Aj=[];
solver.Av=[];

%assemble
[solver]=AssemblyOfEqu_mirror(Mesh,physic,solver);
solver.b=complex(zeros(Mesh.Dof,1));
solver=AssemblyOfInc_mirror(Mesh,physic,solver);
solver=AssemblyOfOut_mirror(Mesh,physic,solver);

%solution
solver.A=sparse(solver.Ai,solver.Aj,solver.Av);
P=speye(Mesh.Dof);
P=P(:,Index);
solver.A=P'*solver.A*P;
solver.b=P'*solver.b;
tic
solver.x=solver.A\solver.b;
disp('sub-task B: time of solution:')
toc

Ez2=P*solver.x;


Ez=Ez1+Ez2;
%plot
PlotE(Mesh.vertex,Mesh.tri,real(Ez),0);
% PlotTri(Mesh.vertex,Mesh.tri);