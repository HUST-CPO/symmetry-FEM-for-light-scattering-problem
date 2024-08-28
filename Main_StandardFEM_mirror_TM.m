clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');

%read mesh data
load mirror-full.mat
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
Mesh.Dof=Mesh.nbrVertex;
%treat mesh data
Mesh=GetEdge(Mesh);
%boundary condition setting
Mesh.out=FindEdge([1:4, 6:72],Mesh);
Mesh.inc=FindEdge([5],Mesh);
physic.E0=1;

%physic model
physic.c_const=299792458;
physic.lam0=1e-6;%wavelength
physic.k0=2*pi/physic.lam0;
%material
physic.epsilonr=ones(86,1);%Corresponding to triID
physic.epsilonr(1)=3.5*3.5;
physic.epsilonr(3:86)=3.5*3.5;
physic.mur=ones(86,1);%Corresponding to triID

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
disp('time of solution:')
toc

%post
Ez=solver.x;

%plot
PlotE(Mesh.vertex,Mesh.tri,real(Ez),0);
% PlotTri(Mesh.vertex,Mesh.tri);