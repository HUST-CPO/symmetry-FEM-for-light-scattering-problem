clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');

%C3 point group  sub-task E2   chi=exp(1i*4*pi/3)
% chi=exp(1i*4*pi/3)
%read mesh data
load rotational-unit.mat
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
%treat mesh data
Mesh=GetEdge(Mesh);
Mesh.Dof=Mesh.nbrEdges;
%boundary condition setting
Mesh.out=FindEdge([10:11],Mesh);
Mesh.SrcID=[1,4];Mesh.DstID=[6,9];
[Mesh.src,Mesh.dst,Mesh.src2dstV]=GetPeriodicBoundIndex(Mesh);

%physic model
physic.c_const=299792458;
physic.lam0=0.9e-6;%wavelength
physic.k0=2*pi/physic.lam0;
physic.Eb=[-1-1i*sqrt(3);sqrt(3)-1i;0];
physic.curlEb=zeros(3,1);
%material
physic.epsilonr=ones(2,1);%Corresponding to triID
physic.epsilonr(2)=(0.088565-1i*5.6272)^2;
physic.mur=ones(2,1);%Corresponding to triID

%init matrix
solver.Ai=[];
solver.Aj=[];
solver.Av=[];
solver.b=complex(zeros(Mesh.Dof,1));

%assemble
solver=AssemblyOfEqu_rotation(Mesh,physic,solver);
solver=AssemblyOfOut_rotation(Mesh,physic,solver);
solver=AssemblyOfBEle_rotation(Mesh,physic,solver);

%solution
solver.A=sparse(solver.Ai,solver.Aj,solver.Av);
P=speye(Mesh.Dof);
for i=1:length(Mesh.src)
    P(Mesh.dst(i),Mesh.src(i))=Mesh.src2dstV(i)*exp(1i*4*pi/3);
end
tempIndex=sort(Mesh.dst);
P(:,tempIndex)=[];
solver.A=P'*solver.A*P;
solver.b=P'*solver.b;
tic
solver.x=solver.A\solver.b;
disp('time of solution:')
toc

%post
Et=P*solver.x;
[Ex,Ey]=GetExEy(Mesh.vertex,Mesh.tri,Mesh.edgesOfTri,Et);


%plot
PlotE(Mesh.vertex,Mesh.tri,real(Ex),0);
% PlotTri(Mesh.vertex,Mesh.tri);
