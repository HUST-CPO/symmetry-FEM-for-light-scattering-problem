function [solver]=AssemblyOfEqu_mirror(mesh,physic,solver)

k0=physic.k0;
%guass points
[xt,yt,pt,IntegralOrder]=GetGuassPoints(4);

%init
NbrNon=mesh.nbrTri*3*3;
NumNonA=0;
Ai=zeros(NbrNon,1);
Aj=zeros(NbrNon,1);
Av=complex(zeros(NbrNon,1));

%loop tri
for n=1:mesh.nbrTri
    %coordinate of tri's vertex
    x=mesh.vertex(mesh.tri(n,:),1);
    y=mesh.vertex(mesh.tri(n,:),2);

    %jac matrix
    Jac=zeros(3,3);JacS=zeros(3,3);
    Jac(1,1)=x(2)-x(1);Jac(1,2)=y(2)-y(1);
    Jac(2,1)=x(3)-x(1);Jac(2,2)=y(3)-y(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    JacS(1,1)=InvJac(2,2);JacS(1,2)=-InvJac(2,1);
    JacS(2,1)=-InvJac(1,2);JacS(2,2)=InvJac(1,1);JacS(3,3)=1;
    DetJac=abs(det(Jac));

    %bf
    Ez=zeros(3,IntegralOrder);curlEz=zeros(3,3,IntegralOrder);
    for i=1:IntegralOrder
        for j=1:3
            Ez(j,i)=BF_Ez(j,xt(i),yt(i));
            temp=zeros(3,1);
            [temp(1),temp(2)]=BF_curlEz(j,xt(i),yt(i));
            curlEz(j,:,i)=(JacS*temp)';
        end
    end

    %material
    domain=mesh.triID(n);
    eps=physic.epsilonr(domain);

    %submatrix
    Se=complex(zeros(3,3));Te=complex(zeros(3,3));
    for i=1:3
        for j=1:3
            for k=1:IntegralOrder
                temp1=curlEz(j,:,k);
                temp2=curlEz(i,:,k);
                Se(i,j)=Se(i,j)+pt(k)*DetJac*dot(temp1,temp2);
                Te(i,j)=Te(i,j)+pt(k)*DetJac*k0*k0*eps*Ez(i,k)*Ez(j,k);
            end
        end
    end

    %mapping
    MappingIndex=zeros(3,1);
    for i=1:3
        MappingIndex(i)=mesh.tri(n,i);
    end

    %assemble
    for i=1:3
        for j=1:3
            MappIndexi=MappingIndex(i);
            MappIndexj=MappingIndex(j);

            NumNonA=NumNonA+1;
            Ai(NumNonA)=MappIndexi;
            Aj(NumNonA)=MappIndexj;
            Av(NumNonA)=Se(i,j)-Te(i,j);
        end
    end

end 


solver.Ai=[solver.Ai;Ai];
solver.Aj=[solver.Aj;Aj];
solver.Av=[solver.Av;Av];

end