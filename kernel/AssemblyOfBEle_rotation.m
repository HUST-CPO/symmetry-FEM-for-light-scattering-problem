function solver=AssemblyOfBEle_rotation(mesh,physic,solver)

k0=physic.k0;

%gauss points
[xt,yt,pt,IntegralOrder]=GetGuassPoints(2);

%loop of tri
for n=1:mesh.nbrTri
    %coordinate of tri's vertex
    xx=mesh.vertex(mesh.tri(n,:),1);
    yy=mesh.vertex(mesh.tri(n,:),2);

    %jac matrix
    Jac=zeros(3,3);
    Jac(1,1)=xx(2)-xx(1);Jac(1,2)=yy(2)-yy(1);
    Jac(2,1)=xx(3)-xx(1);Jac(2,2)=yy(3)-yy(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    JacS(1,1)=InvJac(2,2);JacS(1,2)=-InvJac(2,1);
    JacS(2,1)=-InvJac(1,2);JacS(2,2)=InvJac(1,1);JacS(3,3)=1;
    DetJac=abs(det(Jac));
    DetJacx=det(Jac);
    TJac=Jac'/DetJacx;

    %bf
    Et=zeros(3,3,IntegralOrder);curlEt=zeros(3,3,IntegralOrder);
    for i=1:IntegralOrder
        for j=1:3
            temp=zeros(3,1);
            [temp(1),temp(2)]=BF_Et(j,xt(i),yt(i));
            Et(j,:,i)=InvJac*temp;
            temp=zeros(3,1);
            temp(3)=BF_curlEt(j,xt(i),yt(i));
            curlEt(j,:,i)=TJac*temp;
        end
    end

    %ambient field
    bEt=physic.Eb;bcurlEt=physic.curlEb;

    %material
    domain=mesh.triID(n);
    epsr=physic.epsilonr(domain);

    %submatrix
    St=zeros(3,1);Tt=zeros(3,1);
    for i=1:3
        for k=1:IntegralOrder
            St(i,1)=St(i,1)+pt(k)*DetJac*sum(curlEt(i,:,k)'.*(bcurlEt));
            Tt(i,1)=Tt(i,1)+pt(k)*DetJac*k0*k0*sum(Et(i,:,k)'.*(epsr*bEt));
        end
    end

    %assemble
    for i=1:3
        MappingIndexVi=mesh.edgesOfTri(n,i);

        solver.b(MappingIndexVi,1)=solver.b(MappingIndexVi,1)+Tt(i,1)-St(i,1);
    end
end

end