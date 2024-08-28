function solver=AssemblyOfInc_mirror(mesh,physic,solver)

k0=physic.k0;

%guass points
IntegralEdgeOrder=5;
lx=[-9.061798459386640e-01,-5.384693101056831e-01, 0.000000000000000e+00, 5.384693101056831e-01, 9.061798459386640e-01];
lp=[2.369268850561890e-01, 4.786286704993665e-01, 5.688888888888889e-01, 4.786286704993665e-01,2.369268850561890e-01];

NbrNon=length(mesh.inc);
NumNonA=0;
Ai=zeros(NbrNon*4,1);
Aj=zeros(NbrNon*4,1);
Av=complex(zeros(NbrNon*4,1));

%loop of tri
for n=1:NbrNon
    numEle=mesh.inc(n,1);
    numEdg=mesh.inc(n,2);
    %coordinate
    x=zeros(5,1);y=zeros(5,1);
    x(1:3)=mesh.vertex(mesh.tri(numEle,:),1);
    y(1:3)=mesh.vertex(mesh.tri(numEle,:),2);
    if numEdg==1
        x(4)=x(1);x(5)=x(2);
        y(4)=y(1);y(5)=y(2);
    elseif numEdg==2
        x(4)=x(1);x(5)=x(3);
        y(4)=y(1);y(5)=y(3);
    else
        x(4)=x(2);x(5)=x(3);
        y(4)=y(2);y(5)=y(3);
    end
    avery=(y(4)+y(5))/2;

    phi1=(y(5)-y(4));phi2=(x(4)-x(5));
    if x(5)~=x(4)
        PhysicsFactor=(phi1)/(-phi2);
        PhysicsFactor=abs(sqrt(1+PhysicsFactor*PhysicsFactor)*(-phi2)/2);
    elseif x(5)==x(4)
        PhysicsFactor=abs(phi1)/2;
    end

    u=zeros(IntegralEdgeOrder,1);v=zeros(IntegralEdgeOrder,1);
    if numEdg==1
        for i=1:IntegralEdgeOrder
            u(i)=(lx(i)+1)/2;
            v(i)=0;
        end
    elseif numEdg==2
        for i=1:IntegralEdgeOrder
            u(i)=0;
            v(i)=(lx(i)+1)/2;
        end
    else
        for i=1:IntegralEdgeOrder
            u(i)=(lx(i)+1)/2;
            v(i)=(1-lx(i))/2;
        end
    end

    %mapping
    index=zeros(2,1);
    mappingIndex=zeros(2,1);
    if numEdg==1
        index(1)=1;index(2)=2;
        mappingIndex(1)=mesh.tri(numEle,1);
        mappingIndex(2)=mesh.tri(numEle,2);
    elseif numEdg==2
        index(1)=1;index(2)=3;
        mappingIndex(1)=mesh.tri(numEle,1);
        mappingIndex(2)=mesh.tri(numEle,3);
    else
        index(1)=2;index(2)=3;
        mappingIndex(1)=mesh.tri(numEle,2);
        mappingIndex(2)=mesh.tri(numEle,3);
    end

    %bf
    Ez=zeros(2,IntegralEdgeOrder);
    for i=1:2
        for k=1:IntegralEdgeOrder
            Ez(i,k)=BF_Ez(index(i),u(k),v(k));
        end
    end

    %material
    domain=mesh.triID(numEle);
    eps=physic.epsilonr(domain);
    refra=sqrt(eps);

    %submatrix
    Ae=complex(zeros(2,2));be=complex(zeros(2,1));
    for i=1:2
        for j=1:2
            for k=1:IntegralEdgeOrder
                Ae(i,j)=Ae(i,j)+lp(k)*PhysicsFactor*1i*k0*refra*Ez(i,k)*Ez(j,k);
            end
        end

        for k=1:IntegralEdgeOrder
            be(i)=be(i)+2*lp(k)*PhysicsFactor*1i*k0*refra*Ez(i,k)*physic.E0;
        end
    end

    %assemble
    for i=1:2
        mappingIndexi=mappingIndex(i);
        for j=1:2
            mappingIndexj=mappingIndex(j);
            NumNonA=NumNonA+1;
            Ai(NumNonA)=mappingIndexi;
            Aj(NumNonA)=mappingIndexj;
            Av(NumNonA)=Ae(i,j);
        end
        solver.b(mappingIndexi)=solver.b(mappingIndexi)+be(i);
    end

end

solver.Ai=[solver.Ai;Ai];
solver.Aj=[solver.Aj;Aj];
solver.Av=[solver.Av;Av];

end