# MeshSolve - given spring mesh, solve for equilibrium positions of vertices
# 
# V = # mesh vertices in R^d
# E = # of springs
#
# 'Vertices' - dxV matrix, columns contain vertex positions
#
# 'Incidence' - VxE generalized oriented incidence matrix
#    springs <-> columns
#    intratile spring <-> column containing 1 -1
#    intertile spring <-> column containing w1 w2 w3 -1 where (w1, w2, w3) represents a weighted combination of vertices
#
#  most functions compute Springs=Vertices*Incidences, dxE matrix
#    spring vectors <-> columns

# 'Stiffnesses', 'RestLengths' - 1xE vectors specifying spring properties

# 'Moving' - integer vector containing indices of moving vertices
# could be changed to 1xE binary vector


function Energy( Springs, Stiffnesses, RestLengths)
    # potential energy in springs
    Lengths=sqrt(sum(Springs.^2,1))   # spring lengths (row vector)
    sum(Stiffnesses[:].*(Lengths[:]-RestLengths[:]).^2)/2
end

function Gradient( Springs, Incidence, Stiffnesses, RestLengths)
    # gradient of potential energy with respect to vertex positions
    # returns dxV array, same size as Vertices
    # physically, -gradient is spring forces acting on vertices
    Forces=similar(Springs)
    Lengths=sqrt(sum(Springs.^2,1))   # need fix for divide by zero?
    for a=1:size(Springs,2)
        Forces[:,a]=Stiffnesses[a]*(1-RestLengths[a]/Lengths[a])*Springs[:,a]
    end
    Forces*Incidence'
end

function Hessian( Springs, Incidence, Stiffnesses, RestLengths)
    # Hessian of the potential energy as an Vd x Vd matrix
    # i.e. VxV block matrix of dxd blocks
    # Note: symmetric positive definite
    
    V = size(Incidence,1)
    d = size(Springs,1)
    H = zeros(V*d, V*d)

    Lengths=sqrt(sum(Springs.^2,1))

    for a=1:size(Springs,2)
        # one-step build of dH
        #        dH = (1-RestLengths[a]/Lengths[a])*eye(d)+RestLengths[a]*Springs[:,a]*Springs[:,a]'/Lengths[a]^3
        # two-step build of dH
        dH = eye(d)-Springs[:,a]*Springs[:,a]'/Lengths[a]^2  # projection perpendicular to Springs
        dH = eye(d)-RestLengths[a]/Lengths[a]*dH
        dH = Stiffnesses[a]*dH;
        VertexList=find(Incidence[:,a])    # vertices incident on spring a
        for i=VertexList
            for j=VertexList
                # indices of (i,j) block of Hessian
                irange=(i-1)*d+(1:d)
                jrange=(j-1)*d+(1:d)
                H[ irange, jrange ] += Incidence[i,a]*Incidence[j,a]*dH
            end
        end
    end  
    H
end

using PyPlot
using HDF5
using JLD

data=load("8000x8000_100px.jld")
Vertices=data["v"]
Incidence=data["e"]
RestLengths=data["l"]
Stiffnesses=data["k"]

d=size(Vertices,1)
V=size(Vertices,2)
E=size(Incidence,2)
Lengths=zeros(1,V)

#function MeshSolve(Vertices, Incidence, Stiffnesses, RestLengths, Moving)
Moving=1:7400
Moving2=[]
for i in Moving
    Moving2=[Moving2; (i-1)*d+collect(1:d)]
end
Moving2=convert(Array{Int64,1},Moving2)

Vertices[:,Moving]=Vertices[:,Moving]+30*randn(size(Vertices[:,Moving]))
eta=0.5
niter=100
U=zeros(1,niter)     # energy vs. time
g=similar(Vertices)  # gradient of potential energy

for iter=1:niter
    Springs=Vertices*Incidence
    g=Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    if iter<10
    # gradient descent
        Vertices[:,Moving]=Vertices[:,Moving]-eta*g[:,Moving]
    else
    # Newton's method
        H=Hessian(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta*reshape(sparse(H[Moving2,Moving2])\g[:,Moving][:],2,length(Moving))
    end
    # visualize the dynamics
    subplot(221)
    cla()
#    scatter(Vertices[1,:],Vertices[2,:])
    subplot(222)
    U[iter]=Energy(Springs,Stiffnesses,RestLengths)
    println(iter," ", U[iter])
    plot(1:iter,U[1:iter])
    subplot(223)
    Lengths=sqrt(sum(Springs.^2,1))
    cla()
    plot(1:E,Lengths')
    draw()
end




