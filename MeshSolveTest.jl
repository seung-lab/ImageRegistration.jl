using PyPlot
using HDF5
using JLD

#=
data=load("r4c2-r4c3.jld")
Vertices=data["nodes"]
Incidence=data["edges"]
RestLengths=data["edge_lengths"]
Stiffnesses=data["edge_coeffs"]
Fixed=data["nodes_fixed"]
=#
data=load("r4c2-r4c3.jld")
Vertices=data["v"]
Incidence=data["e"]
RestLengths=data["l"]
Stiffnesses=data["k"]
Fixed=data["f"]
#=
fid = jldopen("./solvedMesh.jld", "w")
fid["n"] = 3317;
fid["m"] = 9781;
fid["nodes"] = Vertices[:, 1:3317];
=#

Vertices=map(Float64,Vertices)
Incidence=sparse(map(Float64,full(Incidence)))
RestLengths=map(Float64,RestLengths)
Stiffnesses=map(Float64,Stiffnesses)

d=size(Vertices,1)
V=size(Vertices,2)
E=size(Incidence,2)
Lengths=zeros(1,V)

Moving = ~Fixed

#function MeshSolve(Vertices, Incidence, Stiffnesses, RestLengths, Moving)
Moving2=[Moving[:]'; Moving[:]'][:]   # double the dimensionality

#Vertices[:,Moving]=Vertices[:,Moving]+20*randn(size(Vertices[:,Moving]))
eta=0.25
niter=250
U=zeros(1,niter)     # energy vs. time
g=similar(Vertices)  # gradient of potential energy

for iter=1:niter
    Springs=Vertices*Incidence
    g=Gradient(Springs, Incidence, Stiffnesses, RestLengths)
    if iter<10
        # gradient descent
        Vertices[:,Moving]=Vertices[:,Moving]-eta*g[:,Moving]
    else
        #  Newton's method
        #=
        H=Hessian(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta*reshape(sparse(H[Moving2,Moving2])\g[:,Moving][:],2,length(find(Moving)))
=#
        H=Hessian2(Springs, Incidence, Stiffnesses, RestLengths)
        Vertices[:,Moving]=Vertices[:,Moving]-eta*reshape(H[Moving2,Moving2]\g[:,Moving][:],2,length(find(Moving)))
    end
    U[iter]=Energy(Springs,Stiffnesses,RestLengths)
    println(iter," ", U[iter])
    #    visualize the dynamics
    #=
    subplot(221)
    cla()
    scatter(Vertices[1,:],Vertices[2,:])
    subplot(222)
    plot(1:iter,U[1:iter])
    subplot(223)
    Lengths=sqrt(sum(Springs.^2,1))
    cla()
    plot(1:E,Lengths')
    draw()
    =#
end

#=
fid["nodes_t"] = Vertices[:, 1:3317];
fid["edges"] = Incidence[1:3317, 1:9781]
=#
close(fid);


