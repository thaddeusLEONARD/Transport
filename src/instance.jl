
using  DataFrames,CSV,DelimitedFiles


# caracteristics of node objects - Data and solution parts
struct Vertex
    vertex_idx::Int
    demand::Int
    vertex_type::String
    lat::Float64
    long::Float64
    TW_min::Int
    TW_max::Int
    route_idx::Int        # index of the route in the solution
    order_in_route::Int    # position of the customer in the route (easy to sort for visualization)
    serv::Int # Time at which the vertex is served in the solution
end
Vertex(id,d,nt,l,L,ai,width, horizon) = Vertex(id,d,nt,l,L,ai, nt in ["S","LS"] ? ai+width : ai+horizon, -1, -1, -1)
import Base.show ;  Base.show(io::IO, n::Vertex) = println(io, "<",lpad(n.vertex_type,2,"-"),">",lpad(n.vertex_idx,2))

mutable struct Instance
    SV_cap::Int
    speed_ratio::Float64
    time_horizon::Int
    tw_width::Int
    service_duration::Int
    dist_matrix::Matrix{Float64}
    nodes::Vector{Vertex}
    D::Vector{Int}  ##depot
    P::Vector{Int}  ##parking ??
    S::Vector{Int}  ## ??
    LS::Vector{Int} ##parking ??
end

## accessor : get the distance between two nodes
getDistance(inst::Instance, v1::Vertex, v2::Vertex) = inst.dist_matrix[v1.vertex_idx,v2.vertex_idx]
getDistance(inst::Instance, i::Int, j::Int) = inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]

#function Instance(instfname, paramfname, distfname)


function parseInstance(paramf,instancef,distmatf)
    # read parameters
    # SMALL VEHICLE CAPACITY	Speed ratio small vs large	Time horizon	TW width	service duration
    param = CSV.File(paramf, delim="\t", types=[Int,Float64,Int,Int,Int]) |> DataFrame
    SV_cap, speed_ratio, time_horizon, tw_width, service_duration = param[1,:]

    # Read the instance file in a data frame ansd store it into a
    problem = CSV.File(instancef,delim='\t', types=[Int,Int,String,Float64,Float64,Int]) |> DataFrame
    nodes = [Vertex(problem[i,:]...,tw_width,time_horizon) for i = 1:size(problem, 1)]
    D, P, S, LS = (findall(x -> x.vertex_type==nodetype, nodes) for nodetype in ["D", "P", "S", "LS"])
    # duplicate the depot :
    @assert length(D) == 1
    push!(nodes, nodes[D[1]])
    push!(D, length(nodes))


    tmpi = 1
    tmp_capa = 0
    for node in nodes
        if node.vertex_type!="S"
            tmp_capa+=node.demand
        end
    end
    a =ceil(tmp_capa/SV_cap)
    a = a+1 ## le nombre de fois que je met chaque parking
    tmp_nodes = deepcopy(nodes)
    for node in tmp_nodes
        if node.vertex_type=="P"
            tmpp=[]
            for i in 1:a
                push!(tmpp,node)
            end
            splice!(nodes,tmpi:tmpi,tmpp)
            tmpi+=Int(a)
        end
        tmpi+=1
    end

    # read the distance matrix
    distanceMatrix = open(distmatf) do f ; readdlm(f) ; end
    # enforce the triangular inequality in the distance matrix

    return Instance(SV_cap, speed_ratio, time_horizon, tw_width, service_duration, distanceMatrix, nodes, D, P, S, LS),a
end
