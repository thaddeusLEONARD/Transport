# Using CBC https://github.com/JuliaOpt/Cbc.jl
using Cbc
#using Gurobi
using JuMP,JSON
using DelimitedFiles

#using GLPKMathProgInterface
#using GLPK
include("instance.jl")

function cycle_solved(m, x,n)
    N = size(x)[1]
    x_val = JuMP.value.(x)
    # find cycle
    pris_i = [0 for i in 1:N]
    cycle_idx = Int[]
    push!(cycle_idx, 1)
    while true
        v, idx = findmax(x_val[cycle_idx[end],1:N])

        if idx == cycle_idx[1] || idx == n
            push!(cycle_idx,idx)
            break
        else
            push!(cycle_idx,idx)
        end
    end

    return cycle_idx
end

function get_ensemble(inst::Instance)
    n = size(inst.nodes)[1]
    Jl = Int[]
    Jp = Int[]
    Js = Int[]
    J = Int[]
    Jnonp = Int[]
    Jnonp2 = Int[]
    JnonpL = Int[]
    for i in 1:n
        if inst.nodes[i].vertex_type!="S"
            append!(Jl,[i])
            if inst.nodes[i].vertex_type!="P"
                append!(JnonpL,[i])
            end
        end
        if inst.nodes[i].vertex_type=="P"
            append!(Jp,[i])
        else
            append!(Jnonp2,[i])
            if inst.nodes[i].vertex_type!="D"
                append!(Jnonp,[i])
            end
        end
        if inst.nodes[i].vertex_type=="S"
            append!(Js,[i])
        end
        if inst.nodes[i].vertex_type!="D"
            append!(J,[i])
        end
    end
    return Jl,Js,Jp,Jnonp,Jnonp2,JnonpL,J
end

function get_parking(inst::Instance,Jp,a)
    n = size(inst.nodes)[1]

    Jpp = [[] for i in 1:a]
    tmpii=1
    i=Jp[1]
    while i < Jp[length(Jp)]
        tmpi=i
        #println(tmpii)
        #println(Jpp)
        while tmpi<= inst.nodes[i].vertex_idx==inst.nodes[tmpi].vertex_idx
            push!(Jpp[tmpii],tmpi)
            tmpi+=1
        end
        i=tmpi
        tmpii+=1
    end
    return Jpp
end

function get_max(inst::Instance)
    n = size(inst.nodes)[1]
    tmp_max = 0
    for i in 1:n
            if inst.nodes[i].TW_min+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[n].vertex_idx]>tmp_max
                tmp_max=inst.nodes[i].TW_min+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[n].vertex_idx]
            end
    end

    return tmp_max
end


function get_min(inst::Instance)
    n = size(inst.nodes)[1]
    T = inst.time_horizon
    tmp_min = [Float64(T) for i in 1:n]
    #println(tmp_min)
    for i in 1:n
        for j in 1:n
            if inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]>0 && inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]<tmp_min[i]
                #println(i)
                tmp_min[i] = inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]
           end
        end

    end

    return tmp_min
end


# buildTSP
# Given a matrix of city locations, build the TSP
# Inputs:
#   n       Number of cities
#   durations  n-by-n matrix of travel times between cities
# Output:
#   m       JuMP model

function buildTSP(inst::Instance,a)
    model = Model(with_optimizer(Cbc.Optimizer, logLevel=0))
    #model = Model(with_optimizer(Gurobi.Optimizer, Presolve=0, OutputFlag=0))
    n = size(inst.nodes)[1]
    T = inst.time_horizon
    TT = 23400
    tmp_max =get_max(inst)
    Q = inst.SV_cap
    tmp_min = get_min(inst)

    Jl,Js,Jp , Jnonp,Jnonp2,JnonpL, J= get_ensemble(inst)
    Jpp    =get_parking(inst,Jp,a)

    @variable(model, x[1:n,1:n] ,Bin);
    @variable(model, y[1:n,1:n] ,Bin);
    @variable(model, z[1:n,1:n] ,Bin);
    @variable(model, T>=qL[1:n]>=0 ); #
    @variable(model, T>=qS[1:n]>=0 ); #
    @variable(model, Q>=Cs[1:n]>=0 ); #
    @variable(model, Re[1:n] ,Bin);
    # constraints -- satisfy the demand exactly
    #@constraint(model, Demand[j=1:nf], sum(m[i,j]*x[i] for i=1:n) == D[j]);

    # objective -- minimize total waste
    @objective(model, Min, sum(sum(x[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx])+y[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio)-z[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio) for i=1:n) for j in 1:n));


    # contraintes de depot entrée sortie 1fois et flots
    #@constraint(model, Cs[1]==Q)
    @constraint(model, qL[n]<=T)
    @constraint(model, qS[n]<=T)
    ## le x peut pas aller dans un Js
    @constraint(model,[i = 1:n,j = Js,i!=j], x[i,j] == 0);

    ## on peut pas aller sur soi meme
    @constraint(model,[i = 1:n], x[i,i] == 0);
    @constraint(model,[i = 1:n], y[i,i] == 0);
    # on sort un fois de 1 et on rentre un fois dans n
    @constraint(model, sum(x[1,j] for j in J) ==1 );
    @constraint(model, sum(y[1,j] for j in J) ==1 );
    @constraint(model, sum(x[j,n] for j in J) ==1 );
    @constraint(model, sum(y[j,n] for j in J) ==1 );
    @constraint(model, sum(z[1,j] for j in J) ==1 );
    @constraint(model, sum(z[j,n] for j in J) ==1 );
    #@constraint(model, sum(sum(y[i,j] for j in 1:n) for i in 1:n) >= length(Js)+3);
    # on rentre jamais en 1 et on sort jamais de n
    @constraint(model, sum(x[n,j] for j in 1+1:n) ==0 );
    @constraint(model, sum(y[n,j] for j in 1+1:n) ==0 );
    @constraint(model, sum(x[j,1] for j in 1:n-1) ==0 );
    @constraint(model, sum(y[j,1] for j in 1:n-1) ==0 );

    for jpp in Jpp
        @constraint(model,[i = jpp,j = jpp], y[i,j] == 0);
        @constraint(model,[i = jpp,j = jpp], x[i,j] == 0);
        for i in jpp
            @constraint(model,[i2 in jpp; i<i2],qL[i]-inst.service_duration-tmp_min[i]*2>=qL[i2]);
            @constraint(model,[i2 in jpp; i<i2],qS[i]-inst.service_duration-tmp_min[i]*2>=qS[i2]);
        end
    end


    # pour tout i,j in J on rentre autant de fois dans un noeud qu'on en sort
    @constraint(model,[i = J], sum(x[i,j] for j in 1:n) == sum(x[j,i] for j in 1:n) );
    @constraint(model,[i = J], sum(y[i,j] for j in 1:n) == sum(y[j,i] for j in 1:n) );
    @constraint(model,[i = Jnonp], sum(z[i,j] for j in 1:n) == sum(z[j,i] for j in 1:n));



    @constraint(model,[i = 1:n], sum(x[i,j] for j in 1:n)<=1)
    @constraint(model,[j = 1:n], sum(x[i,j] for i in 1:n)<=1)
    #@constraint(model,[i = Jl], sum(y[i,j] for j in 1:n)<=1)
    #@constraint(model,[j = Jl], sum(y[i,j] for i in 1:n)<=1)
    @constraint(model,[j = Js], sum(y[i,j] for i in 1:n)==1)
    @constraint(model,[i = Js], sum(y[i,j] for j in 1:n)==1)
    #@constraint(model,[i = 1:n], sum(z[i,j] for j in 1:n)<=1)
    #@constraint(model,[j = 1:n], sum(z[i,j] for i in 1:n)<=1)
    ## tu rentres exactement une fois dans chaque Jnonp
    @constraint(model,[j = Jnonp], sum(x[i,j] +y[i,j]-z[i,j] for i in 1:n j!=i)   == 1);


    #def de z
    @constraint(model,[i = Jl,j = Jl,i!=j], z[i,j] <= (x[i,j]+y[i,j])/2);
    @constraint(model,[j = Jnonp], sum(x[i,j] for i in 1:n j!=i) >= sum(z[i,j] for i in 1:n j!=i));
    #@constraint(model,[j = Js], sum(z[i,j] for i in 1:n)==0)
    #@constraint(model,[i = Js], sum(z[i,j] for j in 1:n)==0)

    #def de rechargement
    @constraint(model,[i = Jp], Re[i] <= sum(x[j,i] for j in 1:n) - (qS[i]-qL[i])/TT) ;
    @constraint(model,[i = Jnonp2], Re[i] == 0);
    ## demande de i vers j pour yij
    @constraint(model,[i=1:n,j = Jnonp], Cs[j] <= Cs[i]-inst.nodes[j].demand*(1-x[i,j])+Q*(1-y[i,j]));
    ## ravitaillement des yij
    @constraint(model,[j = Jp], Cs[j] >= Q*Re[j]);

    #time windows
    #@constraint(model,[j = 1:n,i = 1:n,i!=j],  qL[j]+T*(1-x[i,j])>=qL[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx] );
    @constraint(model,[j = 1:n,i = 1:n-1,i!=j],  qL[j]+TT*(1-x[i,j])>=qL[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx] );
    @constraint(model,[j = 1:n,i = 1:n-1,i!=j],  qS[j]+TT*(1-y[i,j])>=qS[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio );
    @constraint(model,[j = 1:n],  qL[j]>=inst.service_duration+inst.nodes[j].TW_min)
    @constraint(model,[j = 1:n],  qS[j]>=inst.service_duration+inst.nodes[j].TW_min)
    @constraint(model,[j = 1:n],  qL[j]<=inst.nodes[j].TW_max)
    @constraint(model,[j = 1:n],  qS[j]<=inst.nodes[j].TW_max)


    @time optimize!(model)
    #=
    x_val = JuMP.value.(x)
    for i=1:n
        for j=1:n
            if x_val[i,j]>0
                println(" i,j : ",i,",",j," val_x : ",x_val[i,j])
            end
        end

    end
    x_val = JuMP.value.(y)

    for i=1:n
        for j=1:n
            if x_val[i,j]>0
                println(" i,j : ",i,",",j," val_y : ",x_val[i,j])
            end
        end
    end


    x_val = JuMP.value.(z)

    for i=1:n
        for j=1:n
            if x_val[i,j]>0
                println(" i,j : ",i,",",j," val_z : ",x_val[i,j])
            end
        end

    end
    =#
    #=cs = JuMP.value.(Cs)
    qs = JuMP.value.(qS)
    ql = JuMP.value.(qL)
    for i=1:n
        println(" i : ",i," cs : ",cs[i])
        println(" i : ",i," qs : ",qs[i])
        println(" i : ",i," ql : ",ql[i])

    end=#
    if termination_status(model)==MOI.INFEASIBLE
        println("INFEASIBLE")
        return [],[],0
    else
        cyclex = cycle_solved(model,x,n)
        cycley = cycle_solved(model,y,n)
        return cyclex,cycley, objective_value(model)
    end
end # end buildTSP




function main()
    # Main program starting here

    # PARSE iNstance, parameters and distanceMatrix
    instdir = "..\\instances2019\\"
    outdir = "..\\output\\"
    instNames= ["C1-2-8.txt","C2-2-8.txt","C1-3-10.txt","C1-3-12.txt","C2-3-10.txt","C2-3-12.txt",
          "R1-2-8.txt","R2-2-8.txt","R1-3-10.txt","R1-3-12.txt","R2-3-10.txt","R2-3-12.txt"]

    paramf = string(instdir,"parameters.txt")


        for i in 7:7
            instf =string(instdir,instNames[i])
            matf =string(instdir,"distancematrix98.txt")
            inst,a = parseInstance(paramf,instf,matf)
            outsvc=string(i,"svc.csv")
            outspeed=string(i,"speed.csv")
            outhorizon=string(i,"horizon.csv")
            outwidth=string(i,"width.csv")
            outduration=string(i,"duration.csv")
            outRes=string(i,"res.csv")
            #parametre initiaux
            debSvc=inst.SV_cap*1.5
            debSpeed=inst.speed_ratio*1.5
            debHorizon=inst.time_horizon*1.5
            debWidth=inst.tw_width*1.5
            debService=inst.service_duration*1.5
            svc=[]
            speed=[]
            horizon=[]
            width=[]
            duration=[]
            #=for j in 1:11
                #maj parameters
                inst,a = parseInstance(paramf,instf,matf)
                inst.SV_cap=debSvc-floor((j-1)/10*debSvc)
                deb=time()
                @time  cyclex,cycley,opt= buildTSP(inst,a)
                fin=time()-deb
                push!(svc,[opt,fin])

                inst,a = parseInstance(paramf,instf,matf)
                inst.speed_ratio=debSpeed-floor((j-1)/10*debSpeed)
                deb=time()
                @time  cyclex,cycley,opt= buildTSP(inst,a)
                fin=time()-deb
                push!(speed,[opt,fin])

                inst,a = parseInstance(paramf,instf,matf)
                inst.time_horizon=debHorizon-floor((j-1)/10*debHorizon)
                deb=time()
                @time  cyclex,cycley,opt= buildTSP(inst,a)
                fin=time()-deb
                push!(horizon,[opt,fin])

                inst,a = parseInstance(paramf,instf,matf)
                inst.tw_width=debWidth-floor((j-1)/10*debWidth)
                deb=time()
                @time  cyclex,cycley,opt= buildTSP(inst,a)
                fin=time()-deb
                push!(width,[opt,fin])

                inst,a = parseInstance(paramf,instf,matf)
                inst.service_duration=debService-floor((j-1)/10*debService)
                deb=time()
                @time  cyclex,cycley,opt= buildTSP(inst,a)
                fin=time()-deb
                push!(duration,[opt,fin])

            end
            writedlm(outsvc,svc, ", ")
            writedlm(outspeed,speed, ", ")
            writedlm(outhorizon,horizon, ", ")
            writedlm(outwidth,width, ", ")
            writedlm(outduration,duration, ", ")=#
            inst,a = parseInstance(paramf,instf,matf)
            deb=time()
            @time  cyclex,cycley,opt= buildTSP(inst,a)
            fin=time()-deb
            CSV.write(string(i,"resx.csv"), inst.nodes[cyclex],delim='\t')
            CSV.write(string(i,"resy.csv"), inst.nodes[cycley],delim='\t')
            # saving some results to be loaded somewhere else:
            #saving in a CSV File
            CSV.write(string(outdir,i,"routes.csv"), inst.nodes,delim='\t')
            CSV.write(string(outdir,i,"routesx.csv"), inst.nodes[cyclex],delim='\t')
            CSV.write(string(outdir,i,"routesy.csv"), inst.nodes[cycley],delim='\t')
            # saving in JSON
            open(string(outdir,i,"routes.json"),"w") do f
                JSON.print(f, inst.nodes[cyclex], 4)
            end
        end



end
