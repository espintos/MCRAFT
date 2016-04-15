function randrange(max,state)
    return a=(rem(state,max)+1)
end

function removezeros(Vec)
    return Vec[1:find(Vec)[end]]
end

function getmwdw(H)
    """
    Returns the MWD in weight of a distribution.
    H is a vector with the Histogram Representation
    """
    H=removezeros(H)
    MWD=((1:length(H)).*H)./sum([1:length(H)].*H)
    return MWD
end

function get_histogram(D)
    """
    Returns the Histogram Representation of a Distribution.
    D is a vector with the Direct Representation
    """
    if D[1]==0
        return zeros(1)
    else
        D=removezeros(D)
        max_length = maximum(D)
        H=zeros(Int,max_length)
        for i = 1:length(D)
          if D[i] == 0
          else
            H[D[i]] += 1
          end
        end
        return H
    end
end

function equalize_length2(names,v1...)
    max_length=1
    for i = 1:(length(v1))
        if length(v1[i]) > max_length
            max_length = length(v1[i])
        end
    end
    for i = 1:(length(v1))
        append!(v1[i],zeros(max_length-length(v1[i])))
    end

    suma = v1[1]
    for i = 2:(length(v1))
        suma += v1[i]
    end

    DF = DataFrame(x=1:max_length, y=v1[1], N=names[1])

    for i = 2:length(v1)
        DF = vcat(DF,DataFrame(x=1:max_length, y=v1[i], N=names[i]))
    end
    return suma,DF
end

function equalize_length3(names,v1...)
    max_length=1
    for i = 1:(length(v1))
        if length(v1[i]) > max_length
            max_length = length(v1[i])
        end
    end
    for i = 1:(length(v1))
        append!(v1[i],zeros(max_length-length(v1[i])))
    end

    suma = v1[1]
    for i = 2:(length(v1))
        suma += v1[i]
    end

    denom = sum([1:length(suma)].*suma)
    suma = getmwdw(suma)



    DF = DataFrame(x=1:max_length, y=suma, N="FULL")

    for i = 1:length(v1)
        DF = vcat(DF,DataFrame(x=1:max_length, y=([(1:length(v1[i]))].*v1[i])./denom, N=names[i]))
    end
    return DF
end

function Mn(Vec)
    aux=Vec #aux=directtolinear(Vec)
    Mnum=0.0
    for i = 1:length(aux)
        Mnum+=aux[i]*i*104.14
    end
    return Mnum/sum(aux)
end

function Mw(Vec)
    aux=Vec #aux=directtolinear(Vec)
    Mwe=0.0
    denom = 0.0
    for i = 1:length(aux)
        Mwe+=aux[i]*(i*104.14)^2
        denom+=aux[i]*i*104.14
    end
    return Mwe/denom
end
