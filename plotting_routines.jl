function simpleplot(Vec)
    """
    A single Polymer MWDw plot. (dots and smoothed line)
    Vec its a vector with the "Direct representation"
    """
    D66=removezeros(Vec)
    D666=((1:length(D66)).*D66;)./sum([1:length(D66)].*D66)

    plot(

    layer(x=1:length(D66),
    y=D666, Geom.smooth(method=:loess,smoothing=0.20),
    Theme(default_color=color("red"))),

    layer(x=1:length(D666), y=D666,
    Geom.point,
    Theme(default_point_size=0.6mm, default_color=color("black"))),

    Theme(grid_line_width=0mm,
    panel_stroke=color("black")),
    Guide.xticks(ticks=[0:100:length(D666)]),
    Guide.yticks(ticks=[0:0.001:maximum(D666)]),
    Guide.xlabel("chain length"),
    Guide.ylabel("weight franction", orientation=:vertical)

    )
end

function multipleplot(N,v1,v2...)
    """
    Plots multiple Polymer MWDWs with "Direct Representation".
    N is the array with the values of N to be plotted
    v1, v2, etc... are the arrays to plot.
    """
    DF =  DataFrame(x=1:length(getmwdw(v1)), y=getmwdw(v1)*1000, N=string(N[1]))
        for i = 2:(length(v2)+1)
            DF = vcat(DF,DataFrame(x=1:length(getmwdw(v2[i-1])), y=getmwdw(v2[i-1])*1000, N=string(N[i])))
        end

    p = plot(DF, x=:x, y=:y, color=:N, Geom.smooth(method=:loess,smoothing=0.30),

    Theme(grid_line_width=0mm,
    guide_title_position=:center,
    key_position=:right,
    panel_stroke=color("grey"),
    major_label_font_size=10pt),

    Scale.color_discrete_manual("red","orange","green","black","blue","cyan"),
    #Scale.y_continuous(minvalue=0, format=:scientific ),
    Guide.xticks(ticks=[0:100:500]),
    Guide.yticks(ticks=[0:1:7]),
    Guide.xlabel("chain length"),
    Guide.ylabel("weight fraction x 10⁻³", orientation=:vertical)
    )
end

function comparewithgproms(V1,G1,tipo="both")
    #V1 es el vector a comparar, previo calculo de mwdw
    #G1 es el vector resultados de gproms
    #tipo = string options are: dots, line, both

D11=removezeros(V1)
D111=((1:length(D11)).*D11;)./sum([1:length(D11)].*D11);

set_default_plot_size(20cm, 14cm)

if tipo == "both"

p=plot(

layer(x=1:length(D111),y=D111,
    Geom.smooth(method=:loess,smoothing=0.05),
    Theme(default_color=color("red"))),
    Guide.manual_color_key("",["Monte Carlo", "gProms"],["red", "blue"]),
    Guide.xticks(ticks=[0:100:length(D111)]),Guide.yticks(ticks=[0:0.001:maximum(D111)]),
    Guide.ylabel("Chain length"), # label for y-axis
    Guide.xlabel("Weight fraction"),  # label for x-axis
    Guide.title(""),

layer(x=1:length(D111), y=D111,
    Geom.point,
    Theme(default_point_size=0.6mm,default_color=color("red"))),


layer(x=1:length(G1),
    y=G1, Geom.line(),
    Theme(default_color=color("blue")))
)

elseif tipo == "dots"

p=plot(

layer(x=1:length(G1),
    y=G1, Geom.line(),
    Theme(default_color=color("blue"))),

layer(x=1:length(D111), y=D111,
    Geom.point,
    Theme(default_point_size=0.6mm,default_color=color("red"))),
    Guide.manual_color_key("",["Monte Carlo", "gProms"],["red", "blue"]),
    Guide.xticks(ticks=[0:100:length(D111)]),Guide.yticks(ticks=[0:0.001:maximum(D111)]),
    Guide.ylabel("Chain length"), # label for y-axis
    Guide.xlabel("Weight fraction"),  # label for x-axis
    Guide.title("")


)

elseif tipo == "line"

p=plot(

layer(x=1:length(D111),y=D111,
    Geom.smooth(method=:loess,smoothing=0.05),
    Theme(default_color=color("red"))),
    Guide.manual_color_key("",["Monte Carlo", "gProms"],["red", "blue"]),
    Guide.xticks(ticks=[0:500:length(D111)]),Guide.yticks(ticks=[0:0.001:maximum(D111)]),
    Guide.xlabel("Chain length"), # label for y-axis
    Guide.ylabel("Weight fraction"),  # label for x-axis
    Guide.title(""),

layer(x=1:length(G1),
    y=G1, Geom.line(),
    Theme(default_color=color("blue")))
)
end

return p

end

function plotwithribbon(RUNS)
    """
    This will plot the number of molecules for each chain length
    #RUNS is the matrix with the runs, each column a run.
    """

    #convierto todas las corridas al vector ordenado y sin ceros
    n=size(RUNS)[2]
    aux1=aux2=0
    for i = 1:n
        aux1=[removezeros(RUNS[:,i]) for i=1:n] #antes directtolinear
        aux2=[length(aux1[i])for i=1:n]
    end
    RUNSnz = zeros(Int16,maximum(aux2),n)
    for i = 1:n
        while length(aux1[i]) < maximum(aux2)
            push!(aux1[i],0)
        end
    end
    for i = 1:n
        RUNSnz[:,i]=aux1[i]
    end

    #segundo obtengo el promedio, maximo y minimo
    #squeeze() es para convertir la matrix nx1 en vector unidimensional
    D3mean = squeeze(mean(RUNSnz,2),2)
    D3max = float64(squeeze(maximum(RUNSnz,2),2))
    D3min = float64(squeeze(minimum(RUNSnz,2),2))

    #convierto estos valores a MWDw
    D3mean2 = ((1:length(D3mean)).*D3mean;)./sum(D3mean)
    D3max2 = ((1:length(D3max)).*D3max;)./sum(D3max);
    D3min2 = ((1:length(D3min)).*D3min;)./sum(D3min);

    #creo el dataframe para tener mas ordenado los datos y grafico
    df = DataFrames.DataFrame(x=1:length(D3mean2), y=D3mean, ymin=D3min, ymax=D3max)

    Gadfly.plot(df, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Gadfly.Geom.smooth(method=:loess,smoothing=0.3), Gadfly.Geom.ribbon)

end

function plotwithribbonMWDn(RUNS)
    """
    This will plot just the MWDn for each chain length
    #RUNS is the matrix with the runs, each column a run.
    """

    #convierto todas las corridas al vector ordenado y sin ceros
    n=size(RUNS)[2]
    aux1=aux2=0
    for i = 1:n
        aux1=[removezeros(RUNS[:,i]) for i=1:n] #antes directtolinear
        aux2=[length(aux1[i])for i=1:n]
    end
    RUNSnz = zeros(Int16,maximum(aux2),n)
    for i = 1:n
        while length(aux1[i]) < maximum(aux2)
            push!(aux1[i],0)
        end
    end
    for i = 1:n
        RUNSnz[:,i]=aux1[i]
    end

    #primero calculo la MWDn
    #segundo obtengo el promedio, maximo y minimo
    #squeeze() es para convertir la matrix nx1 en vector unidimensional
    RUNSnz2 = zeros(Float64,maximum(aux2),n)
    for i = 1:n
        RUNSnz2[:,i] = RUNSnz[:,i]/sum(RUNSnz[:,i])
    end

    D3mean = squeeze(mean(RUNSnz2,2),2)
    D3max = float64(squeeze(maximum(RUNSnz2,2),2))
    D3min = float64(squeeze(minimum(RUNSnz2,2),2))

    #creo el dataframe para tener mas ordenado los datos y grafico
    df = DataFrames.DataFrame(x=1:length(D3mean), y=D3mean, ymin=D3min, ymax=D3max)

    Gadfly.plot(df, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Gadfly.Geom.smooth(method=:loess,smoothing=0.3), Gadfly.Geom.ribbon)

end

function plotwithribbonMWDw(RUNS)
    """
    This will plot just the MWDw for each chain length
    #RUNS is the matrix with the runs, each column a run.
    """

    #convierto todas las corridas al vector ordenado y sin ceros
    n=size(RUNS)[2]
    aux1=aux2=0
    for i = 1:n
        aux1=[removezeros(RUNS[:,i]) for i=1:n] #antes directtolinear
        aux2=[length(aux1[i])for i=1:n]
    end
    RUNSnz = zeros(Int16,maximum(aux2),n)
    for i = 1:n
        while length(aux1[i]) < maximum(aux2)
            push!(aux1[i],0)
        end
    end
    for i = 1:n
        RUNSnz[:,i]=aux1[i]
    end

    #primer calculo la MWDw
    #segundo obtengo el promedio, maximo y minimo
    #squeeze() es para convertir la matrix nx1 en vector unidimensional

    RUNSnz2 = zeros(Float64,maximum(aux2),n)
    for i = 1:n
        RUNSnz2[:,i] = ((1:length(RUNSnz[:,i])).*RUNSnz[:,i])/sum([1:length(RUNSnz[:,i])].*RUNSnz[:,i])
    end

    D3mean = squeeze(mean(RUNSnz2,2),2)
    D3max = float64(squeeze(maximum(RUNSnz2,2),2))
    D3min = float64(squeeze(minimum(RUNSnz2,2),2))

    #creo el dataframe para tener mas ordenado los datos y grafico
    df = DataFrames.DataFrame(x=1:length(D3mean), y=D3mean, ymin=D3min, ymax=D3max)

    Gadfly.plot(df, x=:x, y=:y, ymin=:ymin, ymax=:ymax,
    Gadfly.Geom.smooth(method=:loess,smoothing=0.3),
    Gadfly.Geom.ribbon,
    Guide.xticks(ticks=[0:50:400]),
    Guide.yticks(ticks=[0:0.005:0.015]))

end

 """
    This will plot just the MWDw for each chain length with the ribbon smoothed.
    #RUNS is the matrix with the runs, each column a run.
"""

function getloess (x,y,smoothfactor)
        x=float64([1:length(y)]);
        y=squeeze(y,2);
        model=Loess.loess(x,y,span=smoothfactor);
        us = Loess.collect(minimum(x):0.1:maximum(x))
        vs = Loess.predict(model, us)
    return us,vs
end

function plotwithribbonMWDw_smoothed(RUNS) #RUNS es el vector directo de corridas, cada columna es una corrida

    #primero convierto todas las corridas al vector ordenado y sin ceros
    n=size(RUNS)[2]
    aux1=aux2=0
    for i = 1:n
        aux1=[removezeros(RUNS[:,i]) for i=1:n] #antes directtolinear
        aux2=[length(aux1[i])for i=1:n]
    end
    RUNSnz = zeros(Int16,maximum(aux2),n)
    for i = 1:n
        while length(aux1[i]) < maximum(aux2)
            push!(aux1[i],0)
        end
    end
    for i = 1:n
        RUNSnz[:,i]=aux1[i]
    end

    #segundo obtengo el promedio, maximo y minimo
    #squeeze() es para convertir la matrix nx1 en vector unidimensional
    RUNSnz2 = zeros(Float64,maximum(aux2),n)
    for i = 1:n
        RUNSnz2[:,i] = ((1:length(RUNSnz[:,i])).*RUNSnz[:,i])/sum([1:length(RUNSnz[:,i])].*RUNSnz[:,i])
    end

    D3mean = squeeze(mean(RUNSnz2,2),2)*1000
    D3max = float64(squeeze(maximum(RUNSnz2,2),2))*1000
    D3min = float64(squeeze(minimum(RUNSnz2,2),2))*1000

    upper = getloess(1:length(D3max),D3max,0.1)
    lower = getloess(1:length(D3min),D3min,0.1)

    #creo el dataframe para tener mas ordenado los datos y grafico
    df = DataFrames.DataFrame(x=1:length(D3mean), y=D3mean, ymin=D3min, ymax=D3max)

    Gadfly.plot(

    Gadfly.layer(df, x=:x, y=:y, Gadfly.Geom.smooth(method=:loess,smoothing=0.1),Theme(default_color=color("grey"))),
    Gadfly.layer(x=upper[1],y=df[:y],ymin=lower[2],ymax=upper[2],Gadfly.Geom.ribbon),

    Theme(grid_line_width=0mm,
    guide_title_position=:center,
    key_position=:none,
    panel_stroke=color("grey"),
    major_label_font_size=10pt),
Guide.xticks(ticks=[0:100:400]),
Guide.yticks(ticks=[0:1:7]),
Guide.xlabel("chain length"),
Guide.ylabel("weight fraction x 10⁻³", orientation=:vertical)

    )
end

function plotwithribbonMWDw(RUNS) #RUNS es el vector directo de corridas, cada columna es una corrida

    #primero convierto todas las corridas al vector ordenado y sin ceros
    n=size(RUNS)[2]
    aux1=aux2=0
    for i = 1:n
        aux1=[removezeros(RUNS[:,i]) for i=1:n] #antes directtolinear
        aux2=[length(aux1[i])for i=1:n]
    end
    RUNSnz = zeros(Int16,maximum(aux2),n)
    for i = 1:n
        while length(aux1[i]) < maximum(aux2)
            push!(aux1[i],0)
        end
    end
    for i = 1:n
        RUNSnz[:,i]=aux1[i]
    end

    #segundo obtengo el promedio, maximo y minimo
    #squeeze() es para convertir la matrix nx1 en vector unidimensional
    RUNSnz2 = zeros(Float64,maximum(aux2),n)
    for i = 1:n
        RUNSnz2[:,i] = ((1:length(RUNSnz[:,i])).*RUNSnz[:,i])/sum([1:length(RUNSnz[:,i])].*RUNSnz[:,i])
    end

    D3mean = squeeze(mean(RUNSnz2,2),2)*1000
    D3max = float64(squeeze(maximum(RUNSnz2,2),2))*1000
    D3min = float64(squeeze(minimum(RUNSnz2,2),2))*1000

    upper = getloess(1:length(D3max),D3max,0.1)
    lower = getloess(1:length(D3min),D3min,0.1)

    #creo el dataframe para tener mas ordenado los datos y grafico
    df = DataFrames.DataFrame(x=1:length(D3mean), y=D3mean, ymin=D3min, ymax=D3max)

    Gadfly.plot(

    Gadfly.layer(df, x=:x, y=:y, Gadfly.Geom.smooth(method=:loess,smoothing=0.1),Theme(default_color=color("grey"))),
    Gadfly.layer(x=1:500,y=df[:y],ymin=D3min,ymax=D3max,Gadfly.Geom.ribbon),

    Theme(grid_line_width=0mm,
    guide_title_position=:center,
    key_position=:none,
    panel_stroke=color("grey"),
    major_label_font_size=10pt),
Guide.xticks(ticks=[0:100:500]),
Guide.yticks(ticks=[0:1:7]),
Guide.xlabel("chain length"),
Guide.ylabel("weight fraction x 10⁻³", orientation=:vertical)

    )
end

function full_mwd_plot(DF)
    """
    """

    p = plot(DF, x=:x, y=:y, color=:N, Geom.smooth(method=:loess,smoothing=0.1),

    Theme(grid_line_width=0mm,
    guide_title_position=:center,
    key_position=:right,
    panel_stroke=color("grey"),
    major_label_font_size=10pt),

    Scale.color_discrete_manual("red","orange","green","black","blue","cyan"),
    #Scale.y_continuous(minvalue=0, format=:scientific ),
    Guide.xticks(ticks=[0:200:maximum(aaaa[:x])]),
    Guide.yticks(ticks=[0:0.001:maximum(DF[:y])]),
    Guide.xlabel("chain length"),
    Guide.ylabel("weight fraction x 10⁻³", orientation=:vertical)
    )
end
