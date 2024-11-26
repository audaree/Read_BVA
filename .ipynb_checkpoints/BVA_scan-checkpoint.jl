using Dates, DataFrames
using NativeFileDialog
using Plots, Printf

function plot_waves(Date, heave, north, west, heave_waves_valid_zero_up, north_waves_valid_zero_up, west_waves_valid_zero_up)
######################################## 

    println("Preparing to plot heave, north, and west time series")
    points_heave_waves = [ x[1] for x in heave_waves_valid_zero_up ];
    points_north_waves = [ x[1] for x in north_waves_valid_zero_up ];
    points_west_waves = [ x[1] for x in west_waves_valid_zero_up ];

    # create plots of heave, north, and west
    title_string = Dates.format(first(Date), "dd/mm/yyyy HH:MM") * " to " * Dates.format(last(Date), "dd/mm/yyyy HH:MM")
    p1_hnw = plot(Date[points_heave_waves], heave, label="Heave", c="#4a536b", lw=0.5, title=title_string, titlefontsize=12) ##last(split(infil,"\\")))
    p2_hnw = plot(Date[points_north_waves], north, label="North", c="#aed6dc", lw=0.5)
    p3_hnw = plot(Date[points_west_waves], west, label="West", c="#ff9a8d", lw=0.5)

    # display plots to screen
    plot_wse = Plots.plot(p1_hnw, p2_hnw, p3_hnw, layout = (3, 1), size = (1400, 900),
        xlim=(first(Date),last(Date)), ylim=(0,Inf), xtickfontsize=7,ytickfontsize=8,
        framestyle = :box,fg_legend=:transparent, legend=:topleft,
        margin = 1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

    display(plot_wse)

    # create a plot file to be saved as a .PNG
##    plt_file = first(infil, length(infil)-4)*"_plot_hnw_"*Date.format(start_date, "yyyy_mm_dd_HHMM")*".png"

    # Save plot to file
##    savefig(plt_file)
##    println("Plot file saved as ",plt_file)
       
    end    # plot_waves()


function get_displacement(Data, start_val, end_val)
################################################
# Decode the real time data to displacements - See DWTP (16 Jan 2019) 2.1.1 p. 19    
    
    arry = collect(Iterators.flatten(zip(SubString.(Data, start_val, end_val),SubString.(Data, start_val+9, end_val+9))));
    
    displacements = []
    
    for i in arry
        append!(displacements,parse(Int, SubString.(i, 1, 1), base=16)*16^2 + parse(Int, SubString.(i, 2, 2), base=16)*16^1 + parse(Int, SubString.(i, 3, 3), base=16)*16^0)
    end

    displacements[findall(>=(2048), displacements)] = displacements[findall(>=(2048), displacements)] .- 4096;
    displacements = 0.457*sinh.(displacements/457);    # see DWTP p.19 (16)
    
    return displacements
    
end    # get_displacement()


function get_hnw(Data,start_val,end_val)
######################################## 
    # get WSEs for desired 30-minute record
    heave = get_displacement(Data[start_val:end_val,:], 1, 3);              
    north = get_displacement(Data[start_val:end_val,:], 4, 6);
    west = get_displacement(Data[start_val:end_val,:], 7, 9);
    
    # Check for missing or extra points in data
    for wse in [heave, north, west]
        wse_length = length(wse)

        if wse_length > 4608
            
            wse = wse[1:4608]    # Only use the first 4608 points
        end

        if wse_length < 4608
            
            for i in wse_length:4608
                push!(wse,0)    # Zero pad WSEs out to 4608 points
            end
            
        end
    end
    
    return (heave, north, west)
    
    end    # get_hnw()


function get_heights_and_periods(wse)
########################################    
    zero_up = []; valid_zero_up = []; crest_points = []; trough_points = []

    for i in 2:length(wse)-1
        if (wse[i]*wse[i+1] < 0 && wse[i+1] > 0) || (wse[i] == 0 && wse[i-1] < 0 && wse[i+1] > 0)
            push!(zero_up,i)
        end
    end

    # wse Threshold set at 10mm. Refer to Section 9 Wave statistics pp. 9-10 in Datawell Library Manual
    threshold = 0.05

    i = 1; j = 2

    while j < length(zero_up)

        crest = maximum(wse[zero_up[i]:zero_up[j]])
        crest_point = zero_up[i] + argmax(wse[zero_up[i]:zero_up[j]]) - 1
        trough = minimum(wse[crest_point:zero_up[j]])

        # Check that crest higher than threshold AND trough less than threshold - Possible Valid Wave!!
        if (crest > threshold) & (trough < -threshold)
            crest_point = zero_up[i] + argmax(wse[zero_up[i]:zero_up[j]]) - 1
            trough_point = crest_point + argmin(wse[crest_point:zero_up[j]]) - 1

            push!(crest_points,crest_point)
            push!(trough_points,trough_point)

            next_crest = maximum(wse[zero_up[j]:zero_up[j+1]])

            # Check that NEXT crest also exceeds threshold (if so then Valid Wave)
            if (next_crest > threshold)
    ##            println("Crest found at ",crest_point," Trough at ",trough_point)
                push!(valid_zero_up,(zero_up[i],zero_up[j]));
                i = j
            end

        end

        j = j+1

    end

    # Process last recorded wave
    crest = maximum(wse[zero_up[i]:zero_up[j]])
    trough = minimum(wse[zero_up[i]:zero_up[j]])

    if (crest > threshold) & (trough < -threshold)
        
        crest_point = zero_up[i] + argmax(wse[zero_up[i]:zero_up[j]]) 
        trough_point = crest_point + argmin(wse[crest_point:zero_up[j]]) 
        push!(valid_zero_up,(zero_up[i],zero_up[j]));

    end

    heights = []

    for i in 1:length(valid_zero_up)

        crest = maximum(wse[valid_zero_up[i][1]:valid_zero_up[i][2]]);
        trough = minimum(wse[valid_zero_up[i][1]:valid_zero_up[i][2]]);
        push!(heights,crest - trough)
    ##    @printf("Wave %d = %2.3f\n",i,crest - trough)

    end 
    
    x_point = []

    for i in 1:length(valid_zero_up)
        push!(x_point,valid_zero_up[i][1] + abs(wse[valid_zero_up[i][1]]) / (wse[valid_zero_up[i][1]+1] - wse[valid_zero_up[i][1]]))
    end
    
    # need to get final zero crossing point - considering case where last valid zero crossing point is also last west value
    if last(valid_zero_up)[2] == length(wse)
        push!(x_point,last(valid_zero_up)[2] - wse[last(valid_zero_up)[2]])
    else
        push!(x_point,last(valid_zero_up)[2] + abs(wse[last(valid_zero_up)[2]]) / (wse[last(valid_zero_up)[2]+1] - wse[last(valid_zero_up)[2]]))
    end
    
    periods = x_point |> diff;
    
#    return(heights, periods, zero_up, valid_zero_up, x_point, crest_points, trough_points)
    return(heights, periods, x_point)
    
    end    # get_heights_and_periods()


function do_plots(df,parameter)
#########################################
## Display suspect values

    p1 = plot(df.Date,df.Height, lw=0.5, c=:blue, alpha=0.75, label=parameter*" heights")
    p1 = hline!([suspect_height],c=:red,label="Suspect height")

    p2 = plot(df.Date,df.Period, lw=0.5, c=:pink, label=parameter*" periods")
    p2 = hline!([suspect_period],c=:red,label="Suspect period")


    plot_heights = plot(p1, p2, layout = (2, 1), size = (1400, 300),framestyle = :box, xlim=(first(heave_params_df.Date),last(heave_params_df.Date)),
        fg_legend=:transparent, bg_legend=:transparent, legend=:topright,
        margin = 1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1, show=true)

    display(plot_heights)

    return()
    
    end


################################################
################################################
################################################
##           START OF MAIN PROGRAM
################################################
################################################
################################################

# Widen screen for better viewing
display("text/html", "<style>.container { width:100% !important; }</style>")

# Select a HVA or BVA file
infil = pick_file("C:\\QGHL\\Wave_data\\Bris\\BVA\\", filterlist="HVA,BVA;hva,bva");
println("Selected ",infil)

if uppercase(split(infil, ".")[end]) == "HVA"
    
    df = DataFrame(CSV.File(infil,header=0, delim=",-"));
    rename!(df,[:Sequence,:Data, :Packet]);

    # create a Packet vector of hex values
    packet = []
    
    for i in 1:length(df.Packet)
        push!(packet,string2hex(SubString(df.Packet[i],1,2)))
        push!(packet,string2hex(SubString(df.Packet[i],3,4)))
        push!(packet,string2hex(SubString(df.Packet[i],5,6)))
    end  
    
    # create a Data vector of hex values
    Data = df.Data
    
elseif uppercase(split(infil, ".")[end]) == "BVA"
    
    #Change the type-interpretation of the binary file data to unsigned integer
    println("Reading BINARY data from ",infil)
    data = reinterpret(UInt8, read(infil));

    # turn the data vector into a matrix of 12 values matching hexadecimal bytes - see DWTP 2.1 p.18
    cols = 12
    rows = Int(length(data) / cols)
    mat = reshape(view(data, :), cols, :);

    # Interleave last 4 matrix columns to form packet vector
    ## based on mschauer @ https://discourse.julialang.org/t/combining-two-arrays-with-alternating-elements/15498/2
    packet = collect(Iterators.flatten(zip(mat[10,:],mat[11,:],mat[12,:])));
    
    ## get data for the Heave, North, and West displacements
    Data = []

    # Convert binary data to hexidecimal vectors
    j = 0
    println("Building displacements vectors - this takes a while!")
    while true

        try
            heave_waves = []

            for i = j*12+1:j*12+12
                push!(heave_waves,string(data[i], base = 16, pad = 2))
            end

            push!(Data,join(heave_waves)[1:18])

        catch

            # escape if something is amiss        
            break

        end
        j = j+1

    end

else
    println("Not able to read this file type at present")
    exit()
end

start_val = 1
end_val = length(Data)

heave, north, west = get_hnw(Data,start_val,end_val);

start_date = split(infil, "\\")[end]
start_time = DateTime(start_date[1:4]*"-"*start_date[5:6]*"-"*start_date[7:8]*" 00:00:00", dateformat"y-m-d H:M:S")

Date = [start_time];

for i in 2:length(heave)
    push!(Date,start_time + Microsecond(390625*i))
end

# get the individual wave heights from Heave, North, and West displacements
heave_waves, heave_periods, heave_x_point = get_heights_and_periods(heave);
north_waves, north_periods, north_x_point = get_heights_and_periods(north);
west_waves, west_periods, west_x_point = get_heights_and_periods(west);

# get the individual wave heights from Heave, North, and West displacements
heave_waves, heave_periods, heave_x_point = get_heights_and_periods(heave);
north_waves, north_periods, north_x_point = get_heights_and_periods(north);
west_waves, west_periods, west_x_point = get_heights_and_periods(west);

# generate df's of heave, north, and west details (Date, Height, Period)
heave_params_df = DataFrame(Date = [], Height = [], Period = [])
north_params_df = DataFrame(Date = [], Height = [], Period = [])
west_params_df = DataFrame(Date = [], Height = [], Period = [])

for i in eachindex(heave_periods)
    push!(heave_params_df,[start_time + Microsecond.(ceil.((heave_x_point[i][1]/2.56) * 1000000)),heave_waves[i],heave_periods[i]/2.56]);
end

for i in eachindex(north_periods)
    push!(north_params_df,[start_time + Microsecond.(ceil.((north_x_point[i][1]/2.56) * 1000000)),north_waves[i],north_periods[i]/2.56])
end

for i in eachindex(west_periods)
    push!(west_params_df,[start_time + Microsecond.(ceil.((west_x_point[i][1]/2.56) * 1000000)),west_waves[i],west_periods[i]/2.56])
end

##############################################################################################################################################
# create plots of heave, north, and west
title_string = Dates.format(first(Date), "dd/mm/yyyy HH:MM") * " to " * Dates.format(last(Date), "dd/mm/yyyy HH:MM")
##p1_hnw = plot(heave_params_df.Date, heave_params_df.Height, label="Heave", c="#4a536b", lw=0.5, title=title_string, titlefontsize=12) ##last(split(infil,"\\")))
##p2_hnw = plot(north_params_df.Date, north_params_df.Height, label="North", c="#aed6dc", lw=0.5)
##p3_hnw = plot(west_params_df.Date, west_params_df.Height, label="West", c="#ff9a8d", lw=0.5)

# display plots to screen
##plot_wse = Plots.plot(p1_hnw, p2_hnw, p3_hnw, layout = (3, 1), size = (1400, 900),
##    xlim=(first(heave_params_df.Date),last(heave_params_df.Date)), ylim=(0,Inf), xtickfontsize=7,ytickfontsize=8,
##    framestyle = :box,fg_legend=:transparent, legend=:topleft,
##    margin = 1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

##display(plot_wse)
##############################################################################################################################################

# get the individual wave heights from Heave, North, and West displacements
heave_waves, heave_periods, heave_x_point = get_heights_and_periods(heave);
north_waves, north_periods, north_x_point = get_heights_and_periods(north);
west_waves, west_periods, west_x_point = get_heights_and_periods(west);

# get the individual wave heights from Heave, North, and West displacements
heave_waves, heave_periods, heave_x_point = get_heights_and_periods(heave);
north_waves, north_periods, north_x_point = get_heights_and_periods(north);
west_waves, west_periods, west_x_point = get_heights_and_periods(west);

# generate df's of heave, north, and west details (Date, Height, Period)
heave_params_df = DataFrame(Date = [], Height = [], Period = [])
north_params_df = DataFrame(Date = [], Height = [], Period = [])
west_params_df = DataFrame(Date = [], Height = [], Period = [])

for i in eachindex(heave_periods)
    push!(heave_params_df,[start_time + Microsecond.(ceil.((heave_x_point[i][1]/2.56) * 1000000)),heave_waves[i],heave_periods[i]/2.56]);
end

for i in eachindex(north_periods)
    push!(north_params_df,[start_time + Microsecond.(ceil.((north_x_point[i][1]/2.56) * 1000000)),north_waves[i],north_periods[i]/2.56])
end

for i in eachindex(west_periods)
    push!(west_params_df,[start_time + Microsecond.(ceil.((west_x_point[i][1]/2.56) * 1000000)),west_waves[i],west_periods[i]/2.56])
end

############################################
# set suspect values
suspect_height = 10
suspect_period = 15
############################################

heave_suspects = findall(( heave_params_df.Height .> suspect_height ) .|| (( heave_params_df.Height .> suspect_height ) .& ( heave_params_df.Period .> suspect_period)))
north_suspects = findall(( north_params_df.Height .> suspect_height ) .|| (( north_params_df.Height .> suspect_height ) .& ( north_params_df.Period .> suspect_period)))
west_suspects = findall(( west_params_df.Height .> suspect_height ) .|| (( west_params_df.Height .> suspect_height ) .& ( west_params_df.Period .> suspect_period)))

@printf("%s\n",title_string)
if isempty(heave_suspects)
    println("No suspect Heave values")
else
    do_plots(heave_params_df,"Heave")
    for i in heave_suspects
        @printf("Suspect Heave at %s %6.2fm %6.2fs\n",Dates.format(heave_params_df[i,:].Date, "yyyy-mm-dd HH:MM"), heave_params_df[i,:].Height, heave_params_df[i,:].Period)
    end
end

if isempty(north_suspects)
    
    println("No suspect North values")
else
    do_plots(north_params_df,"North")
    for i in north_suspects
        @printf("Suspect North at %s %6.2fm %6.2fs\n",Dates.format(north_params_df[i,:].Date, "yyyy-mm-dd HH:MM"), north_params_df[i,:].Height, north_params_df[i,:].Period)
    end
end

if isempty(west_suspects)
    println("No suspect West values")
else
    do_plots(west_params_df,"West")
    for i in west_suspects
        @printf("Suspect West at %s %6.2fm %6.2fs\n",Dates.format(west_params_df[i,:].Date, "yyyy-mm-dd HH:MM"), west_params_df[i,:].Height, west_params_df[i,:].Period)
    end
end

################################################
################################################
##           END OF MAIN PROGRAM
################################################
################################################