function plot_f82(f82_df)
########################################
# Plot current speed and direction from f82_df
    if !isempty(f82_df)
        p1_f82 = plot(f82_df.Date,f82_df.Speed,label="Speed (m/s)", ylabel="Speed (m/s)", yguidefontcolor=:red,c="red", lw=2, fillrange = 0, fillalpha = 0.025, fillcolor = :red, legend=:bottomleft,
            gridstyle=:dot, left_margin = 20Plots.mm, right_margin = 20Plots.mm, foreground_color_grid="red", yforeground_color_text="red", ytickfonthalign=:left)
        
        subplot = Plots.twinx()
        p1_f82 = plot!(subplot,f82_df.Date,f82_df.Direction,label="Direction (True)", ylabel="Direction (ᵒ True)", yguidefontcolor=:blue, c="blue", line=:dot, lw=2, legend=:bottomright,
            gridstyle=:dashdot, ylim=(0,360), yflip=true, yticks = 0:30:360, showaxis=:y, foreground_color_grid="blue", yforeground_color_text="blue", ytickfonthalign=:right)
        
        plot_f82 = Plots.plot(p1_f82, layout = (1, 1), size = (1500, 500),
                xlim=(first(f82_df.Date),last(f82_df.Date)),  xticks = first(f82_df.Date):Hour(12):last(f82_df.Date), xtickfontsize=7, ytickfontsize=8,
                framestyle = :box, fg_legend=:transparent, title = "Current speed and direction", margin=5Plots.mm, 
                grid=true, gridlinewidth=0.5, gridalpha=1)            

            display(plot_f82)
    end

end    # plot_f82()


function plot_f29(f29_df)
########################################
# Plot H3 Hmax TH3 THmax as calculated by Datawell and stored to F29

    if !isempty(f29_df)
    
        # display plots to screen
        tm_tick = range(first(f29_df.Date),last(f29_df.Date),step=Hour(2))
        ticks = Dates.format.(tm_tick,"HH\ndd")
        
        date_range = Dates.format(first(f29_df.Date), "dd/mm/yyyy HH:MM") * " - " * Dates.format(last(f29_df.Date), "dd/mm/yyyy HH:MM")

        p1_f29 = plot(f29_df.Date,f29_df.Hmax,label="Hmax", c="red", lw=2, fillrange = 0, fillalpha = 0.025, fillcolor = :red, title=date_range, titlefontsize=12)
        p1_f29 = plot!(f29_df.Date,f29_df.H3,label="Hsig", c="blue", lw=2, fillrange = 0, fillalpha = 0.05, fillcolor = :blue)

        p2_f29 = plot(f29_df.Date,f29_df.THmax,label="THmax", c="red", lw=2, fillrange = 0, fillalpha = 0.025, fillcolor = :red)
        p2_f29 = plot!(f29_df.Date,f29_df.TH3,label="THsig", c="blue", lw=2, fillrange = 0, fillalpha = 0.05, fillcolor = :blue)
        
        # display plots to screen
        plot_f29 = Plots.plot(p1_f29, p2_f29,layout = (2, 1), size = (1400, 700),
            xlim=(first(f29_df.Date),last(f29_df.Date)),  xticks=(tm_tick,ticks), xtickfontsize=7, ytickfontsize=8,
            framestyle = :box,fg_legend=:transparent, legend=:bottomleft,
            margin=1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

        display(plot_f29)
        
    else
    
        println("ALERT: No F29 data available!")

    end    
       
end    # plot_f29()


function plot_spectra(f23_df,f20_df,Data,found_list)
################################################
# function to plot spectra calculated by Datawell and from WSEs
    println("Preparing to plot spectra")

    # Extract parameters from F23 df
    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)
    heave, north, west = get_hnw(Data,start_val,end_val)

    # DWR4 parameters
    N = 512
    Sample_frequency = 2.56
    dt = 1/Sample_frequency

    # calc window
    w = DSP.Windows.hanning(N)
    w_norm = sqrt(Sample_frequency*sum(w.^2))

    start = 0
    segments = []

    # Do spectral analysis for 17 individual segments of 512 water surface elevations
    for segment in (1:17)
        finish = start+N

        # calculate the Fourier coefficients vide (5.6.2)
        Hl = do_fft(heave[start+1:finish].*w./w_norm, N)

        # calculate power spectral density
        Pden = calc_psd(Hl, N);
        segments = [segments;Pden]

        # set start for next segment
        start = start + 256
    end 

    # convert vector to matrix of 17 individual 256-value spectra
    combined_segments = reshape(segments,256,17);
    Pden = mean.(eachrow(combined_segments))

    freqs = [0:1:N/2-1;]*Sample_frequency/N

    ps_w = welch_pgram(heave, 512, 256; onesided=true, nfft=512, fs=Sample_frequency, window=hanning);
    f2 = freq(ps_w);
    Pden2 = power(ps_w);

    freqs = [1:1:N/2;]*Sample_frequency/N

    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")

    # Calculate the frequency bins for Datawell's datawell_spectra
    datawell_freqs = []

    for i in 0:1:45
        push!(datawell_freqs,0.025 + i*0.005)
    end

    for i in 46:1:79
       push!(datawell_freqs,-0.20 + i*0.01)
    end

    for i in 80:1:99
        push!(datawell_freqs,-0.98 + i*0.02)
    end

    normalized_datawell_spectra = Vector(f20_df[f20_df.Date .== start_date,4:103][1,1:100])
    smax = values(f20_df[f20_df.Date .== start_date, 3])

    datawell_spectra = normalized_datawell_spectra .* smax[1]

    # determing maximum y-axis value for calculated spectra
#    calc_max_y = maximum([maximum(Pden),maximum(Pden2)])

    # determing maximum y-axis value for Datawell spectra
    datawell_max_y = maximum([maximum(datawell_spectra),maximum(datawell_spectra)])

#    max_y = maximum([calc_max_y,datawell_max_y[1]]) * 1.05
    
    # get frequency-domain parameters calculated from the spectra
    Hm0, Hrms, T01, T02, Tc, Tp, Tp5, Skewness = calculate_frequency_domain_parameters(f2, Pden2)
    
    # Calculate representative spectra for P-M and JONSWAP
    Spectra_PM = calc_representative_spectra(f2, Hm0, Tp, 1.0);
    Spectra_JONSWAP = calc_representative_spectra(f2, Hm0, Tp, 3.3);
#    max_JONSWAP = maximum(Spectra_JONSWAP) * 1.05

    # Get the maximum value on Y-axis
    max_y = maximum([maximum(Pden),maximum(Pden2),datawell_max_y[1],maximum(Spectra_JONSWAP)]) * 1.05

    # Plot the representative spectra
    p_spectra = plot(f2,Spectra_JONSWAP, lw=2, c=:lightblue, label="JONSWAP spectrum (γ = 3.3)")
    p_spectra = plot!(f2,Spectra_PM, lw=2, c=:lightgreen, label="Pierson-Moskowitz spectrum (γ = 1.0)")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    # Add frequency-domain parameters to plot
    x_lim = xlims(p_spectra)[1]; y_lim = ylims(p_spectra)[2]
    p_spectra = annotate!(x_lim*-15, y_lim*0.95, "Hm0 = " * string(round(Hm0, digits=2)) * "m") 
    p_spectra = annotate!(x_lim*-15, y_lim*0.90, "Hrms = " * string(round(Hrms, digits=2)) * "m") 
    p_spectra = annotate!(x_lim*-15, y_lim*0.85, "T01 = " * string(round(T01, digits=2)) * "s") 
    p_spectra = annotate!(x_lim*-15, y_lim*0.80, "T02 = " * string(round(T02, digits=2)) * "s") 
    p_spectra = annotate!(x_lim*-15, y_lim*0.75, "Tp = " * string(round(Tp, digits=2)) * "s") 
    p_spectra = annotate!(x_lim*-15, y_lim*0.70, "Tp5 = " * string(round(Tp5, digits=2)) * "s") 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    # plot calculated spectra
    p_spectra = plot!(freqs, Pden, c="blue", lw=3, label="Calc. spectra\n", 
        fillrange = 0, fillalpha = 0.075, fillcolor = :blue)

    # plot Welch's spectra
    p_spectra = plot!(f2, Pden2, c="green", lw=3, label="Welch's spectra\n", 
        fillrange = 0, fillalpha = 0.075, fillcolor = :green)

    # plot Datawell spectra
    p_spectra = plot!(datawell_freqs,datawell_spectra, c="red", lw=3, label="Datawell spectra",
        fillrange = 0, fillalpha = 0.075, fillcolor = :red)

    p_spectra = vline!([1/Tp5; 1/Tp5], lw=:1, ls =:dash, c=:red, label="Tp5")

    plot_spc = Plots.plot(p_spectra, layout = (1, 1), size = (1400, 800),
            xlim=(0,1.0),  xticks = 0:0.05:1.0, xtickfontsize=7, xlabel="Frequency (Hertz)",
            ylim=(0,max_y), yticks=0:.5:max_y, ytickfontsize=8, ylabel="Spectral Density (m²/Hz.)",
            framestyle = :box, fg_legend=:transparent,
            title=title_string, titlefontsize=12,
            grid=true, gridlinewidth=0.5, gridalpha=1, foreground_color_grid="lightgrey", margin = 5Plots.mm)            

    display(plot_spc)

end    # plot_spectra()


function plot_hnw(f23_df,fc1_df,Data,found_list)
######################################## 

    function spike_value(wse)
    #####################################    
        median_value = median(wse)
        std_value = std(wse)
        
        return(median_value + 3*std_value)
        
        end    # spike_value()


    println("Preparing to plot heave, north, and west time series")
    # Extract parameters from F23 df
    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)

    # get WSEs for desired 30-minute record
    heave, north, west = get_hnw(Data,start_val,end_val)

    spike = spike_value(heave)
    heave_spikes = findall(i->(i>=spike), abs.(heave));

    spike = spike_value(north)
    north_spikes = findall(i->(i>=spike), abs.(north));

    spike = spike_value(west)
    west_spikes = findall(i->(i>=spike), abs.(west));

    # time stamp each WSE
    points = collect(0:1:length(heave)-1)/2.56
    times = []

    for i in 1:length(points)
        push!(times,unix2datetime(datetime2unix(start_date) + points[i]))
    end

    # create plots of heave, north, and west
    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")
    p1_hnw = Plots.scatter(times[heave_spikes], heave[heave_spikes], label="", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p1_hnw = plot!(times,heave, label="", c="#4a536b", lw=0.5, title=title_string, titlefontsize=12) ##last(split(infil,"\\")))

    # get plotting limits
    x_lim1 = xlims(p1_hnw)[1]; y_lim1 = ylims(p1_hnw)[1]
    x_lim2 = xlims(p1_hnw)[2]; y_lim2 = ylims(p1_hnw)[2]

    x_pos = x_lim1 + abs(x_lim2-x_lim1)*0.02
    p1_hnw = annotate!(x_pos, y_lim2*1.1, Plots.text("Firmwave ver. = " * fc1_df.Firmware[1], :grey, :left, 7))
    x_pos = x_lim1 + abs(x_lim2-x_lim1)*0.13
    p1_hnw = annotate!(x_pos, y_lim2*1.1, Plots.text("Hatch UID = " * string(fc1_df.Hatch_uid[1]), :grey, :left, 7))
    x_pos = x_lim1 + abs(x_lim2-x_lim1)*0.26
    p1_hnw = annotate!(x_pos, y_lim2*1.1, Plots.text("Hull UID = " * string(fc1_df.Hull_uid[1]), :grey, :left, 7))

    p2_hnw = Plots.scatter(times[north_spikes], north[north_spikes], label="", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p2_hnw = plot!(times,north, label="", c="#aed6dc", lw=0.5)
    p3_hnw = Plots.scatter(times[west_spikes], west[west_spikes], label="", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p3_hnw = plot!(times,west, label="", c="#ff9a8d", lw=0.5)

    hline!(p1_hnw, [0], lw=1, label="")
    hline!(p2_hnw, [0], lw=1, label="")
    hline!(p3_hnw, [0], lw=1, label="")

    # get plotting limits
    x_lim1 = xlims(p1_hnw)[1]; y_lim1 = ylims(p1_hnw)[1]
    x_lim2 = xlims(p1_hnw)[2]; y_lim2 = ylims(p1_hnw)[2]

    # display plots to screen
    plot_wse = Plots.plot(p1_hnw, p2_hnw, p3_hnw, layout = (3, 1), size = (1400, 900),
        xlim=(first(times),last(times)),  xticks = first(times):Minute(5):last(times),xtickfontsize=7,ytickfontsize=8,
        framestyle = :box,fg_legend=:transparent, legend=:bottomleft,
        margin = 1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

    display(plot_wse)

    # create a plot file to be saved as a .PNG
##    plt_file = first(infil, length(infil)-4)*"_plot_hnw_"*Dates.format(start_date, "yyyy_mm_dd_HHMM")*".png"

    # Save plot to file
##    savefig(plt_file)
##    println("Plot file saved as ",plt_file)
       
end    # plot_hnw()


function plot_3d(f23_df,Data,found_list)
######################################## 
 
    println("Preparing 3D plot")

    # Extract parameters from F23 df
    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)

    # get WSEs for desired 30-minute record
    heave, north, west = get_hnw(Data,start_val,end_val)

    # create plots of heave, north, and west
    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")
    p1 = plot(-west, north, heave, zcolor=heave, m=(1, 3, :RdYlGn_11, Plots.stroke(0)), leg=false, 
        cbar=true, c="lightgrey", w=0.5, xlabel="East-West", ylabel="North-South", zlabel="Heave")
    # display plots to screen
    plotty = Plots.plot(p1, layout = (1, 1), size = (1000, 1000), framestyle = :box, title=title_string);         

    display(plotty)
   
end    # plot_3d()


function plot_2d(f23_df,Data,found_list)
######################################## 
# function to do a 2-d plot of WSEs from West and North arrays
    println("Preparing 2d plot")
        
    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)
    
    heave, north, west = get_hnw(Data,start_val,end_val)

    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")

    east = -west

    p1 = plot(east,north, zcolor=heave, m=(1, 3, :RdYlGn_11, Plots.stroke(0)), leg=false, cbar=true, c="lightgrey", colorbar_title="Heave displacement (m)",
            w=0.5, xlabel="East-West displacement (m)", ylabel="North-South displacement (m)")

    plot_2d = Plots.plot(p1, layout = (1, 1), size = (900, 900),  xlim=(minimum(east)*1.1,maximum(east)*1.1), ylim=(minimum(north)*1.1,maximum(north)*1.1), 
        margin = 1Plots.mm, framestyle = :box, title=title_string, titlefontsize=12,
        grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1) 

    display(plot_2d)
    
end    # plot_2d()


function plot_hnw_2d(f23_df,Data,found_list)
########################################
# function to produce 2-d plots of displacements W-N; N-H; W-H

    println("Preparing Heave-North-West 2d plot")


    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)

    heave, north, west = get_hnw(Data,start_val,end_val)

    west = -west

    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")

    p_title = plot( title=title_string, grid = false, showaxis = false,top_margin=20Plots.mm)

    p1 = plot(west,north, size = (450, 450), zcolor=heave, m=(1, 3, :RdYlGn_11, Plots.stroke(0)), leg=false, cbar=false, c="lightgrey", colorbar_title="Heave displacement (m)",
        xlim=(minimum(west)*1.1,maximum(west)*1.1), ylim=(minimum(north)*1.1,maximum(north)*1.1), 
            w=0.5, title="West-North", titlefontsize=10, xlabel="East-West displacement (m)", ylabel="North-South displacement (m)")

    p2 = plot(north,heave, size = (450, 450), zcolor=heave, m=(1, 3, :RdYlGn_11, Plots.stroke(0)), leg=false, cbar=false, c="lightgrey", colorbar_title="Heave displacement (m)",
        xlim=(minimum(north)*1.1,maximum(north)*1.1), ylim=(minimum(heave)*1.1,maximum(heave)*1.1), 
            w=0.5, title="North-Heave", titlefontsize=10, xlabel="North-South displacement (m)", ylabel="Heave displacement (m)")

    p3 = plot(west,heave, size = (450, 450), zcolor=heave, m=(1, 3, :RdYlGn_11, Plots.stroke(0)), leg=false, cbar=false, c="lightgrey", colorbar_title="Heave displacement (m)",
        xlim=(minimum(west)*1.1,maximum(west)*1.1), ylim=(minimum(heave)*1.1,maximum(heave)*1.1), 
            w=0.5, title="West-Heave", titlefontsize=10, xlabel="East-West displacement (m)", ylabel="Heave displacement (m)")

    plots_2d = Plots.plot(p1, p2, p3, size = (1500, 500), layout = @layout([B C D]),plot_title=title_string, 
        ymirror = true, right_margin = 15Plots.mm, bottom_margin = 15Plots.mm, framestyle = :box, titlefontsize=12,
        grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1) 

    display(plots_2d)
    
end    # plot_hnw_2d()  

"""
function plot_wavelet(f23_df,Data,found_list)
######################################## 

    println("Preparing to plot wavelet - this plot takes a while!")
 
    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)
    end_date = start_date + Minute(30)
   
    heave, north, west = get_hnw(Data,start_val,end_val)
    # time stamp each WSE
    points = collect(0:1:length(heave)-1)/2.56
    times = []
    for i in 1:length(points)
        push!(times,unix2datetime(datetime2unix(start_date) + points[i]))
    end

    fs = 2.56

    @suppress begin    ## Suppress wavelet warning messages - see https://github.com/JuliaIO/Suppressor.jl
 
        n = length(heave)
        c = wavelet(Morlet(π), β=2);

        # plotting
        res = ContinuousWavelets.cwt(heave, c)
        freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
        freqs[1] = 0

##        display("text/html", "<style>.container { width:100% !important; }</style>")

        # get title string for plot
        title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")

        p1 = plot(times, heave,legend=false,title=title_string, ylabel= "Heave (m)", xticks=true, framestyle=:box)
        ##vline!(xmin:Dates.Hour(6):xmax,c="lightgrey",lw=0.5)


        p2 = contourf(times, freqs, abs.(res)', c=cgrad(:Spectral,rev=true), colorrange=(0, maximum(real(res))), levels=5,
            xtickfont = :serif,titlefontsize=13, ylim=(0,300), yticks=(0:25:300), ytickfont=:serif,
            xlabel= "Date (EST)", ylabel= "Period (s)", colorbar=false, framestyle = :box)
        p2 = hline!([25:25:275],c="white",lw=0.5,label="")
        p2 = vline!(start_date:Dates.Minute(6):end_date,c="white",lw=0.5,label="")

        l = @layout [a{.3h};b{.7h}]

        plot_heave_and_wavelet = plot(p1, p2, size=(1500, 800), xlim=(start_date,end_date),  
            xticks = first(times):Minute(6):last(times), xminorticks=12, tickdir=:out, xtickfontvalign = :bottom, ytickfonthalign = :left, 
            left_margin = 10Plots.mm, bottom_margin = 10Plots.mm,layout=l)

        display(plot_heave_and_wavelet)
        
    end    

end    # plot_wavelet()
"""

function plot_heave_histogram(f23_df,Data,found_list)

    println("Preparing to do Histogram plot")
       
    start_date, start_val, end_val = get_start_end_dates(f23_df,found_list)

    heave, north, west = get_hnw(Data,start_val,end_val)

    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM")

    # get the lower increment value
    floor_val = floor(minimum(heave)*10)/10;
    if isodd(floor_val*10)
        floor_val = (floor_val*10 - 1) / 10;
    end

    # get the upper increment value
    ceiling_val = ceil(maximum(heave)*10)/10;
    if isodd(ceiling_val*10)
        ceiling_val = (ceiling_val*10 + 1) / 10;
    end
    
    incr = 0.2    # bin width

#    println(floor_val,' ',ceiling_val)

    # Plot a histrgram of the heave with bins that are 0.2m wide
    p_hist = histogram(clamp.(heave, floor_val, ceiling_val), bins = round(Int,(ceiling_val-floor_val)/incr), xlabel = "WSE", leg = false, xlim=[floor_val,ceiling_val],
                title=title_string, titlefontsize=12, 
                bottom_margin = 10Plots.mm, framestyle = :box, size = (900, 700))
    
    # Determine plotting limits of p_hist so we can add annotations in desired position on the plot
    x_lim = xlims(p_hist)[1]; y_lim = ylims(p_hist)[2]
    p_hist = annotate!(x_lim*0.75, y_lim*0.90, "Skewness = " * string(round(skewness(heave), digits=4))) 
    p_hist = annotate!(x_lim*0.75, y_lim*0.85, "Kurtosis = " * string(round(kurtosis(heave), digits=4))) 

    # do a Normal Distribution plot
    L = Normal(mean(heave), std(heave))
    # do a Normal Distribution plot
    p_hist = plot!(Plots.twinx(), x->pdf(L, x), c = "red", lw = 3, framestyle = :box, size = (1400, 800), 
        xlim=[floor_val,ceiling_val], showaxis=:y, 
        label = "Normal distribution", fg_legend=:transparent)

    plot_hist = plot(p_hist, layout = (1, 1), link=:x, xticks = ([floor_val; floor_val:0.2:ceiling_val; ceiling_val], 
        ["< "*string(floor_val); floor_val:0.2:ceiling_val;
        "> "*string(ceiling_val)]), xrot = 90, ylim=[0,Inf], margin=15Plots.mm)

    display(plot_hist)
    
end    # plot_heave_histogram()


function plot_normal(heave)
    
    d = Normal(mean(heave), std(heave))
    lo, hi = minimum(heave),maximum(heave)
    x = range(lo, hi; length = length(heave))

    # get the lower increment value
    floor_val = floor(minimum(heave)*10)/10;
    if isodd(floor_val*10)
        floor_val = (floor_val*10 - 1) / 10;
    end

    # get the upper increment value
    ceiling_val = ceil(maximum(heave)*10)/10;
    if isodd(ceiling_val*10)
        ceiling_val = (ceiling_val*10 + 1) / 10;
    end

    incr = 0.2    # bin width

    @printf("%s\n","      BIN       Count")
    @printf("%s"," ======================")

    for i in floor_val:incr:ceiling_val

        bin_count = count(j->(j>i && j< i+incr), heave)

        if bin_count > 0
            @printf("\n%5.2f to %5.2f  %5d", i, i+incr, bin_count)
        end

    end

    normal_plot = histogram(clamp.(heave, floor_val, ceiling_val), 
        ylabel = "No. of WSEs per bin", bins = round(Int,(ceiling_val-floor_val)/incr), alpha=0.5)
    normal_plot = plot!(Plots.twinx(), x, (pdf.(d, x)), c=:red, lw=4, ticks=false)


    normal_plot = plot!(Plots.twinx(), x, cdf.(d, x) * length(heave), 
        c=:green, lw=4, ls=:dot, label = "Cumulative WSEs", ylabel = "Cumulative WSEs")


    show_normal_plot = plot(normal_plot, title=title_string, 
        framestyle = :box, size = (1400, 800), xlim=(floor_val,ceiling_val), legend=:topleft, 
        leftmargin = 20Plots.mm, rightmargin = 20Plots.mm, fg_legend=:transparent)

    display(show_normal_plot)
    
    return()
    
    end    # plot_normal()


function plot_rayleigh(heights)
    
    d = Normal(mean(heights), std(heights))
    lo, hi = minimum(heights),maximum(heights)
    x = range(lo, hi; length = length(heights))

    # get the lower increment value
    floor_val = floor(minimum(heights)*10)/10;
    if isodd(floor_val*10)
        floor_val = (floor_val*10 - 1) / 10;
    end

    # get the upper increment value
    ceiling_val = ceil(maximum(heights)*10)/10;
    if isodd(ceiling_val*10)
        ceiling_val = (ceiling_val*10 + 1) / 10;
    end

    incr = .2    # bin width

    # display number of waves in each wave-height bin over the range of waves
    @printf("%s\n","      BIN       Count")
    @printf("%s"," ======================")

    max_count = 0

    for i in floor_val:incr:ceiling_val
    
    ##    increasing = length(heights) - count(j->(j>i && j<= i+incr), heights)
        bin_count = count(j->(j>i && j< i+incr), heights)
        
        if bin_count > max_count
            max_count = bin_count
        end
        
        if bin_count>0
            @printf("\n%5.2f to %5.2f  %5d", i, i+incr, bin_count)
        end
    end

    # create some representative heights over range from 0 to max wave heights
    heights1 = lo:0.01:hi

    # get Scale of Rayleigh distribution
    σ = std(heights1)

    # calculate the representative Rayleigh over the range of wave heights
    rayleigh = []
    [push!(rayleigh, x/σ^2 * exp(-x^2/(2*σ^2))) for x in heights1];

    rayleigh_plot = histogram(heights, bins=:sqrt, xlabel = "Wave heights", ylabel = "No. of waves per bin", label="", alpha=0.25)
    rayleigh_plot = vline!([hmean; hmean], lw=1, ls =:dot, c=:red, label="")
    rayleigh_plot = vline!([hs; hs], lw=1, ls =:dot, c=:red, label="")
    rayleigh_plot = vline!([hmax; hmax], lw=1, ls =:dot, c=:red, label="")
    rayleigh_plot = annotate!(hmax, max_count*0.8, "Hmax", :red)
    rayleigh_plot = annotate!(hs, max_count*0.8, "Hs", :red)
    rayleigh_plot = annotate!(hmean, max_count*0.8, "Hmean", :red)

    rayleigh_plot = plot!(Plots.twinx(), heights1, rayleigh, lw=4, 
        c=:red, label = "Rayleigh distribution", ticks=false)

    rayleigh_plot = plot!(Plots.twinx(), x, cdf.(d, x) * length(heights), ylabel = "Cumulative waves", 
        c=:green, lw=4, ls=:dot, label = "Cumulative waves (total = "*string(length(heights))*")")
    
    do_rayleigh_plot = plot(rayleigh_plot, link=:x, title=title_string, titlefontsize=12, 
        xlim=(0,Inf),ylim=(0,Inf),
        bottom_margin = 10Plots.mm, leftmargin = 20Plots.mm, rightmargin = 20Plots.mm, 
        framestyle = :box, size = (1400, 800), fg_legend=:transparent)

    display(do_rayleigh_plot)
    
    return()
    
    end    # plot_rayleigh()
