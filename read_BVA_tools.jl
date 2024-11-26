
function get_window()
################################################
# calculate Hann window
    return [0.5*(1-cos(2*pi*k/N)) for k in (1:N)]
end    # get_window()


function do_fft(heave, N)
################################################
# calculate the Fourier coefficients vide (5.6.2)
    return([sum([heave[k]*exp(2*pi*-1im*k*l/N) for k in (1:N)]) for l in (1:N)])
end    # do_fft()


function calc_psd(Hl, N)
################################################
# The power spectral density is obtained from the Fourier coefficients    
    PSD = zeros(trunc(Int,N/2))

    for l = 1:trunc(Int,N/2)   
        if (l==1) || (l==trunc(Int,N/2)-1)
            PSD[l] = abs(Hl[l])^2
        else
            PSD[l] = abs(Hl[l])^2+abs(Hl[N-l-1])^2
        end
    end

    # Smooth coefficients vide (5.6.6)
    PSD_smooth = PSD
    [PSD_smooth[i] = PSD[i-1]/4 + PSD[i]/2 + PSD[i+1]/4 for i in (2:trunc(Int,N/2)-1)]

    return PSD_smooth    
end    # calc_psd()


function plot_spectra(heave)
################################################
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

    # determing maximum y-axis value for time-series plot
    max_y = maximum([abs(minimum(heave)),maximum(heave)]) * 1.05

    # determing maximum y-axis value for spectral plots
    max_y = maximum([maximum(Pden),maximum(Pden2)]) * 1.05

    # plot calculated spectra
    p_spectra = plot(freqs, Pden, label="Calc", 
        c="blue", lw=3, 
        size = (1200, 500), framestyle = :box, 
        xlim=(0,0.6), xticks=0:.05:1.2, xlabel="Frequency (Hertz)",
        ylim=(0,max_y), yticks=0:.5:max_y, ylabel="Spectral Density (sq.m/Hertz)",
        fillrange = 0, fillalpha = 0.075, fillcolor = :blue)
        
    plot_spc = Plots.plot(p_spectra, layout = (1, 1), size = (1400, 800),
            xlim=(0,1.28),  xticks = 0:0.05:1.28, xtickfontsize=7, ytickfontsize=8,
            framestyle = :box, fg_legend=:transparent, title = " Spectral plot", titlefontsize=12,
            grid=true, gridlinewidth=0.5, gridalpha=1, foreground_color_grid="lightgrey", margin = 5Plots.mm)            

    display(plot_spc)
    
end


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


function process_f20(f20_vals,heave_spectrum)
#######################################
# function to calculate parameters from Upcross wave height quantiles message  (0xf20)
# refer to DWTP (Ver. 16 January2019) Section 4.8 pp.56-61
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f20_vals[3])*bitstring(f20_vals[4])*bitstring(f20_vals[5])*bitstring(f20_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Number of Segments Used
    segments = parse(Int, bitstring(f20_vals[9]); base=2); 
    
    # get Smax
    smax = 5000 * (exp(parse(Int, bitstring(f20_vals[10]) * bitstring(f20_vals[11])[1:4]; base=2) / 200) - 1) / (exp(4094/200) - 1);    # (45) p.37
    # Obtain the heave_spectrum from s0 to s99
    for ii in 12:3:159
        push!(heave_spectrum,(exp(parse(Int, bitstring(f20_vals[ii]) * bitstring(f20_vals[ii+1])[1:4]; base=2) / 200) - 1) / (exp(4094/200) - 1))    # (48) p.38
        push!(heave_spectrum,(exp(parse(Int, bitstring(f20_vals[ii+1])[5:8] * bitstring(f20_vals[ii+2]); base=2) / 200) - 1) / (exp(4094/200) - 1))
    end

    return(timestamp,segments,smax,heave_spectrum)

    end    # process_f20()


function process_f21(f21_vals,direction,spread)
#######################################
# function to calculate parameters from Upcross wave height quantiles message  (0xf21)
# refer to DWTP (Ver. 16 January2019) Section 4.8 pp.56-61
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f21_vals[3])*bitstring(f21_vals[4])*bitstring(f21_vals[5])*bitstring(f21_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Number of Segments Used
    segments = parse(Int, bitstring(f21_vals[9]); base=2); 
    
    # Obtain the directions and spreads from s0 to s99
    for ii in 10:3:307
        push!(direction,rad2deg(parse(Int, bitstring(f21_vals[ii]) * bitstring(f21_vals[ii+1])[1:4]; base=2) / 4095 * 2 * pi))
        push!(spread,rad2deg(parse(Int, bitstring(f21_vals[ii+1])[5:8] * bitstring(f21_vals[ii+2]); base=2) / 4095 * pi / 2))
    end

    return(timestamp,segments,direction,spread)

    end    # process_f21()


function process_f23(f23_vals)
#######################################
# function to calculate selected parameters from Spectrum synchronisation message (0xF23)
# refer to DWTP (Ver. 16 January2019) Section 4.3 pp.43-44
# Called by: decode_packet()   
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f23_vals[3])*bitstring(f23_vals[4])*bitstring(f23_vals[5])*bitstring(f23_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Data Stamp
    data_stamp = parse(Int, bitstring(f23_vals[7])*bitstring(f23_vals[8]); base=2);

    # get Segments Used
    segments_used = parse(Int, bitstring(f23_vals[9])*bitstring(f23_vals[10])*bitstring(f23_vals[11]); base=2);

    # get Sample Number
    sample_number = parse(Int, bitstring(f23_vals[12])*bitstring(f23_vals[13]); base=2);

    match_vector = lpad(string(f23_vals[14], base = 16),2,"0")
    for i in 15:22
        match_vector = match_vector * lpad(string(f23_vals[i], base = 16),2,"0") 
    end   
##    println(timestamp,"UTC ",segments_used,' ',match_vector,' ',sample_number)
    
    return(timestamp,segments_used,match_vector,sample_number)
    
    end    # process_f23()


function process_f25(f25_vals)
#######################################
# function to calculate parameters from Directional spectral parameters message (0xf25)
# refer to DWTP (Ver. 16 January2019) Section 4.5 pp.49-51

    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f25_vals[3])*bitstring(f25_vals[4])*bitstring(f25_vals[5])*bitstring(f25_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Number of Segments used
    segments = parse(Int, bitstring(f25_vals[9]))
        
    # get Hs
    hs = parse(Int, bitstring(f25_vals[10]) * bitstring(f25_vals[11])[1:4]; base=2) / 100;    # (58)
    
    # get Integral Period TI
    ti = parse(Int, bitstring(f25_vals[11])[5:8] * bitstring(f25_vals[12]); base=2) / 100;    # (59)
    
    # get Energy Period TE
    te = parse(Int, bitstring(f25_vals[13]) * bitstring(f25_vals[14])[1:4]; base=2) / 100;    # (60)
    
    # get Mean Period T1
    t1 = parse(Int, bitstring(f25_vals[14])[5:8] * bitstring(f25_vals[15]); base=2) / 100;    # (61)

    # get Average Wave Period Tz
    tz = parse(Int, bitstring(f25_vals[16]) * bitstring(f25_vals[17])[1:4]; base=2) / 100;    # (62)
    
    # get T3
    t3 = parse(Int, bitstring(f25_vals[17])[5:8] * bitstring(f25_vals[18]); base=2) / 100;    # (63)

    # get Tc
    tc = parse(Int, bitstring(f25_vals[19]) * bitstring(f25_vals[20])[1:4]; base=2) / 100;    # (64)
    
    # get Rp
    rp = parse(Int, bitstring(f25_vals[20])[5:8] * bitstring(f25_vals[21]); base=2) / 4094;    # (65)

    # get Tp
    tp = parse(Int, bitstring(f25_vals[22]) * bitstring(f25_vals[23])[1:4]; base=2) / 100;    # (66)
    
    # get Smax
    smax = 5000 * (exp(parse(Int, bitstring(f25_vals[23])[5:8] * bitstring(f25_vals[24]); base=2) / 200) - 1) / (exp(4094/200) - 1);

    # get Theta_p
    theta_p = rad2deg(parse(Int, bitstring(f25_vals[25]) * bitstring(f25_vals[26])[1:4]; base=2) / 4095 * 2 * pi);    # (68)
   
    # get Sigma_p
    sigma_p = rad2deg(parse(Int, bitstring(f25_vals[26])[5:8] * bitstring(f25_vals[27]); base=2) / 4095 * pi / 2);    # (69)

    return(timestamp,segments,hs,ti,te,t1,tz,t3,tc,rp,tp,smax,theta_p,sigma_p)
    
    end    # process_f25()


function process_f26(f26_vals)
#######################################
# function to calculate parameters from Online upcross wave statistics message  (0xf26)
# refer to DWTP (Ver. 16 January2019) Section 4.6 pp.51-56

    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f26_vals[3])*bitstring(f26_vals[4])*bitstring(f26_vals[5])*bitstring(f26_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Hmax
    hmax = parse(Int, bitstring(f26_vals[9]) * bitstring(f26_vals[10])[1:4]; base=2) / 100;    # (70)
    
    # get THmax
    thmax = parse(Int, bitstring(f26_vals[10])[5:8] * bitstring(f26_vals[11]); base=2) / 100;    # (71)
    
    # get Tmax
    tmax = parse(Int, bitstring(f26_vals[12]) * bitstring(f26_vals[13])[1:4]; base=2) / 100;    # (72)
    
    # get HTmax
    htmax = parse(Int, bitstring(f26_vals[13])[5:8] * bitstring(f26_vals[14]); base=2) / 100;    # (73)

    # get Havg
    havg = parse(Int, bitstring(f26_vals[15]) * bitstring(f26_vals[16])[1:4]; base=2) / 100;    (74)
    
    # get Tavg
    tavg = parse(Int, bitstring(f26_vals[16])[5:8] * bitstring(f26_vals[17]); base=2) / 100;    (75)

    # get HsRMS
    hsrms = parse(Int, bitstring(f26_vals[18]) * bitstring(f26_vals[19])[1:4]; base=2) / 100;    # (76)
    
    # get Nw
    nw = parse(Int, bitstring(f26_vals[19])[5:8] * bitstring(f26_vals[20]); base=2);    # (77)

    # get Nc
    nc = parse(Int, bitstring(f26_vals[21]) * bitstring(f26_vals[22])[1:4]; base=2);    # (78)
    
    # get Epsilon
    epsilon = parse(Int, bitstring(f26_vals[22])[5:8] * bitstring(f26_vals[23]); base=2) / 4094;    # (79)

    # get Coverage
    coverage = parse(Int, bitstring(f26_vals[24]) * bitstring(f26_vals[25])[1:4]; base=2) / 4094 * 100;    # 80
   
    return(timestamp,hmax,thmax,tmax,htmax,havg,tavg,hsrms,nw,nc,epsilon,coverage)
    
    end    # process_f26()


function process_f28(f28_vals,m2,n2,k)
#######################################
# function to calculate parameters from Secondary directional spectrum message  (0xf28)
# refer to DWTP (Ver. 16 January2019) Section 4.2.4 pp.41-43
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f28_vals[3])*bitstring(f28_vals[4])*bitstring(f28_vals[5])*bitstring(f28_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Number of Segments Used
    segments = parse(Int, bitstring(f28_vals[9]); base=2); 
    
    # Obtain the directions and spreads from s0 to s99
    for ii in 10:9:451
        push!(m2,parse(Int, bitstring(f28_vals[ii]) * bitstring(f28_vals[ii+1])[1:4]; base=2) / 2047)
        push!(n2,parse(Int, bitstring(f28_vals[ii+1])[5:8] * bitstring(f28_vals[ii+2]); base=2) / 2047)
        push!(k,(exp(parse(Int, bitstring(f28_vals[ii+3]) * bitstring(f28_vals[ii+4])[1:4]; base=2) / 2124.5841)) - 1) / (exp(4094/2124.5841 - 1))    # (56) p.43
        push!(m2,parse(Int, bitstring(f28_vals[ii+4])[5:8] * bitstring(f28_vals[ii+5]); base=2) / 2047)
        push!(n2,parse(Int, bitstring(f28_vals[ii+6]) * bitstring(f28_vals[ii+7])[1:4]; base=2) / 2047)
        push!(k,(exp(parse(Int, bitstring(f28_vals[ii+7])[5:8] * bitstring(f28_vals[ii+8]); base=2) / 2124.5841)) - 1) / (exp(4094/2124.5841 - 1))
    end
    return(timestamp,segments,m2,n2,k)

    end    # process_f28()


function process_f29(f29_vals,hq)
#######################################
# function to calculate parameters from Upcross wave height quantiles message  (0xf29)
# refer to DWTP (Ver. 16 January2019) Section 4.8 pp.56-61
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f29_vals[3])*bitstring(f29_vals[4])*bitstring(f29_vals[5])*bitstring(f29_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Coverage
    coverage = parse(Int, bitstring(f29_vals[9]) * bitstring(f29_vals[10])[1:4]; base=2) / 4094 * 100;    # (81)
    
    # get Nw
    nw = parse(Int, bitstring(f29_vals[10])[5:8] * bitstring(f29_vals[11]); base=2);    # (82)
    
    # get Energy Period Epsilon
    epsilon = parse(Int, bitstring(f29_vals[12]) * bitstring(f29_vals[13])[1:4]; base=2) / 4094;    # (83)
    
    # get Hmax
    hmax = parse(Int, bitstring(f29_vals[13])[5:8] * bitstring(f29_vals[14]); base=2) / 100;    # (84)

    # get THmax
    thmax = parse(Int, bitstring(f29_vals[15]) * bitstring(f29_vals[16])[1:4]; base=2) / 100;    # (85)
    
    # get H10
    h10 = parse(Int, bitstring(f29_vals[16])[5:8] * bitstring(f29_vals[17]); base=2) / 100;    # (86)

    # get TH10
    th10 = parse(Int, bitstring(f29_vals[18]) * bitstring(f29_vals[19])[1:4]; base=2) / 100;    # (87)
    
    # get H3
    h3 = parse(Int, bitstring(f29_vals[19])[5:8] * bitstring(f29_vals[20]); base=2) / 100;    # (88)

    # get TH3
    th3 = parse(Int, bitstring(f29_vals[21]) * bitstring(f29_vals[22])[1:4]; base=2) / 100;    # (89)
    
    # get Havg
    havg = parse(Int, bitstring(f29_vals[22])[5:8] * bitstring(f29_vals[23]); base=2) / 100;    # (90)

    # get Tavg
    tavg = parse(Int, bitstring(f29_vals[24]) * bitstring(f29_vals[25])[1:4]; base=2) / 100;    # (91)
    
    # Obtain the % wave heights Hq0 to Hq21
    for ii in 25:3:55
        push!(hq,parse(Int, bitstring(f29_vals[ii])[5:8] * bitstring(f29_vals[ii+1]); base=2) / 100;)    # (92)
        push!(hq,parse(Int, bitstring(f29_vals[ii+2]) * bitstring(f29_vals[ii+3])[1:4]; base=2) / 100;)
    end

    # process the final % wave height Hq22
    push!(hq,parse(Int, bitstring(f29_vals[58])[5:8] * bitstring(f29_vals[59]); base=2) / 100;)
    
    return(timestamp,coverage,nw,epsilon,hmax,thmax,h10,th10,h3,th3,havg,tavg,hq)

    end    # process_f29()


function process_f80(f80_vals)
#######################################
# function to calculate parameters from Battery life expectancy message (0xf80)
# refer to DWTP (Ver. 16 January2019) Section 4.16 pp.91-92
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f80_vals[3])*bitstring(f80_vals[4])*bitstring(f80_vals[5])*bitstring(f80_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Latitude
    latitude = -(180 - parse(Int, bitstring(f80_vals[9]) * bitstring(f80_vals[10]) * bitstring(f80_vals[11]); base=2) / (2^24-1) * 180);
    
    # get Longitude
    longitude = parse(Int, bitstring(f80_vals[12]) * bitstring(f80_vals[13]) * bitstring(f80_vals[14]); base=2) / (2^24-1) * 360;

##    println(timestamp,"UTC ",' ',latitude,' ',longitude)
    
    return(timestamp,latitude,longitude)
    
    end    # process_f80()


function process_f81(f81_vals)
#######################################
# function to calculate selected parameters from SST message (0xf81)
# refer to DWTP (Ver. 16 January2019) Section 4.11 pp.67-68

    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f81_vals[3])*bitstring(f81_vals[4])*bitstring(f81_vals[5])*bitstring(f81_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get SST
    sst = parse(Int, bitstring(f81_vals[9])*bitstring(f81_vals[10]); base=2)/100 − 273.15;

##    println(timestamp,"UTC ",' ',sst)
    
    return(timestamp,sst)
    
    end    # process_f81()


function process_f82(f82_vals)
#######################################
# function to calculate selected parameters from Acoustic current meter message (0xf82)
# refer to DWTP (Ver. 16 January2019) Section 4.12 pp.68-71

    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f82_vals[3])*bitstring(f82_vals[4])*bitstring(f82_vals[5])*bitstring(f82_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get ACM Firmware Version
    firmware_version = parse(Int, bitstring(f82_vals[9])*bitstring(f82_vals[10])*bitstring(f82_vals[11])
        *bitstring(f82_vals[12])*bitstring(f82_vals[13])*bitstring(f82_vals[14])*bitstring(f82_vals[15])
        *bitstring(f82_vals[16]); base=2);

    # get Current Speed
    speed = parse(Int, bitstring(f82_vals[17])*bitstring(f82_vals[18])[1:4]; base=2)/1000;

    # get Current Direction
    direction = rad2deg(parse(Int, bitstring(f82_vals[18])[5:8]*bitstring(f82_vals[19]); base=2)/4095*2*pi);

    # get SST
    sst = parse(Int, bitstring(f82_vals[25])*bitstring(f82_vals[26]); base=2)/100 − 273.15;

##    println(timestamp,"UTC ",firmware_version,' ',speed,' ',direction,' ',sst)
    
    return(timestamp,firmware_version,speed,direction,sst)
    
    end    # process_f82()


function process_fc1(fc1_vals)
#######################################
# function to calculate parameters from System message for the DWR4   (0xfC1)
# refer to DWTP (Ver. 16 January2019) Section 4.15.2 pp.85-91
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(fc1_vals[3])*bitstring(fc1_vals[4])*bitstring(fc1_vals[5])*bitstring(fc1_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get Firmware Version Number
    mybytes = UInt8.(fc1_vals[9:18]);
    firmware = strip(String(mybytes),['\0'])
    
    # get Hatch UID
    hatch_uid = parse(Int, bitstring(fc1_vals[24]) * bitstring(fc1_vals[23]) * bitstring(fc1_vals[22])
        * bitstring(fc1_vals[21]) * bitstring(fc1_vals[20]) * bitstring(fc1_vals[19]); base=2) 
    
    # get Hull UID
    hull_uid = parse(Int, bitstring(fc1_vals[30]) * bitstring(fc1_vals[29]) * bitstring(fc1_vals[28])
        * bitstring(fc1_vals[27]) * bitstring(fc1_vals[26]) * bitstring(fc1_vals[25]); base=2) 
    
    # get Uptime
    uptime = parse(Int, bitstring(fc1_vals[34]) * bitstring(fc1_vals[33]) * bitstring(fc1_vals[32]) * bitstring(fc1_vals[31]); base=2) / 3600 /24  # time in days
    
    # get Energy used from batteries
    battery_energy = parse(Int, bitstring(fc1_vals[38]) * bitstring(fc1_vals[37]) * bitstring(fc1_vals[36]) * bitstring(fc1_vals[35]); base=2) 

    # get Energy to boostcaps
    boostcaps_energy = parse(Int, bitstring(fc1_vals[42]) * bitstring(fc1_vals[41]) * bitstring(fc1_vals[40]) * bitstring(fc1_vals[39]); base=2) 

    # get Hatch electronics temperature
    hatch_temp = (200 + parse(Int, bitstring(fc1_vals[42]); base=2)) - 273.15    # (121) p.79
    
    # get Battery Voltage
    battery_voltage = 6 + parse(Int, bitstring(fc1_vals[44]); base=2)/10    # (122) p.79

    # get Batteries per section
    batteries_per_section = parse(Int, bitstring(fc1_vals[45]); base=2) 
    
    # Number of battery sections
    battery_section_number = parse(Int, bitstring(fc1_vals[46]); base=2)

    # Initial energy in a battery
    initial_battery_energy = 36000 * parse(Int, bitstring(fc1_vals[47]); base=2)    # (126) p.80
    
    # get vertical accelerometer offset Ov
    ov = parse(Int, bitstring(fc1_vals[48]) * bitstring(fc1_vals[49])[1:4]; base=2) / 800    # (136) p.88

    # get the number of times the output of the vertical accelerometer reached its maximum value Cv
    cv = parse(Int, bitstring(fc1_vals[50]); base=2)
    
    # get X-axis accelerometer offset Ox
    ox = parse(Int, bitstring(fc1_vals[51]) * bitstring(fc1_vals[52])[1:4]; base=2) / 800    # (138) p.88

    # get y-axis accelerometer offset Oy
    oy = parse(Int, bitstring(fc1_vals[52])[5:8] * bitstring(fc1_vals[53]); base=2) / 800    # (139) p.88 

    # get the number of times the output of the X-axis accelerometer reached its maximum value Cx
    cx = parse(Int, bitstring(fc1_vals[54]); base=2)
    
    # get the number of times the output of the Y-axis accelerometer reached its maximum value Cy
    cy = parse(Int, bitstring(fc1_vals[54]); base=2)
    
    # get the average orientation of the buoy during last 30 minutes Mu0
    mu0 = rad2deg(parse(Int, bitstring(fc1_vals[56]) * bitstring(fc1_vals[57])[1:4]; base=2) / 4095 * 2 * pi)    # (142) p.89
 
    # get the standard deviation of the orientation of the buoy during last 30 minutes Sigma0
    sigma0 = rad2deg(parse(Int, bitstring(fc1_vals[57])[5:8] * bitstring(fc1_vals[58]); base=2) /4095 * pi / 2)    # (143) p.89
 
     # get the average inclination of the earth magnetic field of the buoy during last 30 minutes Mui
    mui = rad2deg(parse(Int, bitstring(fc1_vals[59]) * bitstring(fc1_vals[60])[1:4]; base=2) / 4095 * pi)    # (144) p.90
 
    # get the standard deviation of the earth magnetic field of the orientation of the buoy during last 30 minutes Sigmai
    sigmai = rad2deg(parse(Int, bitstring(fc1_vals[60])[5:8] * bitstring(fc1_vals[61]); base=2) / 4095 * pi / 16)    # (145) p.90
 
    # get the average length of the earth magnetic field vector during last 30 minutes MuH
    muh = parse(Int, bitstring(fc1_vals[62]) * bitstring(fc1_vals[63])[1:4]; base=2) / 4095 * 128 * 10^-6    # (146) p.90
 
    # get the standard deviation of the length of the earth magnetic field vector during last 30 minutes SigmaH
    sigmah = parse(Int, bitstring(fc1_vals[63])[5:8] * bitstring(fc1_vals[64]); base=2) / 4095 * 256 * 10^-9    # (146) p.90 
 
    # get the number of times the pitch angle reached its maximum value during last 30 minutes Cpitch
    cpitch = parse(Int, bitstring(fc1_vals[65]); base=2)
    
    # get the number of times the roll angle reached its maximum value during last 30 minutes Croll
    croll = parse(Int, bitstring(fc1_vals[66]); base=2)
    
    # get the temperature of the acceleration sensor Tsensor
    tensor = 200 + parse(Int, bitstring(fc1_vals[67]); base=2) - 273.15    # (150) p.91
    
    return(timestamp,firmware,hatch_uid,hull_uid,uptime,battery_energy,boostcaps_energy,hatch_temp,battery_voltage,batteries_per_section,
        battery_section_number,initial_battery_energy,ov,cv,ox,oy,cx,cy,mu0,sigma0,mui,sigmai,muh,sigmah,cpitch,croll,tensor)

    end    # process_fc1()


function process_fc3(fc3_vals)
#######################################
# function to calculate parameters from Battery life expectancy message (0xfc3)
# refer to DWTP (Ver. 16 January2019) Section 4.16 pp.91-92

    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(fc3_vals[3])*bitstring(fc3_vals[4])*bitstring(fc3_vals[5])*bitstring(fc3_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(10)

    # get ble
    ble = parse(Int, bitstring(fc3_vals[9]); base=2);

##    println(timestamp,"UTC ",' ',ble)
    
    return(timestamp,ble)
    
    end    # process_fc3()


function plot_hnw(f23_df,Data,found_list)
######################################## 
    
    # Extract parameters from F23 df
    start_date = f23_df[found_list[1],:].Date - Minute(30) # <------- NOTE subtracted 30min from start_date to match Waves4 results
    segments = f23_df[found_list[1],:].Segments
    match_vector = f23_df[found_list[1],:].Match_vector
    sample_nos = f23_df[found_list[1],:].Sample_number
    data_vector = f23_df[found_list[1],:].Data_vector
    start_val = data_vector - Int(sample_nos/2) + 1
    end_val = data_vector
##    println(start_date,' ',segments,' ',match_vector,' ',sample_nos,' ',start_val,' ',end_val,' ',Data[data_vector])

    # get WSEs for desired 30-minute record
    heave = get_displacement(Data[start_val:end_val,:], 1, 3);              
    north = get_displacement(Data[start_val:end_val,:], 4, 6);
    west = get_displacement(Data[start_val:end_val,:], 7, 9);

    # time stamp each WSE
    points = collect(0:1:sample_nos-1)/2.56
    times = []
    for i in 1:length(points)
        push!(times,unix2datetime(datetime2unix(start_date) + points[i]))
    end
    
    # create plots of heave, north, and west
    title_string = Dates.format(start_date, "dd/mm/yyyy HH:MM") * " (UTC)"
    p1_hnw = plot(times,heave, label="Heave", c="blue", lw=0.5, title=title_string, titlefontsize=12) ##last(split(infil,"\\")))
    p2_hnw = plot(times,north, label="North", c="red", lw=0.5)
    p3_hnw = plot(times,west, label="West", c="green", lw=0.5)

    hline!(p1_hnw, [0], lw=1, label="")
    hline!(p2_hnw, [0], lw=1, label="")
    hline!(p3_hnw, [0], lw=1, label="")

    # display plots to screen
    plot_wse = Plots.plot(p1_hnw, p2_hnw, p3_hnw, layout = (3, 1), size = (1400, 1300),
        xlim=(first(times),last(times)),  xticks = first(times):Minute(5):last(times),xtickfontsize=7,ytickfontsize=8,
        framestyle = :box,fg_legend=:transparent, legend=:bottomleft,
        margin = 15Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

    display(plot_wse)

    # create a plot file to be saved as a .PNG
    plt_file = first(infil, length(infil)-4)*'_'*Dates.format(start_date, "yyyy_mm_dd_HHMM")*".png"

    # Save plot to file
##            savefig(plt_file)
##            println("Plot file saved as ",plt_file)
    
    return(heave)
    
end    # plot_hnw()


function plot_f82(f82_df)
    # Plot current speed and direction from f82_df
    if !isempty(f82_df)
        p1_f82 = plot(f82_df.Date,f82_df.Speed,label="Speed (m/s)", c="red", lw=2, fillrange = 0, fillalpha = 0.025, fillcolor = :red, legend=:bottomleft,
            gridstyle=:dot, left_margin = 5Plots.mm, right_margin = 15Plots.mm, foreground_color_grid="red", yforeground_color_text="red", ytickfonthalign=:left)
        
        subplot = twinx()
        p1_f82 = plot!(subplot,f82_df.Date,f82_df.Direction,label="Direction", c="blue", line=:dot, lw=2, legend=:bottomright,
            gridstyle=:dashdot, ylim=(0,360), yflip=true, yticks = 0:30:360, showaxis=:y, margin = 15Plots.mm, foreground_color_grid="blue", yforeground_color_text="blue", ytickfonthalign=:right)
        
        plot_f82 = Plots.plot(p1_f82, layout = (1, 1), size = (1400, 800),
                xlim=(first(f82_df.Date),last(f82_df.Date)),  xticks = first(f82_df.Date):Hour(12):last(f82_df.Date), xtickfontsize=7, ytickfontsize=8,
                framestyle = :box, fg_legend=:transparent, title = "Current speed and direction",
                grid=true, gridlinewidth=0.5, gridalpha=1)            

            display(plot_f82)
    end
end


function plot_f29(f29_df)
# Plot H3 Hmax TH3 THmax as calculated by Datawell and stored to F29

    if !isempty(f29_df)
        p1_f29 = plot(f29_df.Date,f29_df.Hmax,label="Hmax", c="red", lw=2, fillrange = 0, fillalpha = 0.025, fillcolor = :red)
        p1_f29 = plot!(f29_df.Date,f29_df.H3,label="Hsig", c="blue", lw=2, fillrange = 0, fillalpha = 0.05, fillcolor = :blue)

        p2_f29 = plot(f29_df.Date,f29_df.THmax,label="THmax", c="red", lw=2, fillrange = 0, fillalpha = 0.025, fillcolor = :red)
        p2_f29 = plot!(f29_df.Date,f29_df.TH3,label="THsig", c="blue", lw=2, fillrange = 0, fillalpha = 0.05, fillcolor = :blue)
        
        # display plots to screen
        plot_f29 = Plots.plot(p1_f29, p2_f29,layout = (2, 1), size = (1400, 700),
            xlim=(first(f29_df.Date),last(f29_df.Date)),  xticks = first(f29_df.Date):Hour(12):last(f29_df.Date), xtickfontsize=7, ytickfontsize=8,
            framestyle = :box,fg_legend=:transparent, legend=:bottomleft,
            margin=1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

        display(plot_f29)
        
    else
    
        println("ALERT: No F29 data available!")

    end    
       
end
