function string2hex(str)
################################################    
    parsed_val = (parse(Int, str[1], base=16) * 16 + parse(Int, str[2], base=16));
    if parsed_val == 0
        hex = 0x00
    else
        hex = [(parsed_val>>((i-1)<<3))%UInt8 for i in 1:sizeof(parsed_val)-leading_zeros(parsed_val)>>3][1]
    end
    
    return(hex)
    
    end    # string2hex()


function handle_gaps(df)
################################################    
# a function to identify where gaps found in sequence numbers
#    where gaps occur a dummy record is inserted into the df
#    and both Status1 and Status2 are flagged with a '!' 
    
    nums = []

    # Convert sequence number value from Hex to Integer
    for i in 1:nrow(df)
        str = df.Sequence[i]
        push!(nums,parse(Int, str[1], base=16) * 16 + parse(Int, str[2], base=16))
    end

    # Determine number of gaps in df rows
    counter = diff(nums)
    gaps = findall(counter.<1)
    counter[gaps] .+= 256;

    # now find where a counter is > 1 indicating number of gaps in transmission
    gaps = findall(counter.>1)
    ll = length(gaps)
        
    if ll > 0
        
        println(cumsum(gaps)," found in file!")
        println("File contains ",nrow(df)," records")
    

##    df1 = deepcopy(df)

        # need to work through the gaps in reverse order in order to preserve row numbers
        for i in ll:-1:1
            # for each gap, get sequence number of last valid row
            sequence_number = df[gaps[i],:].Sequence

            # for number of gaps, insert "FFF" values into df and '!' into the status indicators
            for j in 1:counter[gaps[i]]-1

                sequence_number_hex = uppercase(string((parse(Int, sequence_number[1], base=16) * 16 + parse(Int, sequence_number[2], base=16)) + j, base=16, pad=2))
        ##        println(i,' ',j,' ',gaps[i]+j,' ',counter[gaps[i]],' ',sequence_number,' ',sequence_number_hex,' ',"!",' ',"FFFFFFFFFFFFFFFFFF",' ',"!",' ',"FFFFFF")
                insert!(df, gaps[i]+j, [sequence_number_hex, '!', "FFFFFFFFFFFFFFFFFF", '!', "FFFFFF"])

            end
                
        end

        println("File now contains ",nrow(df)," records")
                
    else
                
        println("No gaps in record")
                
    end

    return (df)
        
    end    # handle_gaps()


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


function get_start_end_dates(f23_df,found_list)   
    start_date = f23_df[found_list[1],:].Date - Minute(30) # <------- NOTE subtracted 30min from start_date to match Waves4 results
    segments = f23_df[found_list[1],:].Segments
#   match_vector = f23_df[found_list[1],:].Match_vector
    sample_nos = f23_df[found_list[1],:].Sample_number
    data_vector = f23_df[found_list[1],:].Data_vector
    start_val = data_vector - Int(sample_nos/2) + 1
    end_val = data_vector
    
    return(start_date,start_val, end_val)
    
end    #(get_start_end_dates)
    
    
function get_displacement(Data, start_val, end_val)
################################################
# Decode the real time data to displacements - See DWTP (16 Jan 2019) 2.1.1 p. 19    
    
    arry = collect(Iterators.flatten(zip(SubString.(Data, start_val, end_val),SubString.(Data, start_val+9, end_val+9))));
    displacements = [parse(Int, SubString.(i, 1, 1), base=16)*16^2 + parse(Int, SubString.(i, 2, 2), base=16)*16^1 + parse(Int, SubString.(i, 3, 3), base=16)*16^0 for i in arry]    
    
    displacements[findall(>=(2048), displacements)] = displacements[findall(>=(2048), displacements)] .- 4096;
    displacements = 0.457*sinh.(displacements/457)    # see DWTP p.19 (16)
   
    return(displacements)
    
end    # get_displacement()


function process_f20(f20_vals,heave_spectrum)
#######################################
# function to calculate parameters from Upcross wave height quantiles message  (0xf20)
# refer to DWTP (Ver. 16 January2019) Section 4.8 pp.56-61
    
    # get Timestamp in UTC - refer Section 3.2 HF link header pp. 25-26
    timestamp = unix2datetime.(parse(Int, bitstring(f20_vals[3])*bitstring(f20_vals[4])*bitstring(f20_vals[5])*bitstring(f20_vals[6]); base=2));
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(0)

    # get Number of Segments Used
    segments = parse(Int, bitstring(f20_vals[9]); base=2); 
    
    # get Smax
    smax = 5000 * (exp(parse(Int, bitstring(f20_vals[10]) * bitstring(f20_vals[11])[1:4]; base=2) / 200) - 1) / (exp(4094/200) - 1);    # (45) p.37
    # Obtain the heave_spectrum from s0 to s99
    for ii in 12:3:159
        try
            push!(heave_spectrum,(exp(parse(Int, bitstring(f20_vals[ii]) * bitstring(f20_vals[ii+1])[1:4]; base=2) / 200) - 1) / (exp(4094/200) - 1))    # (48) p.38
            push!(heave_spectrum,(exp(parse(Int, bitstring(f20_vals[ii+1])[5:8] * bitstring(f20_vals[ii+2]); base=2) / 200) - 1) / (exp(4094/200) - 1))
        catch
            push!(heave_spectrum,missing)
            push!(heave_spectrum,missing)
        end
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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
        try
            push!(hq,parse(Int, bitstring(f29_vals[ii])[5:8] * bitstring(f29_vals[ii+1]); base=2) / 100;)    # (92)
            push!(hq,parse(Int, bitstring(f29_vals[ii+2]) * bitstring(f29_vals[ii+3])[1:4]; base=2) / 100;)
        catch
            push!(hq,missing)
            push!(hq,missing)
        end
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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

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
    timestamp = timestamp + Hour(0)

    # get Firmware Version Number
    mybytes = UInt8.(fc1_vals[9:18]);
    firmware = strip(String(mybytes),['\0'])
    
    # get Hatch UID
#    hatch_uid = parse(Int, bitstring(fc1_vals[24]) * bitstring(fc1_vals[23]) * bitstring(fc1_vals[22])
#        * bitstring(fc1_vals[21]) * bitstring(fc1_vals[20]) * bitstring(fc1_vals[19]); base=2) 
    hatch_uid = string(fc1_vals[19], base=16, pad=2) * string(fc1_vals[20], base=16, pad=2) * 
    string(fc1_vals[21], base=16, pad=2) * string(fc1_vals[22], base=16, pad=2) * 
    string(fc1_vals[23], base=16, pad=2) * string(fc1_vals[24], base=16, pad=2)    
    # get Hull UID
#    hull_uid = parse(Int, bitstring(fc1_vals[30]) * bitstring(fc1_vals[29]) * bitstring(fc1_vals[28])
#        * bitstring(fc1_vals[27]) * bitstring(fc1_vals[26]) * bitstring(fc1_vals[25]); base=2) 
    hull_uid = string(fc1_vals[25], base=16, pad=2) * string(fc1_vals[26], base=16, pad=2) * 
    string(fc1_vals[27], base=16, pad=2) * string(fc1_vals[28], base=16, pad=2) * 
    string(fc1_vals[29], base=16, pad=2) * string(fc1_vals[30], base=16, pad=2)    
    # get Uptime
    uptime = parse(UInt, bitstring(fc1_vals[31]) * bitstring(fc1_vals[32]) * bitstring(fc1_vals[33]) * bitstring(fc1_vals[34]); base=2) / 3600 / 24  # time in days
    
    # get Energy used from batteries
    battery_energy = parse(Int, bitstring(fc1_vals[35]) * bitstring(fc1_vals[36]) * bitstring(fc1_vals[37]) * bitstring(fc1_vals[38]); base=2) 

    # get Energy to boostcaps
    boostcaps_energy = parse(Int, bitstring(fc1_vals[39]) * bitstring(fc1_vals[40]) * bitstring(fc1_vals[41]) * bitstring(fc1_vals[42]); base=2)

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
    timestamp = timestamp + Hour(0)

    # get ble
    ble = parse(Int, bitstring(fc3_vals[9]); base=2);

##    println(timestamp,"UTC ",' ',ble)
    
    return(timestamp,ble)
    
end    # process_fc3()


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

            # truncate if too long
            wse = wse[1:4608]
            
        else

            # zero pad if too short (leave it unchanged if right length)
            append!(wse,zeros(4608-wse_length))
            
        end      

    end
    
    return (heave, north, west)
    
end    # get_hnw()


function calc_hm0(Sf,freq)
##########################################    
ax1 = (last(freq) - first(freq)) / (length(freq)-1)

# calc spectral moments m0, m1, m2, m3, and m4
s00 = 0; m0 = 0

for ii in 1:128

    s00 += f2[ii]^0 * Sf[ii];

end

m0 = 0.5*ax1*(first(f2)^0*first(Sf) + 2*s00 + last(f2)^0*last(Sf))

return(4 * m0^0.5)

end    # calc_hm0()


function calc_representative_spectra(frequency,Hm0,Tp,gamma)
##########################################    
    """
    function to calculate representative spectrum based on the Jonswap formula in Tucker and Pitt p.339 (10.3-9a)

    inputs:
        frequency - array of spectral frequencies
        Hm0 - floating point value
        Tp - floating point value
        gamma - floating point value - Peak ehhancement factor (Enter 1 for PM, or 3.3 for Jonswap)

        Typical calls:
        Spectra_PM = calc_representative_spectra(f2, Hm0, Tp, 1.0)
        Spectra_JONSWAP = calc_representative_spectra(f2, Hm0, Tp, 3.3)

    returns:
        Sf - array of representative spectra        
    """

    alpha = 1    # initial Philips constant (will decrease for each iteration required)
    g = 9.81
    fp = 1/Tp    # peak frequency

    hm0 = 99.    # set this to large value (so it will change on first iteration)

    Sf = [];

    while((Hm0 - hm0) <= 0.0005)

        Sf = vcat([alpha*g^2 * (2*pi)^-4 * ff^-5 * exp(-1.25 * (ff/fp)^-4) * gamma^exp(-(ff-fp)^2/(2*0.07^2 * fp^2)) for ff in frequency[findall(<=(fp), frequency)]],
                [alpha*g^2 * (2*pi)^-4 * ff^-5 * exp(-1.25 * (ff/fp)^-4) * gamma^exp(-(ff-fp)^2/(2*0.09^2 * fp^2)) for ff in frequency[findall(>(fp), frequency)]]);
        Sf[1] = 0;

###################################################################################################################################            
###  See discussion at https://stackoverflow.com/questions/44915116/how-to-decide-between-scipy-integrate-simps-or-numpy-trapz  ###
###################################################################################################################################

        ##        hm0 = 4*(np.trapz(Sf, frequency))^0.5    # calculate new Hm0 based on Sf values
        hm0 = calc_hm0(Sf,frequency);   # see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.integrate.simps.html        
        alpha *= 0.95;    # reduce alpha by 5% so iterations approach a limit of 0.0005

    end

    return(Sf)    # calc_representative_spectra()
        
    end    # calc_representative_spectra()


function calc_tp5(f2,Sf)
##########################################
# Calculate Tp5 via Read method
    Sf_max = maximum(Sf)

    numerator = 0; denominator = 0

    Sf_sum = cumsum(Sf.*Sf_max).^5

    for i in 1:length(f2)
        w = Sf[i] / Sf_max
        numerator +=  f2[i] * w^5
        denominator += w^5
    end

    Fp5 = numerator / denominator
    
    return(1/Fp5)    # calc_tp5()

    end    # calc_tp5()


function calculate_frequency_domain_parameters(f2, spectra)
##########################################
# Calculate frequency-domain parameters    
# Calls: calc_tp5()
    ax1 = (last(f2) - first(f2)) / (length(f2)-1)

    # calc spectral moments m0, m1, m2, m3, and m4
    s00 = 0; s01 = 0; s02 = 0; s03 = 0; s04 = 0;
    m0 = 0; m1 = 0; m2 = 0; m3 = 0; m4 = 0

    for ii in 1:128

        s00 += f2[ii]^0 * spectra[ii]
        s01 += f2[ii]^1 * spectra[ii]
        s02 += f2[ii]^2 * spectra[ii]
        s03 += f2[ii]^3 * spectra[ii]
        s04 += f2[ii]^4 * spectra[ii]

    end

    m0 = 0.5*ax1*(first(f2)^0*first(spectra) + 2*s00 + last(f2)^0*last(spectra))
    m1 = 0.5*ax1*(first(f2)^1*first(spectra) + 2*s01 + last(f2)^1*last(spectra))
    m2 = 0.5*ax1*(first(f2)^2*first(spectra) + 2*s02 + last(f2)^2*last(spectra))
    m3 = 0.5*ax1*(first(f2)^3*first(spectra) + 2*s03 + last(f2)^3*last(spectra))
    m4 = 0.5*ax1*(first(f2)^4*first(spectra) + 2*s04 + last(f2)^4*last(spectra))

    ##println("m0 = ",m0," m1 = ",m1, " m2 = ",m2, " m3 = ",m2, " m4 = ",m4)

    # calc wave parameters Hm0, Hrms, T01, T02, Tc
    Hm0 = 4*sqrt(m0)     # Tucker & Pitt p.32 (2.2-6b)
    Hrms = sqrt(8*m0)    # Goda 2nd. Edition p.262 (9.15)
    T01 = m0/m1          # Tucker & Pitt p.41 Table 2.2 
    T02 = sqrt(m0/m2)    # Tucker & Pitt p.40 (2.3-2)
    Tc = sqrt(m2/m4)     # Tucker & Pitt p.41 Table 2.2 - also see Notes

    # identify spectral peak and frequency as peak
    Fp = f2[argmax(spectra)]
    Tp = 1/Fp
    Tp5 = calc_tp5(f2, spectra)

    # calculate spectral width vide Tucker and Pitt p.85 (5.2-8)
    # Note: for JONSWAP, v = 0.39; for PM, v = 0.425
    v = (m0*m2 / m1^2 - 1)^0.5

    # calculate Skewness vide Tucker and Pitt p.109 (5.5-17)
    Skewness = (m0^2 * m3/m1^3 - 3*v^2 - 1) / v^3;
    
    return(Hm0, Hrms, T01, T02, Tc, Tp, Tp5, Skewness)
    
    @printf("%s; Hm0 = %5.2fm; Hrms = %5.2fm; T01 = %5.2fs; T02 = %5.2fs; Tc = %5.2fs; Tp = %5.2fs; Tp5 = %5.2fs; Skewness = %5.4f",
        Dates.format(start_date, "yyyy-mm-dd HH:MM"),Hm0, Hrms, T01, T02, Tc, Tp, Tp5, Skewness)
    
    end    # calculate_frequency_domain_parameters()


function calc_hm0(Sf,freq)
##########################################    
    ax1 = (last(freq) - first(freq)) / (length(freq)-1)

    # calc spectral moments m0, m1, m2, m3, and m4
    s00 = 0; m0 = 0

    for ii in 1:128

        s00 += freq[ii]^0 * Sf[ii];

    end

    m0 = 0.5*ax1*(first(freq)^0*first(Sf) + 2*s00 + last(freq)^0*last(Sf))

    return(4 * m0^0.5)

    end    # calc_hm0()

