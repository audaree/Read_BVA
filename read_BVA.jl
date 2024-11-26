## Julia program to read a selected .BVA file and display 30-minute time series plots
## JW October 2022

using Dates
using DataFrames
using DSP
#using Gtk
using LaTeXStrings
using NativeFileDialog
using Plots
using Statistics

include("./read_BVA_processing_tools.jl")
include("./read_BVA_plotting_tools.jl")

################################################
################################################
##           START OF MAIN PROGRAM
################################################
################################################

# Widen screen for better viewing
##display("text/html", "<style>.container { width:100% !important; }</style>") 

# Select a HVA daily .CSV file
infil = pick_file("C:\\QGHL\\Wave_data\\Bris\\BVA\\", filterlist="*BVA");
println("Selected ",infil)
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

# find all occurrences of 0x7e in packet vector
aa = findall(x->x.==0x7e, vec(packet));

# Create the df's to hold the processed data
f20_vals = []; f21_vals = []; f23_vals = []; f25_vals = []; f26_vals = []; f28_vals = []; f29_vals = [];
    f80_vals = []; f81_vals = []; f82_vals = []; fc1_vals = []; fc3_vals = []

f20_df = DataFrame(Date = [], Segments = [], Smax = [])
for i in 0:99
    col_name = "S$i"
    f20_df[!,col_name] = []
end
f21_df = DataFrame(Date = [], Segments = [])

for i in 0:99
    col_name = "Dir$i"
    f21_df[!,col_name] = []
end
for i in 0:99
    col_name = "Spread$i"
    f21_df[!,col_name] = []
end
f23_df = DataFrame(Date = [], Segments = [], Match_vector = [], Sample_number = [])
f25_df = DataFrame(Date = [], Segments = [], Hs = [], Ti = [], Te = [], T1 = [], Tz = [], T3 = [], Tc = [], Rp = [], Tp = [], Smax = [], Theta_p = [], Sigma_p = [])
f26_df = DataFrame(Date = [], Hmax = [], Thmax = [], Tmax = [], Htmax = [], Havg = [], Tavg = [], Hsrms = [], Nw = [], Nc = [], Epsilon = [], Coverage = [])
f28_df = DataFrame(Date = [], Segments = [])
for i in 0:99
    col_name = "m2_$i"
    f28_df[!,col_name] = []
end
for i in 0:99
    col_name = "n2_$i"
    f28_df[!,col_name] = []
end
for i in 0:99
    col_name = "k$i"
    f28_df[!,col_name] = []
end
f29_df = DataFrame(Date = [], Coverage = [], Nw = [], Epsilon = [], Hmax = [], THmax = [], H10 = [], TH10 = [], H3 = [], TH3 = [], Havg = [], Tavg = [])
for i in 0:22
    col_name = "Hq$i"
    f29_df[!,col_name] = []
end
f80_df = DataFrame(Date = [], Latitude = [], Longitude = [])
f81_df = DataFrame(Date = [], SST = [])
f82_df = DataFrame(Date = [], Firmware_Version = [], Speed = [], Direction = [], SST = [])
fc1_df = DataFrame(Date = [], Firmware = [], Hatch_uid = [], Hull_uid = [], Uptime = [], Battery_energy = [], Boostcaps_energy = [],
    Hatch_temp = [], Battery_voltage = [], Batteries_per_section = [], Battery_section_number = [], Initial_battery_energy = [], Ov = [], Cv = [],
    Ox = [], Oy = [], Cx = [], Cy = [], Mu0 = [], Sigma0 = [], Mui = [], Sigmai = [], Muh = [], Sigmah = [], Cpitch = [], Croll = [], Tensor = [])
fc3_df = DataFrame(Date = [], Battery_life = [])

# determine number of records
max_val = length(aa)-1

# Decode the packet data to messages
# refer to Section 2.1.2 Decoding the packet data to messages p. 20
messages = []

println("Processing the Packet vectors")
for i in 1:max_val
    # determine packet length
    first = aa[i]+1
    last = aa[i+1]
    
    if (last-first > 1)
        decoded = []
        decoded = packet[first:last-1]
        
        bb = findall(x->x.==0x7d, vec(decoded));
            
        if bb != []

            # do an xor of elements with 0x7d
            for ii in bb
                decoded[ii+1] = decoded[ii+1] ⊻ 0x20 # set the xor value as 0x20 vide 2.1.2 p.20
            end

            # remove the 0x7d
            deleteat!(decoded::Vector, bb)

        end
        
        # look for vectors of the spectrum synchronisation message (0xF23)
        if decoded[2] == 0x20
            
            heave_spectrum = []
            append!(f20_vals,decoded)
            timestamp,segments,smax,heave_spectrum = process_f20(decoded,heave_spectrum)         
            list_1 = [timestamp,segments,smax]
            push!(f20_df, [list_1; heave_spectrum])
            
        elseif decoded[2] == 0x21
            
            direction = []
            spread = []
            append!(f21_vals,decoded)
            timestamp,segments,direction,spread = process_f21(decoded,direction,spread)
                        
            list1 = [timestamp,segments]
            list2 = [direction; spread]
            
            push!(f21_df, [list1; list2])
            
        elseif decoded[2] == 0x23
                  
            append!(f23_vals,decoded)
            timestamp,segments_used,match_vector,sample_number = process_f23(decoded)
            push!(f23_df, [timestamp,segments_used,match_vector,sample_number])

        elseif decoded[2] == 0x25

            append!(f25_vals,decoded)
            timestamp,segments,hs,ti,te,t1,tz,t3,tc,rp,tp,smax,theta_p,sigma_p = process_f25(decoded)
            push!(f25_df, [timestamp,segments,hs,ti,te,t1,tz,t3,tc,rp,tp,smax,theta_p,sigma_p])
          
        elseif decoded[2] == 0x26
                  
            append!(f26_vals,decoded)
            timestamp,hmax,thmax,tmax,htmax,havg,tavg,hsrms,nw,nc,epsilon,coverage = process_f26(decoded)
            push!(f26_df, [timestamp,hmax,thmax,tmax,htmax,havg,tavg,hsrms,nw,nc,epsilon,coverage])

        elseif decoded[2] == 0x28
            
            m2 = []
            n2 = []
            k = []
            append!(f28_vals,decoded)
            timestamp,segments,m2,n2,k = process_f28(decoded,m2,n2,k)
                        
            list1 = [timestamp,segments]
            list2 = [m2; n2; k]
            
            push!(f28_df, [list1; list2])
            
        elseif decoded[2] == 0x29
            hq = []
            append!(f29_vals,decoded)
            timestamp,coverage,nw,epsilon,hmax,thmax,h10,th10,h3,th3,havg,tavg,hq = process_f29(decoded,hq)         
            list_1 = [timestamp,coverage,nw,epsilon,hmax,thmax,h10,th10,h3,th3,havg,tavg]
            push!(f29_df, [list_1; hq])

        elseif decoded[2] == 0x80
            
            append!(f80_vals,decoded)
            timestamp,latitude,longitude = process_f80(decoded)
            push!(f80_df, [timestamp,latitude,longitude])
            
        elseif decoded[2] == 0x81
            
            append!(f81_vals,decoded)
            timestamp,sst = process_f81(decoded)
            push!(f81_df, [timestamp,sst])
            
        elseif decoded[2] == 0x82
            
            append!(f82_vals,decoded)
            timestamp,firmware_version,speed,direction,sst = process_f82(decoded)
            push!(f82_df, [timestamp,firmware_version,speed,direction,sst])
      
        elseif decoded[2] == 0xc1
            
            append!(fc1_vals,decoded)
            timestamp,firmware,hatch_uid,hull_uid,uptime,battery_energy,boostcaps_energy,hatch_temp,battery_voltage,batteries_per_section,
                battery_section_number,initial_battery_energy,ov,cv,ox,oy,cx,cy,mu0,sigma0,mui,sigmai,muh,sigmah,cpitch,croll,tensor = process_fc1(decoded)   
            
            push!(fc1_df, [timestamp,firmware,hatch_uid,hull_uid,uptime,battery_energy,boostcaps_energy,hatch_temp,battery_voltage,batteries_per_section,
                battery_section_number,initial_battery_energy,ov,cv,ox,oy,cx,cy,mu0,sigma0,mui,sigmai,muh,sigmah,cpitch,croll,tensor])

        elseif decoded[2] == 0xc3
            
            append!(fc3_vals,decoded)
            timestamp,ble = process_fc3(decoded)
            push!(fc3_df, [timestamp,ble])

         end
        
    end
    
end
    
# remove duplicates from dataframes
f20_df = unique(f20_df)
f21_df = unique(f21_df)
f23_df = unique(f23_df)
f25_df = unique(f25_df)
f26_df = unique(f26_df)    
f28_df = unique(f28_df)
f29_df = unique(f29_df)
f80_df = unique(f80_df)
f81_df = unique(f81_df)
f82_df = unique(f82_df)
fc1_df = unique(fc1_df)
fc3_df = unique(fc3_df)

## Calculate the Heave, North, and West displacements
Data = []

# Convert binary data to hexidecimal vectors
local j = 0
println("Building displacements vectors - this takes a while!")
while true

    local alternating
    
    try
        local aa = []
        
        for i = j*12+1:j*12+12
            push!(aa,string(data[i], base = 16, pad = 2))
        end
        
        push!(Data,join(aa)[1:18])
        
    catch
        
        # escape is something is amiss        
        break
        
    end
    j = j+1
    
end

println("All file data read!")

# remove those vectors from F23 df that are not located in the Data vector df
f23_df[!,"Data_vector"] = [findfirst(x->x==i, Data) for i in f23_df.Match_vector];

# create a vector of dates from the F23 df
vector = Dates.format.(f23_df.Date, "yyyy-mm-dd HH:MM:SS");
println("Preparing to plot data")

# Do time-series plot of available data
plot_f29(f29_df)

# Plot current speed and direction
plot_f82(f82_df)
################################################
################################################
##           END OF MAIN PROGRAM
################################################
################################################

## Allow user to get a 30-minute record and do plots
using Gtk

cb = GtkComboBoxText()
choices = vector

for choice in choices
    push!(cb,choice)
end

set_gtk_property!(cb,:active,1)

signal_connect(cb, "changed") do widget, others...

    # get the active index
    idx = get_gtk_property(cb, "active", Int) + 2
  
    # get the active string 
    str = Gtk.bytestring( GAccessor.active_text(cb) ) 
    
    heave = plot_hnw(f23_df,Data,idx)
    plot_spectra(f23_df,heave,idx)
    plot_3d(f23_df,Data,idx)

end

win = GtkWindow("Select Date",200,200)
Gtk.GAccessor.position(win, Gtk.GtkWindowPosition.CENTER)
push!(win, cb)
showall(win)
