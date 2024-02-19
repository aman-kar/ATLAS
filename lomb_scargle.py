import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle


"""
#Read in csv list for target info (target_name,TIC ID)
ids = np.loadtxt("RTWINS_TIC.txt", delimiter=',', dtype=str)

apertures_used = {}
stars_unpop_data = {}
stars_unpop_stitched_data = {}


#Loop over my stars, running Aman's unpopular analysis code

for star in ids:
    print("!!! --- Starting process for " + star[0] + " --- !!!")

...

    #Store the target name and the RECONS/Display Name
    target_ = [{'tic_id' : star[1],
                'recons_name' : star[0]}]

...

    #Loop through target list
    for j in target_:

        ...

        for idx,i in enumerate(path_to_FFIs):

            ...

            #Save aperture used for this star+sector combination
            nam = str(star[0]) + "_S" + str(sector)
            apertures_used[nam] = [ap_xrange, ap_yrange]

            ...

            #Save unpopular data for this star+sector combination: times and delta mags
            stars_unpop_data[nam] = [s.time, delta_mag]

            ...

        if stitch:
            ...

            #Save unpopular stitched data for this star: times and fluxes
            stars_unpop_stitched_data[star[0]] = [st_time,st_lc]

        ...


#Save the output unpopular data lists to doc: times and delta mags for star+sector combos

with open("unpop_data.txt", 'w') as f:
    for star in stars_unpop_data:
        f.write(star + '\n')
        for i in range(len(stars_unpop_data[star][0])):
            f.write(str(stars_unpop_data[star][0][i]) + "," + str(stars_unpop_data[star][1][i]))
            f.write('\n')

"""




#Read in unpopular data file for each star+sector combination
# =============================================================================
#
# unpop_dat = {}
# file1 = open('unpop_data.txt', 'r')
# Lines = file1.readlines()
# for line in Lines:
#     if line.strip() == '': #skip empty lines (.strip removes newline characters from the string too)
#         continue
#     elif ',' in line.strip(): #If has a comma in it, is a data line
#         unpop_dat[current_name][0].append( float( line.strip().split(',')[0] ) ) #times
#         unpop_dat[current_name][1].append( float( line.strip().split(',')[1] ) ) #mags
#     else: #so is a new name line, make a new star entry and set current star to new one and data follows
#         unpop_dat[line.strip()] = [[], []]
#         current_name = line.strip()
#
#
# =============================================================================

#Merge all sectors of data (delta mags) for each star into single data list for each star
#Note: I do not normalize any of the data myself here.
#Note: When they were converted to delta mags during Aman / unpopular code, that effectively normalized each sector to 0.



def custom_LombScargle(unpop_dat,merge=False,max_period=55):

    if merge:

        unpop_dat_SectorsMerged = {}

        for star_sect in unpop_dat:
            star = star_sect.split('Sector ')[0]
            sect = star_sect.split('Sector ')[1]

            if star not in unpop_dat_SectorsMerged.keys(): #star is NOT in dict yet, add entry for it
                unpop_dat_SectorsMerged[star] = [[],[],[]] #data lists in order of Time,Mag,Sector

            #Add data
            unpop_dat_SectorsMerged[star][0].extend(unpop_dat[star_sect][0]) #times
            unpop_dat_SectorsMerged[star][1].extend(unpop_dat[star_sect][1]) #mags
            unpop_dat_SectorsMerged[star][2].extend( [sect]*len(unpop_dat[star_sect][1]) ) #sector for data point

        unpop_dat = unpop_dat_SectorsMerged

    #Set plotting parameters
    plt.rc("font", size=14, weight='normal')
    plt.rc("axes", labelsize=16, titlesize=10)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.rc("legend", fontsize=10)
    plt.rcParams['font.family'] = 'serif'
    temp = "\ "



    #Generate Master Plot for each individual star+sector
    #! - Assumes delta magnitudes, not fluxes! (easy to adapt though)
    #! - Replace "unpop_dat" with "unpop_dat_SectorsMerged" and it runs the same with sectors data merged now - !#

    save = False
    show = True

    for star in unpop_dat:

        #!! --- Setup data lists --- !!#

        #! - Can change the input data lists here and everything else will proceed consistently
        times_tmp = unpop_dat[star][0].copy() #times in days
        mags_tmp = unpop_dat[star][1].copy() #delta magnitudes in mmag

        # - Remove any points *brighter* than 5x stdev of the data (large flares, bad points, etc)
        # - I don't remove things *dimmer* than 5x stdev, in order to keep in any large transits just in case
        stdmag = np.std(mags_tmp) #standard deviation of the delta magnitudes
        times = [times_tmp[i] for i in range(len(times_tmp)) if mags_tmp[i] > (-5*stdmag)]
        mags = [mags_tmp[i] for i in range(len(mags_tmp)) if mags_tmp[i] > (-5*stdmag)]


        #!! --- Plot light curve (first panel) --- !!#
        # Create 3-panel plot:
        # -first panel is normal light curve
        # -second panel is Lomb-Scargle periodogram
        # -third panel is phase-folded light curve

        plt.figure(0)
        plt.clf()
        fig, (ax) = plt.subplots(3,1,figsize=(9,9))

        a = 0.85 #set alpha for plotting
        p_size = 1 #set point size for plotting

        ax[0].plot(times, mags, 'ko' , alpha=a, markersize=p_size, zorder=1)

        #ticks
        ax[0].minorticks_on()
        ax[0].tick_params(axis='both', which='minor', right=True, top=True, direction='in')
        ax[0].tick_params(axis='both', which='minor', bottom=True, left=True, direction='in')
        ax[0].tick_params(axis='both', which='major', right=True, top=True, direction='in')

        #center line
        ax[0].axhline(np.mean(mags), color='black', linestyle = '-', linewidth=0.85, zorder=0, alpha=1)

        #limits
        ax[0].set_ylim(stdmag*-5, stdmag*5) #set based on stdmag so scales to the data. I've found x5 gives good bounds.
        ax[0].set_xlim(min(times)-0.5, max(times)+0.5)

        ax[0].invert_yaxis() #because magnitudes, invert so brighter is up

        #axis labels, text, title
        ax[0].set_xlabel("BTJD - 2457000 (days)", fontsize=15)
        ax[0].set_ylabel(r'$\Delta$ mmag', fontsize=15, labelpad=0)

        stdmag_new = np.std(mags) #get new standard deviation of data (after removing some outlier points in initial setup)
        ## ax[0].text(0.985, 0.85, r"$\sigma$ = " + '%.2f' % stdmag_new + " mmag", fontsize = 15, color='grey', ha='right', zorder=3, transform = ax[0].transAxes)

        elapsed = max(times) - min(times)
        ax[0].set_title(star + " | " + "unpopular" + " | " + "Timespan = " + '%.2f' % elapsed + "d", fontsize=11)


        #!! --- Compute and plot Lomb-Scargle Periodogram (second panel) --- !!#

        #! - Periodogram samples 100000 frequencies within corresponding period range of 0.05days to 200days
        #! - (can tweak as desired)
        freq_range = np.linspace(1/max_period, 1/0.05,100000) #must pass values as frequencies
        period_range = [1/f for f in freq_range]

        ls = LombScargle(times[:], mags[:], fit_mean=True) #if have error bars on individual data points, can pass them here to incorperate errors into the lomb-scargle model
        power = ls.power(freq_range)
        best_freq = freq_range[np.argmax(power)] #find frequency corresponding to the larget peak in LS periodogram
        best_fit_period = 1/best_freq #get corresponding period for the strongest peak in periodogram

        #Plot Periodogram in period space
        ax[1].plot(period_range, power, 'k-' ,linewidth=1.5, zorder=1)

        #axes
        ax[1].ticklabel_format(useOffset=False, style='plain') #prevents the offset scientific notation on plot labels
        ax[1].set_xscale('log') #plot periodogram in log scale for periods
        ax[1].set_xlabel('Period (Days)', fontsize=15, labelpad=0)
        ax[1].set_ylabel('LS Power', fontsize=15, labelpad=0)
        ax[1].set_xlim(min(period_range)*0.93, max(period_range)*1.07) #set plot bounds slightly above and below limits of points
        #let the ylim autoscale

        #vertical line at strongest peak in periodogram
        ax[1].axvline(best_fit_period, linewidth=1, color='red', zorder=2)

        #Compute false alarm probabilities. Plot lines at 5% and 1% FAP
        FAP = ls.false_alarm_probability(power.max())
        FAP_levels = ls.false_alarm_level([0.05, 0.01])
        ax[1].axhline(FAP_levels[0], linewidth=1, color='green', linestyle='--', zorder=0)
        ax[1].axhline(FAP_levels[1], linewidth=1, color='green', zorder=0)

        #text
        ## ax[1].text(0.03, 1.05, r"P$_{\rm rot}$ = " + '%.2f' % best_fit_period + "d", fontsize = 18, color='red', ha='left', zorder=3, transform = ax[1].transAxes)
        ax[1].text(0.97, 1.05, "FAP = " + '%.3f' % FAP, fontsize = 18, color='red', ha='right', zorder=3, transform = ax[1].transAxes)


        #!! --- Plot phase-folded light curve (third panel) --- !!#
        # - (use period from strongest peak in periodogram)

        #Compute phase-folded info
        per=best_fit_period
        n_phases = 1 #! - set to plot one-phase plot, 2-phase plot, etc
        x = times[:]
        y = mags[:]

        t_first = min(x) #arbitrarily set start of phase to time of first data point
        x_pf = [((n-t_first)/per) % n_phases for n in x] #the modulo trick with % returns just the decimal phase part as needed

        # - Set grid of values for plotting sin wave fit from lomb-scargle
        # If have large time range of data much longer than the period, then time_grid could get huge and slows things down
        # To avoid this, have time grid just span out to n_phases number of periods, so it only samples what is needed to plot phase diagram
        t_max = t_first + (n_phases*per)
        time_grid = np.arange(t_first,t_max, 0.05) #sample point every 0.05 days (! can set to 0.01 if short periods !)

        sin_fit_vals = ls.model(time_grid, 1/per) #need to give freq here. Gives values from sin model at grid points.
        sin_fit_phase = [((n-min(time_grid))/per) % n_phases for n in time_grid] #convert grid times to phases.
        sin_fit_phase, sin_fit_vals = zip(*sorted(zip(sin_fit_phase, sin_fit_vals))) #sorts lists to have points consistent and matching by phase

        #plot phase-folded data points and sin wave fit
        ax[2].plot(x_pf, y, 'ko' , alpha=a, markersize=p_size, zorder=1)
        ax[2].plot(sin_fit_phase, sin_fit_vals, 'r-', linewidth=3.5, zorder=1)

        #! - Compute and plot binned points in phase-folded light curve. (I do this in a very not efficient way...)
        bins = np.linspace(0, n_phases, 20 + 1) #do +1 to get bin at final edge proper
        binned_mags = []
        binned_phases = []
        binned_errs = []
        for i in range(len(bins)): #combine times and mags in each bin
            tmplistmags = []
            tmplisterrs = []
            tmplistphases = []
            if i == len(bins)-1: #check if at end of list, which is final bin edge at the end
                break
            for j in range(len(x_pf)): #loop over phase-folded data points
                if x_pf[j] >= bins[i] and x_pf[j] < bins[i+1]:
                    tmplistmags.append(y[j])
                    tmplistphases.append(x_pf[j])
                else:
                    continue
            if len(tmplistmags) == 0: #so no data points in this phase bin, then skip and go to next phase bin
                continue
            binned_mags.append(np.average(tmplistmags)) #if have error bars on each individual data point, can do weighted mean here
            binned_phases.append(np.mean(tmplistphases))
            binned_errs.append(np.std(tmplistmags)) #use the stdev for scatter of data points in bin as the error on the bin point

        ax[2].errorbar(binned_phases, binned_mags, yerr=binned_errs, fmt='sc', capsize=0, alpha=1.0, markersize=5.5, zorder=3)

        #ticks
        ax[2].minorticks_on()
        ax[2].tick_params(axis='both', which='minor', right=True, top=True, direction='in')
        ax[2].tick_params(axis='both', which='minor', bottom=True, left=True, direction='in')
        ax[2].tick_params(axis='both', which='major', right=True, top=True, direction='in')

        #center line
        ax[2].axhline(np.mean(mags), color='black', linestyle = '-', linewidth=0.85, zorder=0, alpha=1)

        #limits: sometimes the sin wave fit spans to values much beyond that of the data points (for poor fits) so handle that.
        #Default is to set bounds to +-5 * stdev of the data (similar limits used before for plotting and outlier removal)
        if max(sin_fit_vals) > (stdmag*5):
            y_up_lim = max(sin_fit_vals) + (0.07*max(sin_fit_vals))
        else:
            y_up_lim = stdmag*5

        if min(sin_fit_vals) < (stdmag*-5):
            y_low_lim = min(sin_fit_vals) + (0.07*min(sin_fit_vals))
        else:
            y_low_lim = stdmag*-5

        ax[2].set_ylim(y_low_lim, y_up_lim)
        ax[2].set_xlim(-(0.025*n_phases), n_phases*1.025)

        ax[2].invert_yaxis() #because magnitudes, invert so brighter is up

        #axis labels, text
        ax[2].set_xlabel('Phase', fontsize=15)
        ax[2].set_ylabel(r'$\Delta$ mmag', fontsize=15, labelpad=0)
        peak2peakDelta = max(sin_fit_vals) - min(sin_fit_vals)
        ## ax[2].text(0.97, 1.05, r"$\Delta$ = " + '%.2f' % peak2peakDelta + " mmag", fontsize = 18, color='red', ha='right', zorder=2, transform = ax[2].transAxes)
        ax[2].text(0.03, 1.05, r"P$_{\rm rot}$ = " + '%.2f' % best_fit_period + " d", fontsize = 18, color='red', ha='left', zorder=2, transform = ax[2].transAxes)


        plt.tight_layout() #makes all the plot peices fit together and display nicely

        #save plot / wrap-up
        # savename =  "plots\master_unpopular" + temp[0] + star + "_" + "unpopular" + "_MASTER" + ".png"
        # if save: plt.savefig(savename, format='png', dpi=250)
        if show: plt.show()
        plt.close(fig)

    print("-----Done-----")



"""
#! - For merged sectors version, can add text to the title listing the unique sectors included - !#

    #get list of included unique sectors
    unique_sectors = list( np.unique( np.array(unpop_dat_SectorsMerged[star][2]) ) )
    sec_string = ""
    for s in unique_sectors:
        sec_string = sec_string + s + " "

    ax[0].set_title(star + " | " + "unpopular" + " | " + "Timespan = " + '%.2f' % elapsed + "d" + '\nSectors Included: ' + sec_string[:-1], fontsize=11)
"""