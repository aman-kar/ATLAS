#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 10:47:09 2021


Overview: A modified tool the produces SAP, PDCSAP and unpopular lightcurves for a given target.


Description: A script that create lightcurves of targets using the unpopular/tess_cpm package
where flux from all neighbouring targets (pixels) are found and removed from the target pixel.
This is then overlayed with PDCSAP & SAP lightcurves from SPOC/TESS-SPOC to visually appreciate
the difference between the two pipelines.

All ignore comments imply the code will automatically ignore the given lines and you may simply
skip those lines and read the rest of the script to understand this code.

@author: Aman Kar
"""

from astropy.coordinates import Angle
from astroquery.mast import Tesscut
from astroquery.mast.utils import parse_input_location
from library import pretty_plot,dictfetchall
from astropy import units as u
from scipy.stats import iqr
from lomb_scargle import custom_LombScargle

import os, re
import tess_cpm
import psycopg2
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk

def check_before_download(coordinates=None, size=5, sector=None, path=".", inflate=True,
                          objectname=None, force_download=False):
    """ Function to download Tess FFI cutouts of given size.
    This will check if *.fits file has already been download previously.

    Parameters
    ----------
    force_download : bool, false (default)
        If set to True will result in redownloading the fiile. (Must be done manually)

    Returns
    -------
    list
        A list of path name to the location of the *.fits for every of given objectname/target
    """

    coords = parse_input_location(coordinates, objectname) #Locate coordinates for given objectname
    ra = f"{coords.ra.value:.6f}" #Parse out only ra
    matched = [path+m if path!='.' else m for m in os.listdir(path) if ra in m]
    if (len(matched) != 0) and (force_download == False):
        print(f"Found the following FITS files in the \"{path}/\" directory with matching RA values.")
        print(matched)
        print("If you still want to download the file, set the \'force_download\' keyword to True.")
        return matched
    else:
        #Download FFI cutout for the first time | size -> cutout size (e.g. 100 -> 100x100 FFI cutout)
        path_to_FFIs = Tesscut.download_cutouts(coordinates=coordinates, size=size,
                                                sector=sector, path=path, inflate=inflate, objectname=objectname)
        print(path_to_FFIs)
        return list(path_to_FFIs["Local Path"])

""" Setup for local PostgreSQL server connection
Ignore the following section if no database available
"""
print('Connecting to PostgreSQL table \n')
try:
    db = psycopg2.connect("dbname='postgres'")
    curs = db.cursor()

except:
    print("Unable to connect to the database")
    db = 0 #0 if no database connection achieved

#Store the target name and the RECONS/Display Name
print('Follow next prompt only if reducing a single target.\nOtherwise for ATLAS SQL database, enter blank for both.')
target_ = [{'tic_id' : input("Target Name? (TIC XXXXXXXXXXX) : "),
            'recons_name' : input("RECONS Name? (name to display) : ")}]

#PostgreSQL query to fetch target list from database.
#Ignore if executing just for single target.
if not target_[0].get('tic_id'):
    query  = ("select distinct nmpp.tic_id, av.\"Recons_Name\" as recons_name "
             "from atlas_variability av "
             "left join nasa_mdwarfs_planets_25pc nmpp on nmpp.hostname = av.\"Host Name\" "
             "where av.stddev is not NULL "
             "order by recons_name")
    curs.execute(query)
    target_ = dictfetchall(curs) #Store target list
    db.commit()

#Read in all previously saved apertures from apertures.dat
#Or create file if does not exists
apertures = {}
with open('apertures.dat','a+') as file:#Open file in append/read(+) mode
    file.seek(0) #Move stream position to the beginning of the file
    for line in file.readlines():
        t = line.split() #Format should be TIC_XXXX_S1 X1 X2 Y1 Y2
        apertures[t[0]] = list(map(int,t[1:])) #Maps each line onto a dict

#Loop through target list
for j in target_:

    print("Fetching lightcurves for {} ".format(j.get('recons_name') if j.get('recons_name') else j.get('tic_id')))

    #Initialize parameters
    target = j.get('tic_id') if j.get('tic_id') else j.get('recons_name')
    recons_name = j.get('recons_name')
    sap_stat = {"stddev":[],"idr":[]} #Array to store a list of stat for SAP FLux
    pdc_stat = {"stddev":[],"idr":[]} #Array to store a list of stat for PDCSAP FLux
    unpopular_stat = {"stddev":[],"idr":[]} #Array to store a list of stat for unpopular FLux
    sectors = [] #Array list to store every sector being processed
    stitch = False #Defaults to False for stitching lightcurves
    detail = False #If True, produces more detailed plots when performing the CPM Model Analysis
    cutout_size = 64  #(pixel x pixel)
    half_cutout = int(cutout_size/2)
    rowlims = [half_cutout - 5,half_cutout + 5] #X pixels for cutout
    collims = [half_cutout - 5,half_cutout + 5] #Y Pixels for cutout
    bad = 0 #Tracks number of bad sectors
    data = {} #Array containing the FINAL normalized unpopular timeseries data
    path_to_FFIs = sorted(check_before_download(size=cutout_size, #path='/Volumes/batcave/tess_cutouts/',#force_download=True,
                                                objectname=target)) #Fetch pathnames for each sector FFI stack

    print(f'\n{len(path_to_FFIs)} sectors found.\n') #Print number of sectors for which FFIs were found

    if len(path_to_FFIs) > 1 and input(f"Stitch lightcurves from all {len(path_to_FFIs)} sectors?(y/n) : ") == 'y':

        #Initialize stitching arrays
        stitch = True
        st_lc = []
        st_detrended = []
        st_time = []
        st2_time = []
        final_lc = []
        final_time = []

        print("\n")
        print(f"At the end of this script, a stitched lightcurve from all {len(path_to_FFIs)} sectors will be plotted.")

    for idx,i in enumerate(path_to_FFIs):

        plt.style.use('default')

        sector = int(re.search('-s([0-9]{4,4})-',i).group(1)) #Extract Sector number from path to FFI stack
        print(f'\nExtracting Sector {sector}')

        #Get TPF for the target and plot the optimal aperture mask from TESS-SPOC (PDCSAP data)
        tpf = lk.search_targetpixelfile(target,author='TESS-SPOC',mission='TESS',sector=sector).download()
        if tpf:
            tpf.plot(aperture_mask='pipeline', show_colorbar=False, style='bmh', title='Target Pixel File') #Plot TPF Aperture
            plt.show(block=False)
            #Get lightkurve object for given sector from TESS-SPOC (PDCSAP data)
            lc = lk.search_lightcurve(target,author='TESS-SPOC',mission='TESS',sector=sector).download()
        else:
            lc = False #No Data
            print(f'No TESS-SPOC data available for Sector {sector}')

        #Using the CPM method to find your target and download the .fits file locally
        s = tess_cpm.Source([s for s in path_to_FFIs if "{:04d}".format(sector) in s][0], remove_bad=True,bkg_subtract=True)
        s.plot_cutout(l=0, h=98,projection="wcs") #Visually inspect if the correct target has been chosen

        if tpf:
            print("Check if TessCut match with target Coordinates: \n"
                  "RA = {} | Dec = {} ".format(Angle(tpf.ra,unit=u.degree).to_string(unit=u.hour),
                                              Angle(tpf.dec,unit=u.degree).to_string(unit=u.degree))) #Double check field

        if input("Type N if sector appears corrupted else hit Enter to continue: ")=='N':#Bad Data Check/Skip
            bad+=1 #Increase counter by 1
            continue #Skip everything and proceed to next sector

        sectors += [sector] #Append to list of sectors being reduced
        _ = s.plot_cutout(rowlims = rowlims,
                          collims = collims,
                          l=0, h=98) #Plot a mini cutout of just cutoutsize/2 x cutoutsize/2 for visual

        #Create dict key in the format: TIC_XXXXXX_S1
        ap_name = target.replace(' ','_')+'_S'+str(sector)
        ap_input = False

        #Plot previous apertures (if exists) and ask user if its desired or needs replaced
        if ap_name in apertures.keys():
            ap_xrange,ap_yrange = apertures[ap_name][:2],apertures[ap_name][2:]

            #Set apertures
            s.set_aperture(collims=ap_xrange,rowlims=ap_yrange)

            #Plot aperture on mini cutout
            s.plot_cutout(show_aperture=True, rowlims = rowlims, collims = collims, l=0, h=98)
            ap_input = input('Use previous aperture?(y/n) : ')

        if ap_input == 'y':
            pass
        else:
            print('Previous Aperture Not Found/Desired')
            #Define rectangular aperture
            print("\nDraw a rectangular aperture where input values correspond to pixel position. (Must be atleast 2x2)\n")
            print("Example - To draw a 3x4 aperture from X = 32 to 34 and Y = 25 to 28, enter "\
                  "32 34 25 28 for the following 4 prompts.")

            while(True):
                try:
                    ap_xrange = [int(input('X Start? : ')),int(input('X End? : '))] #Aperture X length
                    ap_yrange = [int(input('Y Start? : ')),int(input('Y End? : '))] #Aperture Y length
                    break
                except:
                    print("Wrong Input. Try Again!\n")

            #Store apertures in dict that will be writen to apertures.dat
            apertures[ap_name] = ap_xrange + ap_yrange

            #Set apertures
            s.set_aperture(collims=ap_xrange,rowlims=ap_yrange)

            #Plot aperture on mini cutout
            s.plot_cutout(show_aperture=True, rowlims = rowlims, collims = collims, l=0, h=98)

        ###########################
        #Structuring the CPM Model#
        ###########################

        #Hyperparameters used for the CPM method
        L = 256 #Number of predictor pixels
        k = 150 #Number of contiguous sections to break the lightcurve for train-and-test
        lam_L = 0.01 #Regularization Term for the systematics model
        lam_P = 0.01 #Regularization Term for the polynomial model

        #For more details follow unpopluar tutorial at https://github.com/soichiro-hattori/unpopular/blob/master/intro-tess_cpm.ipynb

        s.plot_cutout(show_aperture=True, l=0, h=98) #plot aperture on FFI cutout
        #s.plot_pix_by_pix(data_type="normalized_flux") #Plots light curve grid for every predictor pixel
        s.add_cpm_model(exclusion_size=10, n=L, predictor_method="similar_brightness") #10+1 x 10+1 exclusion grid

        # This method allows us to see the predictor pixels chosen | Predictor pixel size - size of red dots
        #half_model_length = int(len(s.models)/2)
        #_ = s.models[half_model_length][half_model_length].plot_model(size_predictors=5)

        s.add_poly_model(scale=2, num_terms=4) #Add polynomial model to capture long-term trends

        s.set_regs([lam_L,lam_P])  # The first term is the CPM regularization while the second term is the polynomial regularization value.(lamL, lamP)
        s.holdout_fit_predict(k=k);  # When fitting with a polynomial component, we've found it's better to increase the number of sections.

        #s.plot_pix_by_pix(data_type="poly_model_prediction", split=True) #Plot pix by pixel
        #s.plot_pix_by_pix(data_type="cpm_subtracted_flux") #Plots pix by pix

        s_aperture_normalized_flux = s.get_aperture_lc(data_type="normalized_flux") #Normalized FFI Data
        s_aperture_cpm_prediction = s.get_aperture_lc(data_type="cpm_prediction")

        rescaled_cpm = s.get_aperture_lc(data_type="rescaled_cpm_subtracted_flux") #From Hattori email

        s_aperture_poly_prediction = s.get_aperture_lc(data_type="poly_model_prediction")
        s_aperture_detrended_flux = s.get_aperture_lc(data_type="cpm_subtracted_flux")

        if detail:

            pretty_plot(plt)
            plt.plot(s.time, s_aperture_normalized_flux, ".", c="k", ms=8, label="Normalized FFI Flux")
            plt.plot(s.time, s_aperture_cpm_prediction, "-", lw=3, c="C3", alpha=0.8, label="CPM Systematics Prediction")
            plt.plot(s.time, s_aperture_poly_prediction, "-", lw=3, c="C0", alpha=0.8, label="Polynomial Prediction")

            plt.xlabel("Time - 2457000 [Days]", fontsize=30)
            plt.ylabel("Normalized Flux", fontsize=30)
            plt.tick_params(labelsize=20)
            plt.legend(fontsize=30)
            plt.show(block=False)

            plt.plot(s.time, s_aperture_detrended_flux, "k-")
            plt.xlabel("Time - 2457000 [Days]", fontsize=30)
            plt.ylabel("CPM Subtracted Flux", fontsize=30)
            plt.tick_params(labelsize=20)
            plt.show(block=False)

        #Overlay unpopular with PDCSAP & SAP
        pretty_plot(plt)
        fig, ax = plt.subplots(1,1,figsize=(10,5))
        legend=[] #Legend for given sector's lightcurve

        if lc:#Executes only if TESS-SPOC data is available

            #SAP Flux -> Magnitude, Stddev, IDR
            lc = lc.remove_nans(column='sap_flux') #Remove emtpy values
            sap_flux = lc.sap_flux.value
            sap_stat['btjd'] = lc.time.value #BTJD for X axis
            sap_stat['delta_mag'] = (-2.5*np.log10(sap_flux/np.median(sap_flux)))*1e3 #magnitude formula directly used on normalized flux | units of mmag
            sap_stat['stddev'].append(np.std(sap_stat.get('delta_mag'))) #Stddev for SAP Flux
            sap_stat['idr'].append(iqr(sap_stat.get('delta_mag'),rng=(10,90))) #IDR for SAP Flux

            ax.plot(sap_stat.get('btjd'), sap_stat.get('delta_mag'),lw=0.5)
            legend.append(r'SAP Flux | $\sigma$ = {:.2f} mmag | IDR = {:.2f} mmag'.format(sap_stat.get('stddev')[-1],sap_stat.get('idr')[-1]))

            #PDCSAP Flux -> Magnitude, Stddev, IDR
            lc = lc.remove_nans(column='pdcsap_flux') #Remove emtpy values
            pdcsap_flux = lc.pdcsap_flux.value
            pdc_stat['btjd'] = lc.time.value #BTJD for X axis
            pdc_stat['delta_mag'] = (-2.5*np.log10(pdcsap_flux/np.median(pdcsap_flux)))*1e3 #magnitude formula directly used on normalized flux | units of mmag
            pdc_stat['stddev'].append(np.std(pdc_stat.get('delta_mag'))) #Stddev for PDCSAP Flux
            pdc_stat['idr'].append(iqr(pdc_stat.get('delta_mag'),rng=(10,90))) #IDR for PDCSAP Flux

            ax.plot(pdc_stat.get('btjd'), pdc_stat.get('delta_mag'),lw=0.5)
            legend.append(r'PDCSAP Flux | $\sigma$ = {:.2f} mmag | IDR = {:.2f} mmag'.format(pdc_stat.get('stddev')[-1],pdc_stat.get('idr')[-1]))

        else:#Otherwise skips two color cycles to maintain the same color code for unpopular flux
            next(ax._get_lines.prop_cycler) #Skip SAP and cycle to the next color
            next(ax._get_lines.prop_cycler) #Skip PDCSAP and cycle to the next color

        #unpopular Flux -> Magnitude, Stddev, IDR
        unpopular_stat['delta_mag'] = (-2.5*np.log10(rescaled_cpm/np.median(rescaled_cpm)))*1e3 #magnitude formula directly used on normalized flux | units of mmag
        unpopular_stat['stddev'].append(np.std(unpopular_stat.get('delta_mag')[~np.isnan(unpopular_stat.get('delta_mag'))])) #Stddev for PDCSAP Flux
        unpopular_stat['idr'].append(iqr(unpopular_stat.get('delta_mag')[~np.isnan(unpopular_stat.get('delta_mag'))],rng=(10,90))) #IDR for PDCSAP Flux
        ax.plot(s.time, unpopular_stat.get('delta_mag'),lw=1)
        legend.append(r'unpopular Flux | $\sigma$ = {:.2f} mmag | IDR = {:.2f} mmag'.format(unpopular_stat.get('stddev')[-1],unpopular_stat.get('idr')[-1]))

        ax.legend(legend)
        ax.set_title(f'{target} - Sector {sector}')
        ax.set_xlabel('Time - 2457000 [BTJD days]')
        ax.set_ylabel(r'$\Delta$ (mmag)')
        plt.gca().invert_yaxis()
        plt.show()

        #Plot a LombScargle Periodogram for the given sector ONLY (Courtesy - Andrew Couperus)
        custom_LombScargle({f'{target} Sector {sector}' : [s.time,unpopular_stat.get('delta_mag')]})

        print(f'------------------ {target} at Sector {sector} completed ------------------')

        #Dirty hack to stitch unpopular fluxes
        #Should be thoroughly examined if used for publication!
        if stitch: #Executes only if user inputs y at the beginning

            if idx - bad > 0:#If the second sector and beyond are being processed
                diff, coeff, st_time, st_lc = tess_cpm.utils.stitch_sectors(st_time, s.time, st_lc, rescaled_cpm) #CPM Stitching Method

            else: #If the first sector is being reduced
                st_lc = rescaled_cpm
                st_time = s.time

            final_lc.append(rescaled_cpm)
            final_time.append(s.time)
            data[f'{target} Sector {sector}'] = [s.time,unpopular_stat.get('delta_mag')]

    #Plot the stitched lightcurve and its LS periodogram
    if stitch:
        pretty_plot(plt)
        plt.title(f'{target}')
        plt.xlabel('Time - 2457000 [BTJD days]')
        plt.ylabel('Flux (e- / s)')
        plt.plot(st_time,st_lc,'.')
        custom_LombScargle(data,merge=True)


    #Insert/update database only if database connection was successful and RECONS name was found/entered
    if db and recons_name and input('Would you like to update PostGreSQL? (y/n) : ') == 'y':

        query  = ("UPDATE atlas_variability SET "
                  "pdcsap_stddev = %s,"
                  "pdcsap_idr = %s,"
                  "avg_pdcsap_stddev = %s,"
                  "avg_pdcsap_idr = %s,"
                  "unpopular_stddev = %s,"
                  "unpopular_idr = %s,"
                  "avg_unpopular_stddev = %s,"
                  "avg_unpopular_idr = %s,"
                  "tess_sectors = %s"
                 "WHERE \"Recons_Name\" = %s ")
        curs.execute(query,(pdc_stat.get('stddev')[0],
                            pdc_stat.get('idr')[0],
                            np.mean(pdc_stat.get('stddev')),
                            np.mean(pdc_stat.get('idr')),
                            unpopular_stat.get('stddev')[0],
                            unpopular_stat.get('idr')[0],
                            np.mean(unpopular_stat.get('stddev')),
                            np.mean(unpopular_stat.get('idr')),
                            str(sectors),
                            recons_name))
        db.commit()

#Truncates apertures file and write new 'updated' apertures dict
with open('apertures.dat','w') as file:
    for i in apertures.keys():
        #Format: TIC_XXXX_S1 X1 X2 Y1 Y2
        file.write(f'{i} {apertures[i][0]} {apertures[i][1]} '\
                   f'{apertures[i][2]} {apertures[i][3]}\n')

if db:
    db.close()



