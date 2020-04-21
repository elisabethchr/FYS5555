import numpy as np
import ROOT
from array import array
import sys
import os
from os import listdir
from os.path import isfile, join
from datetime import datetime
import math
import statistics
import glob
ROOT.TH1.SetDefaultSumw2()

# Function to retrieve difference in seconds between two time stamps
def getDuration(start,stop):
    timestamp_start = "%d %06d" %(t.GetDate(True,start),t.GetTime(True,start))
    timestamp_stop  = "%d %06d" %(t.GetDate(True,stop),t.GetTime(True,stop))

    datetime_object_start = datetime.strptime(timestamp_start, '%Y%m%d %H%M%S')
    datetime_object_stop  = datetime.strptime(timestamp_stop,  '%Y%m%d %H%M%S')
    difference = datetime_object_stop - datetime_object_start
    seconds_in_day = 24 * 60 * 60
    if difference.days < 0:
        return difference.seconds
    else:
        return difference.days * seconds_in_day + difference.seconds


# Folder where the data is stored
indir = sys.argv[1]
directories = next(os.walk(indir))[1]

dir = []
for i in sorted(directories):
    dir.append(indir+i)


# The three trees containing different data
c_trend = ROOT.TChain("Trending")
c_headr = ROOT.TChain("Header")
c_weath = ROOT.TChain("Weather")

# Define time span (default: defined as the expedition with Nanuq (PolarquEEEst)
start_time = ROOT.TTimeStamp(2018,7,21,0,0,0)
stop_time  = ROOT.TTimeStamp(2018,9,5,0,0,0)

x_rawrate_all = []
y_rawrate_all = []
y_rawrate_p_all = []    # raw rate as function of pressure
y_rawrate_tin_all = []  # raw rate as function of indoor temperature
y_rawrate_tout_all = [] # raw rate as function of outdoor temperature


x_event_all = []
y_pressure_all = []
y_temp_in_all = []
y_temp_out_all = []

# Loop over files in input directory and select those
# which are inside the time span defined above
for k in xrange(len(dir)):
    indir = dir[k]
    detector = indir.split('/')[-1]
    print "detector = ", detector
    onlyfiles = [f for f in listdir(indir) if isfile(join(indir, f))]
    files_to_use = []
    for fname in sorted(onlyfiles):
        if not fname.endswith("summary.root"): continue

        # Get start date and stop date from file name
        st_y,st_m,st_d = (fname.split("_")[1]).split("-")
        start_dt = ROOT.TTimeStamp(int(st_y),int(st_m),int(st_d),0,0,0).AsDouble()
        en_y,en_m,en_d = (fname.split("_")[2]).split("-")
        stop_dt  = ROOT.TTimeStamp(int(en_y),int(en_m),int(en_d),0,0,0).AsDouble()

        # check if within the specified time period
        if start_dt > start_time.AsDouble() and stop_dt < stop_time.AsDouble():
            files_to_use.append(join(indir,fname))
            #print fname

        # Print summary and add files to TChains
    print "INFO \t Found %i files in wanted time span from %s to %s" %(len(files_to_use),start_time.AsString("s"),stop_time.AsString("s"))
    for ftu in files_to_use:
        c_trend.Add(ftu)
        c_headr.Add(ftu)
        c_weath.Add(ftu)

    # The time offset used in the POLA data is 1st of January 2007
    t = ROOT.TTimeStamp(2007,1,1,0,0,0)
    da = ROOT.TDatime(2007,1,1,0,0,0)
    ROOT.gStyle.SetTimeOffset(da.Convert());

    # Arrays to store values
    x_rawrate = array( 'd' )
    y_rawrate = array( 'f' )

    x_event = array('d')
    y_pressure = array('d')
    y_temp_in = array('d')
    y_temp_out = array('d')
    y_latitude = array('d')
    y_longitude = array('d')
    # raw rate as function of pressure (p)
    y_rawrate_p = array('d')
    # raw rate as function of tempetarue (t)
    y_rawrate_tin = array('d')  # Inddor
    y_rawrate_tout = array('d') # Outdoor

    # Define some histograms
    h_nev_temp = ROOT.TH1F("h_nev_temp","",60,0,60)
    h_time_temp = ROOT.TH1F("h_time_temp","",60,0,60)
    h_err_temp = ROOT.TH1F("h_err_temp","",60,0,60)

    h_nev_pressure = ROOT.TH1F("h_nev_pressure","",1500,0,1500)
    h_time_pressure = ROOT.TH1F("h_time_pressure","",1500,0,1500)
    h_err_pressure = ROOT.TH1F("h_err_pressure","",1500,0,1500)

    # Get the number of entries
    nentries_headr = c_headr.GetEntries()
    nentries_weath = c_weath.GetEntries()
    nentries_trend = c_trend.GetEntries()

    print "Number of runs   :", nentries_headr
    print "Number of events :", nentries_trend

    # Check that number of events corresponds with the number of runs (run = event)
    if nentries_headr != nentries_weath:
        sys.exit("ERROR \t Number of events in header and weather do not match!. Exiting")


    # Set vb to non-zero value if you want some extra print-out
    vb = 0
    storeStartTime_run = True
    tot_numevents = 0.0
    tot_duration = 0
    tot_latitude = 0.0
    tot_longitude = 0.0

    # Loop over runs
    for i in range(nentries_headr):
#        if i%100 == 0: print "Run %i/%i" %(i,nentries_headr)

        # Get entry i from TTrees (important!)
        c_headr.GetEntry(i)
        c_weath.GetEntry(i)

        # Exclude events with few events
        if c_headr.NumEvents < 12000: continue

        # Get duration, start and end time for run (in seconds from 1/1/2007)
        start = long(c_headr.RunStart)
        stop  = long(c_headr.RunStop)

        run_duration  = getDuration(start,stop)

        # Some runs seem to have some bad start and stop values,
        # exclude those
        if run_duration > 1000000:
            print "skipping", fname
            continue
        # Skip short runs (less than 10 minutes)
        if (run_duration/60.) < 10: continue

        # store start of time bin
        if storeStartTime_run:
            start_new_block = start
            storeStartTime_run = False

        if vb > 0:
            print "Run duration = %d" %(run_duration)
            print "Start: %d %d" %(t.GetDate(True,start),t.GetTime(True,start))
            print "Stop : %d %d" %(t.GetDate(True,stop),t.GetTime(True,stop))


        # Accumulate number of events and time
        tot_duration  += run_duration
        tot_numevents += c_headr.NumEvents

        # If reached the wanted time span for adding a new point (default: 12 hours)
        if tot_duration >= 12.*60.*60.:
            # Fill x-values with the time (in seconds) in the middle of the range
            x_rawrate.append((start_new_block+(tot_duration/2.)))
            # Fill the raw rate in y-value
            y_rawrate.append(tot_numevents/tot_duration)
            y_rawrate_p.append(c_weath.Pressure)
            y_rawrate_tin.append(c_weath.IndoorTemperature)
            y_rawrate_tout.append(c_weath.OutdoorTemperature)
            # if k==0:
            #     print "tot_latitude = ", tot_latitude
            #     print "tot_duration = ", tot_duration
            #     print "c_headr latitude = ", c_headr.latitude
            #     print "y_latitude = ", tot_latitude/float(tot_duration), "\n"
            y_latitude.append(c_headr.latitude)
            y_longitude.append(c_headr.longitude)
            # reset variables
            tot_duration = 0.0
            tot_numevents = 0.0

            storeStartTime_run = True

        # Fill histograms for every event
        h_nev_temp.Fill(c_weath.IndoorTemperature,c_headr.NumEvents)    # indoor temperature recorded for each event
        h_time_temp.Fill(c_weath.IndoorTemperature,run_duration)        # indoor temperature recorded over time

        # This histogram is needed to get the proper statistical uncertainty
        # (i.e. same histograms as above, but storing the number of entries, N)
        # Error is given by sqrt(N)
        h_err_temp.Fill(c_weath.IndoorTemperature)

    # Loop over all the events in this run
    tot_duration_ev = 0
    tot_intemp_ev = 0.0
    tot_outtemp_ev = 0.0
    tot_pressure_ev = 0.0
    nev_in_block = 0
    storeStartTime_ev = True
    for j in range(nentries_trend):
        c_trend.GetEntry(j)

#        if j%1000 == 0: print "Event %i/%i" %(j,nentries_trend)

        # Start and stop time of events
        ev_start = long(c_trend.BinStart)
        ev_stop  = long(c_trend.BinEnd)

        # If start of a new time interval
        if storeStartTime_ev:
            ev_start_block = ev_start
            storeStartTime_ev = False

        # Get length of run (in seconds)
        ev_duration = getDuration(ev_start,ev_stop)

        # Some runs seem to have invalid start/stop times
        # Exclude those which are outside the time range defined at the beginning of the program
        # (i.e. the time of the PolarquEEEst experiment)
        if ev_stop > stop_time.AsDouble():
            continue
        # Skip short runs
        if ev_duration > 61: continue

        tot_duration_ev  += ev_duration
        tot_intemp_ev    += c_trend.IndoorTemperature
        tot_outtemp_ev   += c_trend.OutdoorTemperature

        tot_pressure_ev  += c_trend.Pressure
        nev_in_block += 1

        # If total duration is larger then the chosen time interval (default: 10 min)
        # Store the average values within the chosen time period (default: 10 min)
        if tot_duration_ev >= 10.*60.:
            if vb > 0:
                print "Tot duration for data period ", (tot_duration_ev/(60.))

            # Fill x-axis (time in middle of interval), i.e. total time duration for an event (?)
            x_event.append((ev_start_block+(tot_duration_ev/2.)))

            # Fill y-values (temperature, pressure)
            y_pressure.append(tot_pressure_ev/nev_in_block)
            #average indoor temperature recorded per event in a block
            y_temp_in.append(tot_intemp_ev/nev_in_block)
            #average outdoor temperature recorded per event in a block
            y_temp_out.append(tot_outtemp_ev/nev_in_block)
#            print tot_outtemp_ev/nev_in_block
            #avverage latitude recorded per event in a block
            #y_latitude.append(tot_lat_ev/nev_in_block)
            #avverage longitude recorded per event in a block
            #y_longitude.append(tot_lon_ev/nev_in_block)


            # Reset counters
            storeStartTime_ev = True
            tot_duration_ev = 0.0
            tot_intemp_ev = 0.0
            tot_outtemp_ev = 0.0
            tot_pressure_ev = 0.0
            tot_latitude = 0.0
            tot_longitude = 0.0
            nev_in_block = 0

    x_rawrate_all.append(x_rawrate)
    y_rawrate_all.append(y_rawrate)

    y_rawrate_p_all.append(y_rawrate_p)
    y_rawrate_tin_all.append(y_rawrate_tin)
    y_rawrate_tout_all.append(y_rawrate_tout)

    x_event_all.append(x_event)
    y_pressure_all.append(y_pressure)
    y_temp_in_all.append(y_temp_in)
    y_temp_out_all.append(y_temp_out)

    # # Get mean values and standard deviations of histogram values
    # # Indoor temperature
    # print "Mean indoor temperature "+detector+" = ", statistics.mean(y_temp_in)
    # print "sample std indoor temperature "+detector+" = ", statistics.stdev(y_temp_in), "\n"
    #
    #
    # # Outdoor temperature
    # print "Mean outdoor temperature "+detector+" = ", statistics.mean(y_temp_out)
    # print "sample std outdoor temperature "+detector+" = ", statistics.stdev(y_temp_out), "\n"
    #
    # # Pressure
    # print "Mean pressure "+detector+" = ", statistics.mean(y_pressure)
    # print "sample std pressure "+detector+" = ", statistics.stdev(y_pressure)
    #
    # print "len(y_pressure) = ", len(y_pressure)

    # Reset arrays and data
    c_trend.Reset()
    c_headr.Reset()
    c_weath.Reset()

    # Plot longitude and latitude of POLA-01:

    if k==0:
        """
        C_lat = ROOT.TCanvas("c_lat", "Latitude as a function of time "+sorted(directories)[k], 1200, 600)
        g_latitude = ROOT.TGraph(len(x_rawrate), x_rawrate, y_latitude)
        # The default marker in ROOT is horrible so let's change
        g_latitude.SetMarkerStyle(20)
        g_latitude.SetMarkerColor(4)
        # Get x-axis in date and time instead of seconds since 2007
        g_latitude.GetXaxis().SetTimeDisplay(1)
        # Format the axis (default is date and time split in two lines)
        g_latitude.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
        g_latitude.GetXaxis().SetLabelOffset(0.03)
        g_latitude.GetXaxis().SetTitle("Time")
        g_latitude.GetXaxis().SetTitleSize(0.35)
        g_latitude.GetYaxis().SetTitle("Latitude")
        #g_RawRate.GetYaxis().SetLabelSize(0.35)
        g_latitude.SetTitle("Latitude as a function of time")
        # Draw with axis and points
        g_latitude.Draw("AP")
        C_lat.Draw()
        C_lat.Print("../data/plots/%s/Latitude_%s.pdf" %(sorted(directories)[k], sorted(directories)[k]))

        C_lon = ROOT.TCanvas("c_lon", "Longitude as a function of time "+sorted(directories)[k], 1200, 600)
        g_longitude = ROOT.TGraph(len(x_rawrate), x_rawrate, y_longitude)
        # The default marker in ROOT is horrible so let's change
        g_longitude.SetMarkerStyle(20)
        g_longitude.SetMarkerColor(2)
        # Get x-axis in date and time instead of seconds since 2007
        g_longitude.GetXaxis().SetTimeDisplay(1)
        # Format the axis (default is date and time split in two lines)
        g_longitude.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
        g_longitude.GetXaxis().SetLabelOffset(0.03)
        g_longitude.GetXaxis().SetTitle("Time")
        g_longitude.GetXaxis().SetTitleSize(0.35)
        g_longitude.GetYaxis().SetTitle("Longitude")
        #g_RawRate.GetYaxis().SetLabelSize(0.35)
        g_longitude.SetTitle("Longitude as a function of time")
        # Draw with axis and points
        g_longitude.Draw("AP")
        C_lon.Draw()
        C_lon.Print("../data/plots/%s/Longitude_%s.pdf" %(sorted(directories)[k], sorted(directories)[k]))

        C_voyage = ROOT.TCanvas("c_voyage", "Voyage of "+sorted(directories)[k], 1200, 600)
        g_voyage = ROOT.TGraph(len(y_longitude), y_longitude, y_latitude)
        # The default marker in ROOT is horrible so let's change
        g_voyage.SetMarkerStyle(20)
        g_voyage.SetMarkerColor(9)
        g_voyage.GetXaxis().SetLabelOffset(0.03)
        g_voyage.GetXaxis().SetTitle("Time")
        g_voyage.GetXaxis().SetTitleSize(0.35)
        g_voyage.GetYaxis().SetTitle("Voyage")
        #g_RawRate.GetYaxis().SetLabelSize(0.35)
        g_voyage.SetTitle("Longitude as a function of time")
        # Draw with axis and points
        g_voyage.Draw("AP")
        C_voyage.Draw()
        C_voyage.Print("../data/plots/%s/Voyage_%s.pdf" %(sorted(directories)[k], sorted(directories)[k]))
        """

    """
    C = ROOT.TCanvas("c", "Raw rate as a function of time", 1200, 600)
    g_RawRate = ROOT.TGraph(len(x_rawrate),x_rawrate,y_rawrate)
    # The default marker in ROOT is horrible so let's change
    g_RawRate.SetMarkerStyle(20)
    g_RawRate.SetMarkerColor(4)
    # Get x-axis in date and time instead of seconds since 2007
    g_RawRate.GetXaxis().SetTimeDisplay(1)
    # Format the axis (default is date and time split in two lines)
    g_RawRate.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
    g_RawRate.GetXaxis().SetLabelOffset(0.03)
    g_RawRate.GetXaxis().SetTitle("Time")
    g_RawRate.GetXaxis().SetTitleSize(0.35)
    g_RawRate.GetYaxis().SetTitle("Raw rate")
    #g_RawRate.GetYaxis().SetLabelSize(0.35)
    g_RawRate.SetTitle("Raw rate as a function of time")
    # Draw with axis and points
    g_RawRate.Draw("AP")
    #C.Draw()
    C.Print("../data/plots/%s/RawRate_%s.pdf" %(sorted(directories[k]), sorted(directories[k])))

    C1 = ROOT.TCanvas("c1", "Pressure as a function of time", 1200, 600)
    g_pressure = ROOT.TGraph(len(x_event),x_event,y_pressure)
    # The default marker in ROOT is horrible so let's change
    g_pressure.SetMarkerStyle(20)
    g_RawRate.SetMarkerColor(4)
    # Get x-axis in date and time instead of seconds since 2007
    g_pressure.GetXaxis().SetTimeDisplay(1)
    # Format the axis (default is date and time split in two lines)
    g_pressure.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
    g_pressure.GetXaxis().SetLabelOffset(0.03)
    g_pressure.GetXaxis().SetTitle("Time")
    g_pressure.GetYaxis().SetTitle("Pressure")
    g_pressure.SetTitle("Pressure as a function of time")
    # Draw with axis and points
    g_pressure.Draw("AP")
    #C1.Draw()
    C1.Print("../data/plots/%s/pressure_%s.pdf" %(sorted(directories[k]), sorted(directories[k])))


    C3 = ROOT.TCanvas("c3", "c3", 1200, 600)
    h_div_temp = h_nev_temp.Clone("h_div_temp")
    h_div_temp.SetMarkerStyle(20)
    # Set the errors correctly
    for i in range(1,h_div_temp.GetNbinsX()+1):
        h_nev_temp.SetBinError(i,h_err_temp.GetBinError(i))
        h_time_temp.SetBinError(i,h_err_temp.GetBinError(i))
        h_div_temp.Divide(h_nev_temp,h_time_temp)
        h_div_temp.SetTitle("")
        h_div_temp.Draw("e")
        C3.Draw()
        #C3.Print("../data/plots/%s/IndoorTemp_EventsOverTime_%s.pdf" %(sorted(directories[i]), sorted(directories[i])))


    C4 = ROOT.TCanvas("c4", "Indoor temperature over time", 1200, 600)
    g_indoortemp = ROOT.TGraph(len(x_event), x_event, y_temp_in)
    # The default marker in ROOT is horrible so let's change
    g_indoortemp.SetMarkerStyle(20)
    # Get x-axis in date and time instead of seconds since 2007
    g_indoortemp.GetXaxis().SetTimeDisplay(1)
    # Format the axis (default is date and time split in two lines)
    g_indoortemp.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
    g_indoortemp.GetXaxis().SetLabelOffset(0.03)
    g_indoortemp.GetXaxis().SetTitle("Time")
    g_indoortemp.GetYaxis().SetTitle("Temperature")
    g_indoortemp.SetTitle("Indoor temperature over time")
    # Draw with axis and points
    g_indoortemp.Draw("AP")
    C4.Draw()
    #C4.Print("../data/plots/%s/IndoorTemp_%s.pdf" %(sorted(directories[i]), sorted(directories[i])))

    C5 = ROOT.TCanvas("c5", "Outdoor temperature over time", 1200, 600)
    g_outdoortemp = ROOT.TGraph(len(x_event), x_event, y_temp_out)
    # The default marker in ROOT is horrible so let's change
    g_outdoortemp.SetMarkerStyle(20)
    # Get x-axis in date and time instead of seconds since 2007
    g_outdoortemp.GetXaxis().SetTimeDisplay(1)
    # Format the axis (default is date and time split in two lines)
    g_outdoortemp.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
    g_outdoortemp.GetXaxis().SetLabelOffset(0.03)
    g_outdoortemp.GetXaxis().SetTitle("Time")
    g_outdoortemp.GetYaxis().SetTitle("Temperature")
    g_outdoortemp.SetTitle("Outdoor temperature over time")
    # Draw with axis and points
    g_outdoortemp.Draw("AP")
    C5.Draw()
    #C5.Print("../data/plots/%s/OutdoorTemp_%s.pdf" %(sorted(directories[i]), sorted(directories[i])))
    """

directories = sorted(directories)

"""
C1 = ROOT.TCanvas("c1","Raw rate",1200, 600)
mg1 = ROOT.TMultiGraph("mg1","Raw rate")

markercolors = [1, 2, 4]
for i in xrange(len(directories)):
    g_RawRate = ROOT.TGraph( len(x_rawrate_all[i]), x_rawrate_all[i], y_rawrate_all[i] )
    g_RawRate.SetName("gr1")
    g_RawRate.SetTitle(directories[i])
    g_RawRate.SetMarkerStyle(3)
    g_RawRate.SetMarkerColor(markercolors[i])
    g_RawRate.SetDrawOption("AP")
    g_RawRate.SetFillStyle(0)
    g_RawRate.GetYaxis().SetLimits(0, 40)

    mg1.Add(g_RawRate)

mg1.GetXaxis().SetTimeDisplay(1)
# Format the axis (default is date and time split in two lines)
mg1.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
mg1.GetXaxis().SetLabelOffset(0.03)
#mg1.GetXaxis().SetTitle("Time")
#g_RawRate.GetXaxis().SetLabelSize(0.35)
mg1.GetYaxis().SetTitle("Rate [Hz]")
mg1.GetYaxis().SetLimits(0, 40)
mg1.Draw("AP")
C1.Draw()
C1.BuildLegend()
C1.Print("../data/plots/RawRate_all.pdf")

C2 = ROOT.TCanvas("c2","Pressure",1200, 600)
mg2 = ROOT.TMultiGraph("mg2","Pressure as a function of time")

markercolors = [1, 2, 4]
for i in xrange(len(directories)):
    g_pressure = ROOT.TGraph( len(x_event_all[i]), x_event_all[i], y_pressure_all[i] )
    g_pressure.SetName("gr2")
    g_pressure.SetTitle(directories[i])
    g_pressure.SetMarkerStyle(7)
    g_pressure.SetMarkerColor(markercolors[i])
    g_pressure.SetDrawOption("AP")
    g_pressure.SetFillStyle(0)
    g_pressure.GetYaxis().SetLimits(0, 40)
    mg2.Add(g_pressure)

mg2.GetXaxis().SetTimeDisplay(1)
# Format the axis (default is date and time split in two lines)
mg2.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
mg2.GetXaxis().SetLabelOffset(0.03)
mg2.GetYaxis().SetTitle("Pressure")
#mg2.GetYaxis().SetTitleOffset(0.2);
#mg2.GetYaxis().SetTitleSize(0.7)
print "Y axis title size = ", mg2.GetYaxis().GetTitleSize()
mg2.GetYaxis().SetLimits(0, 40)
mg2.Draw("AP")
C2.Draw()
C2.BuildLegend()
C2.Print("../data/plots/pressure_all.pdf")
"""

#**************************************************************
#   Plot raw rate as function of temperature and pressure
#**************************************************************

"""
C3 = ROOT.TCanvas("c3","Raw rate; pressure",1200, 600)
mg3 = ROOT.TMultiGraph("mg3","Raw rate as a function of pressure")

markercolors = [1, 2, 4]
for i in xrange(len(directories)):
    g_RawRate_p = ROOT.TGraph( len(y_rawrate_p_all[i]), y_rawrate_p_all[i], y_rawrate_all[i] )
    g_RawRate_p.SetName("gr3")
    g_RawRate_p.SetTitle(directories[i])
    g_RawRate_p.SetMarkerStyle(20)
    g_RawRate_p.SetMarkerColor(markercolors[i])
    g_RawRate_p.SetDrawOption("AP")
    g_RawRate_p.SetFillStyle(0)
    g_RawRate_p.GetYaxis().SetLimits(0, 40)
    mg3.Add(g_RawRate_p)

mg3.GetXaxis().SetTitle("Pressure")
#mg3.GetYaxis().SetTitleOffset(0.2);
#mg3.GetYaxis().SetTitleSize(0.7)
mg3.GetYaxis().SetLimits(0, 40)
mg3.GetYaxis().SetTitle("Rate [Hz]")
mg3.Draw("AP")
C3.Draw()
C3.BuildLegend()
C3.Print("../data/plots/RawRate_pressure_all.pdf")


C4 = ROOT.TCanvas("c4","Raw rate; temp indoor",1200, 600)
mg4 = ROOT.TMultiGraph("mg4","Raw rate as a function of indoor temperature")

markercolors = [1, 2, 4]
for i in xrange(len(directories)):
    g_RawRate_tin = ROOT.TGraph( len(y_rawrate_tin_all[i]), y_rawrate_tin_all[i], y_rawrate_all[i] )
    g_RawRate_tin.SetName("gr4")
    g_RawRate_tin.SetTitle(directories[i])
    g_RawRate_tin.SetMarkerStyle(20)
    g_RawRate_tin.SetMarkerColor(markercolors[i])
    g_RawRate_tin.SetDrawOption("AP")
    g_RawRate_tin.SetFillStyle(0)
    g_RawRate_tin.GetYaxis().SetLimits(0, 40)
    mg4.Add(g_RawRate_tin)

mg4.GetXaxis().SetTitle("Temperature Indoor")
#mg4.GetYaxis().SetTitleOffset(0.2);
#mg4.GetYaxis().SetTitleSize(0.7)
mg4.GetYaxis().SetLimits(0, 40)
mg4.GetYaxis().SetTitle("Rate [Hz]")
mg4.Draw("AP")
C4.Draw()
C4.BuildLegend()
C4.Print("../data/plots/RawRate_IndoorTemp_all.pdf")


C5 = ROOT.TCanvas("c5","Raw rate; temp outdoor",1200, 600)
mg5 = ROOT.TMultiGraph("mg5","Raw rate as a function of outdoor temperature")

markercolors = [1, 2, 4]
legend = ROOT.TLegend(0.1,0.7,0.48,0.9);
#legend.SetHeader("The Legend Title","C"); # option "C" allows to center the header
for i in xrange(len(directories)):
    g_RawRate_tout = ROOT.TGraph( len(y_rawrate_tout_all[i]), y_rawrate_tout_all[i], y_rawrate_all[i] )
    g_RawRate_tout.SetName("gr5")
    g_RawRate_tout.SetTitle(directories[i])
    g_RawRate_tout.SetMarkerStyle(20)
    g_RawRate_tout.SetMarkerColor(markercolors[i])
    g_RawRate_tout.SetDrawOption("AP")
    g_RawRate_tout.SetFillStyle(0)
    g_RawRate_tout.GetYaxis().SetLimits(0, 40)
    mg5.Add(g_RawRate_tout)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.9)
    legend.AddEntry(mg5,directories[i],"p")

mg5.GetXaxis().SetTitle("Temperature Outdoor")
#C5.SetTitleSize(0.3, "t")
#mg5.GetYaxis().SetTitleOffset(0.2);
#mg5.GetYaxis().SetTitleSize(0.7)
mg5.GetYaxis().SetLimits(0, 40)
mg5.GetYaxis().SetTitle("Rate [Hz]")
mg5.Draw("AP")
C5.Draw()
legend.Draw();
#C5.BuildLegend()
#C5.Print("../data/plots/RawRate_OutdoorTemp_all.pdf")
"""

"""
# Variations in raw rate vs latitude
y_rawrate_var01 = array('d')
y_rawrate_var02 = array('d')
variations = array('d')
y_rawrate_all_true = array('d')

# Excluding possible outliers
for i in xrange(len(y_rawrate_all[0])):
    if y_rawrate_all[0][i] >= 30 and y_rawrate_all[0][i] <= 40:
        y_rawrate_all_true.append(y_rawrate_all[0][i])

print "length including outliers: ", len(y_rawrate_all[0])
print "length excluding outliers: ", len(y_rawrate_all_true)

mu1 = statistics.mean(y_rawrate_all_true) #mean POLA-01
mu2 = statistics.mean(y_rawrate_all[1]) #mean POLA-02, i.e. reference point
for i in xrange(len(y_rawrate_all_true)):
    y_rawrate_var01.append(y_rawrate_all_true[i]/mu1)
    y_rawrate_var02.append(y_rawrate_all[1][i]/mu2)
    variations.append(y_rawrate_var01[i]/y_rawrate_var02[i])

print "Mean (POLA-01/<POLA-01>)/(POLA-02/<POLA-02>): ", statistics.mean(variations)
print "Standard devitation (POLA-01/<POLA-01>)/(POLA-02/<POLA-02>): ", statistics.stdev(variations)


# Find gradients of raw rate as funciton of pressure and temperature
# excluding possible outliers:

y_rawrate_all_true = []
y_rawrate_p_all_true = []
y_rawrate_tout_all_true = []
for i in xrange(len(y_rawrate_all)):
    y_rawrate_p_true = array('d')
    y_rawrate_tout_true = array('d')
    y_rawrate_true = array('d')
    print "length pressure including outliers 0"+str(i)+": ", len(y_rawrate_p_all[i])
    print "length temperature including outliers 0"+str(i)+": ", len(y_rawrate_tout_all[i])

    for j in xrange(len(y_rawrate_all[i])):
        if i==0 and y_rawrate_all[i][j] >= 30 and y_rawrate_all[i][j] <= 40:
            y_rawrate_true.append(y_rawrate_all[i][j])
            y_rawrate_p_true.append(y_rawrate_p_all[i][j])
            y_rawrate_tout_true.append(y_rawrate_tout_all[i][j])

        if i==1 and y_rawrate_all[i][j] >= 30 and y_rawrate_all[i][j] <= 40:
            y_rawrate_true.append(y_rawrate_all[i][j])
            y_rawrate_p_true.append(y_rawrate_p_all[i][j])
            y_rawrate_tout_true.append(y_rawrate_tout_all[i][j])

        if i==2 and y_rawrate_all[i][j] >= 25 and y_rawrate_all[i][j] <= 30:
            y_rawrate_true.append(y_rawrate_all[i][j])
            y_rawrate_p_true.append(y_rawrate_p_all[i][j])
            y_rawrate_tout_true.append(y_rawrate_tout_all[i][j])

    print "length pressure excluding outliers 0"+str(i)+": ", len(y_rawrate_p_true)
    print "length temperature excluding outliers 0"+str(i)+": ", len(y_rawrate_tout_true)

    y_rawrate_all_true.append(y_rawrate_true)
    y_rawrate_p_all_true.append(y_rawrate_p_true)
    y_rawrate_tout_all_true.append(y_rawrate_tout_true)

print len(y_rawrate_all_true)
for i in xrange(len(y_rawrate_all_true)):
    rr_pres = np.polyfit(y_rawrate_p_all_true[i], y_rawrate_all_true[i], 1)
    rr_temp = np.polyfit(y_rawrate_tout_all_true[i], y_rawrate_all_true[i], 1)

    print "pressure linear regression 0"+str(i)+" = ", rr_pres
    print "outdoor temperature linear regression 0"+str(i)+" = ", rr_temp
"""


#**************************************************************
#   Applying barometric correction coefficient
#**************************************************************

# Calculate barometric correction coefficient gamma(p)
def gamma_fit(p, p_ref, rawrate):
    x = array('d')
    y = array('d')
    for j in xrange(len(p)):
        x.append(p[j] - p_ref)
        y.append(np.log(rawrate[j]))
    rr_pres = np.polyfit(x, y, 1)

    beta = rr_pres[0]
    alpha = rr_pres[1]

    return alpha, beta

alpha = []
beta = []
gamma_p_all = []
y_correctedrate_all = []
p_ref = [1011.85, 1008.53, 985.87]

C6 = ROOT.TCanvas("c6","Corrected rate [mb]",1200, 600)
mg6 = ROOT.TMultiGraph("mg6","Corrected rate")

markercolors = [1, 2, 4]
#legend = ROOT.TLegend(0.1,0.7,0.48,0.9)

plain  = ROOT.TStyle("Plain","Plain Style (no colors/fill areas)")
plain.SetCanvasBorderMode(0)
plain.SetPadBorderMode(0)
plain.SetPadColor(0)
plain.SetCanvasColor(0)
plain.SetTitleColor(0)
plain.SetStatColor(0)
plain.SetTextSize(0.3)
for i in xrange(len(directories)):
    gamma_p = array('d')
    gamma_t = array('d')
    y_correctedrate = array('d')

    # obtain constant (alpha) and gradient (beta) for correction coefficient
    a, b = gamma_fit(y_rawrate_p_all[i], p_ref[i], y_rawrate_all[i])
    alpha.append(a)
    beta.append(b)

    # apply the barometic correction coefficients to the raw rate
    for j in xrange(len(y_rawrate_all[i])):
        gamma_p.append(np.exp(0 + beta[i]*(y_rawrate_p_all[i][j] - p_ref[i])))
        y_correctedrate.append(y_rawrate_all[i][j]/gamma_p[j])

    gamma_p_all.append(gamma_p)
    y_correctedrate_all.append(y_correctedrate)

    # plot corrected rate vs. time
    g_correctedRate = ROOT.TGraph( len(y_correctedrate_all[i]), x_rawrate_all[i], y_correctedrate_all[i] )
    g_correctedRate.SetName("gr6")
    g_correctedRate.SetTitle(directories[i])
#    g_correctedRate.SetTitleSize(0.3, "t")
    g_correctedRate.SetMarkerStyle(20)
    g_correctedRate.SetMarkerColor(markercolors[i])
    g_correctedRate.SetDrawOption("AP")
    g_correctedRate.SetFillStyle(0)
    g_correctedRate.GetYaxis().SetLimits(0, 40)
    g_correctedRate.GetYaxis().SetTitleSize(0.9)
#    g_correctedRate.Draw()
    mg6.Add(g_correctedRate)

mg6.GetXaxis().SetTimeDisplay(1)
# Format the axis (default is date and time split in two lines)
mg6.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M}");
mg6.GetXaxis().SetTitle("Time")
mg6.GetXaxis().SetLabelOffset(0.03)
mg6.GetYaxis().SetLimits(0, 40)
mg6.GetYaxis().SetTitle("Rate [Hz]")
mg6.Draw("AP")
C6.UseCurrentStyle()
C6.Draw()
#legend.Draw();
C6.BuildLegend()
#C6.Print("../data/plots/CorrectedRate_all.pdf")
