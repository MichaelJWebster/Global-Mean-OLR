'''
OlrGlobal provides an interface to the NOAA interpolated once daily OLR
measurements. The data is described here:

    http://www.esrl.noaa.gov/psd/data/gridded/data.interp_OLR.html

Created on 21/04/2014

@author: michaelwebster
'''
import sys
from netCDF4 import *
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import os

class OlrGlobal(object):
    '''
    OlrGlobal produces various data dumps from the NOAA inerpolated once daily
    OLR netcdf file.
    
    The data is indexed by the three dimensions, longitude, latitude and time.
    
    - Longitude takes values in [0.0, 357.5] in increments of 2.5 degrees. The
    unit of Longitude is degrees East.
    - Latitude takes values in [90, -90] in increments of 2.5 degrees. The
    units of Latitude is degrees North.
    - Time takes values in [1528872, 1875144] in increments of the number
    of hours in a month. The units of Time are hours since 12:00AM 01/01/1800.
    
    
    '''
    dset = None
    lats = None
    longs = None
    times = None
    values = None

    offset = 327.65
    #
    # Mapping from years to a dictionary  containing a list of months, and the
    # year index.
    #
    yearsToMonths = None
    ml = "monthList"
    yi = "yearIndex"
    mi = "monthIndex"
    
    #
    # Mapping from months to the year they are in.
    #
    monthsToYear  = None
    
    #
    # Mapping between hours from 01/01/1800 at 12:00 AM to year values.
    #
    hoursToDate = None
    
    #
    # Mapping between latitudes and indexes.
    #
    latsToIndexes = None
 
    #
    # This is our transformed map from monthly recordings of OLR, to yearly
    # averaged OLR for each (lat, long) coordinate.
    #
    olrByYear = None
    
    #
    # An array of OLR averages over all months for each year y, over all
    # longitudes for latitude l, indexed as:
    #
    # OLR value for latitude x in year y is olrYearAvgsByLatidude[x][y]
    #
    olrYearAvgsByLatitude = None

    #
    # Store intermediate results. If we can find them, we don't have to do
    # recalculate.
    #
    processedDataFname = "./olrAverages.npz"

    def __init__(self, filename="/Users/michaelwebster/Desktop/climate/OLRMeasurements/olr.mon.mean.nc"):
    #def __init__(self, filename="/Users/michaelwebster/Desktop/climate/OLRMeasurements/olr.mon.ltm.nc"):
        '''
        Constructor
        '''
        self.dset = Dataset(filename, 'r')
        self.latDim = self.dset.variables['lat']
        self.longDim = self.dset.variables['lon']
        self.timeDim = self.dset.variables['time']
        self.values = self.dset.variables['olr']
        
        #
        # Setup the latitude to index mappings.
        #
        self.latsToIndexes = dict()
        lat_length = self.latDim.shape[0]
        idx = 0
        while idx < lat_length:
            self.latsToIndexes[self.latDim[idx]] = idx
            idx = idx+1 
        self.setupYearMappings()
       
        if not self.loadArrays():
            self.setupYearlyData()
            self.getLatitudeYearlyAvgs()
            self.saveArrays()
    
    def getLats(self):
        return self.latDim[:]
    
    def getLongs(self):
        return self.longDim[:]
    
    def getTimes(self):
        return self.timeDim[:]
    
    def getValues(self):
        return self.values[:]

    def setupYearMappings(self):
        units = self.timeDim.units
        times_dim = self.getTimes()
        num_months = len(times_dim)
        self.yearsToMonths = dict()
        self.monthsToYear  = dict() 
        self.hoursToDate   = dict() 
        for i in range(0, num_months):
            monthNumber = times_dim[i]
            m_date = num2date(monthNumber, units)
            year = m_date.year
            if not (year in self.yearsToMonths.keys()):
                self.yearsToMonths[year] = dict()
                self.yearsToMonths[year][self.yi] = year - 1974
                self.yearsToMonths[year][self.ml] = list()
                self.yearsToMonths[year][self.mi] = i
            self.yearsToMonths[year][self.ml].append(monthNumber)
            self.monthsToYear[monthNumber] = year
            self.hoursToDate[monthNumber] = m_date
    
    def setupYearlyData(self):
        numYears = len(self.yearsToMonths.keys())
        
        # create olrByYear as a 3 dimensional array with dimensions of
        # latitude, longitude, and year
        current_times, current_lats, current_longs = self.values.shape
        self.olrByYear = np.zeros((numYears, current_lats, current_longs), np.float32)

        for year in self.yearsToMonths.keys():
            year_dict = self.yearsToMonths[year]
            year_index = year_dict[self.yi]
            monthList = year_dict[self.ml]
            monthIndex = year_dict[self.mi]
            print("Getting Year Averages for Year = %d" % year)
            self.getYearAvgs(self.olrByYear[year_index,:], monthList, monthIndex)

    def getYearAvgs(self, oBY, month_list, monthIndex):
        lats = self.getLats()
        longs = self.getLongs()
        for lat in range(0, oBY.shape[0], 1):
        #for lat in range(37, 38):
            for long in range(0, oBY.shape[1], 1):
                monthData = self.values[:,lat,long]
                oBY[lat,long] = np.average(monthData[monthIndex:(monthIndex + len(month_list))])
    
    def getLatitudeYearlyAvgs(self):
        if self.olrByYear is None:
            self.setupYearlyData()

        # setup the new olrYearAvgsByLatitude array
        numYears = len(self.yearsToMonths.keys())
        numLats = self.latDim.shape[0]
        self.olrYearAvgsByLatitude = np.zeros((numLats, numYears), dtype=np.float64)
        for lat in range(0, numLats, 1):
        #for lat in range(37, 38):
            for yr in range(0, numYears):
                self.olrYearAvgsByLatitude[lat,yr] = self.getLatAvg(lat, yr)

    def plotOlrForYearsAndLatitudes(self, year_list, lat_list):
        years_sorted = sorted(year_list)
        if self.olrYearAvgsByLatitude is None:
            self.getLatitudeYearlyAvgs()
        olrs = list()
        starting_idx = self.latsToIndexes[lat_list[0]]
        ending_idx = self.latsToIndexes[lat_list[1]]
        if starting_idx > ending_idx:
            ei = starting_idx
            starting_idx = ending_idx
            ending_idx = ei 
        print("Starting index is %d, ending index is %d" % (starting_idx, ending_idx))
        for year in years_sorted:
            print("Calculating OLR for year %d" % year)
            yr_idx = self.yearsToMonths[year][self.yi]
            new_avg = 0.0
            idx = starting_idx    
            while idx <= ending_idx:
                print("Calculatinging avg for (year,lat) = (%d, %d)" % (year, self.getLats()[idx]))
                new_avg = self.olrYearAvgsByLatitude[idx, yr_idx] + new_avg
                print("New Avg is %d" % new_avg)
                idx = idx + 1
            new_avg = new_avg/float(ending_idx - starting_idx + 1)
            print ("ending_idx - starting_idx + 1 = %d" % (ending_idx - starting_idx + 1))
            print("Final Avg is %d" % new_avg)
            olrs.append(new_avg)
        plt.plot(years_sorted, olrs)
        plt.plot(years_sorted, olrs, 'g^')
        plt.axis([years_sorted[0], years_sorted[len(years_sorted)-1], new_avg - 10, new_avg + 5])
        plt.title(r'Average Yearly OLR for Latitudes %d$^{\circ} \rightarrow$ %d$^{\circ}$ Degrees North' \
                  % (self.getLats()[starting_idx],self.getLats()[ending_idx]))
        plt.xlabel("Year")
        plt.ylabel(r'Average OLR $Wm^{-2}$')
        m = np.polyfit(years_sorted, olrs, 1)
        yfits = np.polyval(m, years_sorted)
        plt.plot(years_sorted, yfits)
        plt.show()
        return
    
    def getLatAvg(self, lat, yr):
        byLong = self.olrByYear[yr,lat,:]
        latAvg = np.average(byLong)
        print("Lat average for year = %d, latitude = %f is %f" % (self.yearsToMonths.keys()[yr], lat, latAvg))
        return latAvg
    
    def printLatYearAvgs(self):
        if self.olrYearAvgsByLatitude is None:
            self.getLatitudeYearlyAvgs()
        for lat in range(0, self.latDim.shape[0],1):
        #for lat in range(37, 38):
            for year in self.yearsToMonths.keys():
                yi = self.yearsToMonths[year][self.yi]
                avg = self.olrYearAvgsByLatitude[lat, yi]
                print("The average for latitude(%d) for year %d is: %f" % (self.getLats()[lat], year, avg))
        
    def getYearlyData(self):
        return self.olrByYear

    def getHoursToDate(self):
        return self.hoursToDate
    
    def getYearsToMonths(self):
        return self.yearsToMonths
    
    def getMonthsToYear(self):
        return self.monthsToYear

    def saveArrays(self):
        if (self.olrByYear is not None) and (self.olrYearAvgsByLatitude is not None):
            np.savez(self.processedDataFname, olrByYear=self.olrByYear, olrYearAvgsByLatitude=self.olrYearAvgsByLatitude)
        
    def loadArrays(self):
        if not os.path.exists(self.processedDataFname):
            return False
        else:
            data = np.load(self.processedDataFname)
            self.olrByYear = data['olrByYear']
            self.olrYearAvgsByLatitude = data['olrYearAvgsByLatitude']
            return True
        
if __name__ == '__main__':
    dset = OlrGlobal()
    for i in dset.yearsToMonths.keys():
        print "Year is: %d" % i
        print "YearIndex is: %d" % dset.yearsToMonths[i][dset.yi]
        print "Months are: %s" % str(dset.yearsToMonths[i][dset.ml])
    
    for i in sorted(dset.monthsToYear.keys()):
        print("Month(%d) is in Year(%d)" % (i, dset.monthsToYear[i]))
        print("Hours(%d) = date(%s)" % (i, dset.hoursToDate[i]))
    #sys.exit(0)
    ytom = dset.getYearsToMonths()
    htod = dset.getHoursToDate()
    for year in ytom:
        monthList = ytom[year][dset.ml]
        print "<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>"
        print("Year %d has dates:" % year)
        for m in monthList:
            print(htod[m])
    #dset.printLatYearAvgs()
    dset.plotOlrForYearsAndLatitudes(range(1980,2014), [-90, 90])
    
