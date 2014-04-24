#!/usr/bin/env python
'''
OlrGlobal calculates and plots Outgoing Longwave radiation from the the NOAA
interpolated once daily OLR measurements. The data is described here:

    http://www.esrl.noaa.gov/psd/data/gridded/data.interp_OLR.html

The motivation for writing the program is learning about data science, and
learning about climate science. I'll admit to being surprised about the
results. I am interested in investigating the following areas:

1. Can a corresponding natural variation in incoming RF explain these findings,
2. If I add an observed change to incoming RF over the period, will I be able to
show an energy imbalance of the Earth,
3. The averages are a very gross measure, what finer/better measurements can I
do.
4. Is it possible to use multi-threading to speed up the number crunching of
this program, or is there a problem with the algorithm that adds significantly
to the time taken to compute the averages.

This is a first attempt at using numpy and matplotlib.

Created on 21/04/2014

@author: michaelwebster
'''
import sys
import argparse
from netCDF4 import *
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import os
import os.path

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

    The data itself is the interpolated monthly means of satellite measurements
    of OLR over a grid of latitutes and longitudes.

    The Program performs the following on the data:
    - For latitude and longitude combination, take the monthly average data, and
    average it into yearly data.
    - For each latitude and each year, produce the average OLR over all
    longitudes for that data.
    - Plot the average of the averaged latitude/year averages, over some period
    of time, and for a range of latitudes.
    '''
    #
    # A netcdf4 Dataset object.
    #
    dset = None

    #
    # The latitude dimension read from the NOAA netcdf file.
    #
    lats = None

    #
    # The longitude dimension read from the NOAA netcdf file.
    #
    longs = None

    #
    # The time dimension read from the NOAA netcdf file.
    #
    times = None

    #
    # The units read from the time dimension of the netcdf file, and passed
    # as a parameter to the netcdf4.num2date method.
    #
    units = None

    #
    # The 'olr' variable read from the NOAA netcdf file.
    #
    values = None

    #
    # FIXME: Not sure what to do with this number. I need to investigate.
    #
    offset = 327.65

    #
    # Mapping from years to a dictionary  containing a list of months recorded
    # in hours since x where x is a value encoded for the time dimension in the
    # netcdf file. For the olr.mon.mean.nc file, the start time is 01-01-1800.
    #
    yearsToMonths = None

    #
    # List of months in hours since the start point for the given year.
    #
    ml = "monthList"

    #
    # The index of year x in the array created from the monthly data.
    #
    yi = "yearIndex"

    #
    # A list of month indexes for the months in the yearsToMonths[yr][self.ml]
    # list.
    #
    mi = "monthIndex"
    
    #
    # Mapping from months to the year they are in.
    #
    monthsToYear  = None

    #
    # startingYear is the first year in the series. We subtract this value from
    # a year value to get the year index value we use to index our arrays.
    #
    startingYear = 0

    #
    # Mapping between hours from the starting date read from the time dimension
    # netcdf file, to years.
    #
    hoursToDate = None
    
    #
    # Mapping between latitudes and indices in the lat dimension for that
    # latitude.
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
    processedDataFname = None #"./olrAverages.npz"

    def __init__(self, filename):
        '''
        The Constructor:
        - Loads the netcdf dataset from the file with name == filename.
        - Read in the dimensions and variables from the netcdf file.
        - Setup a mapping between latitudes and latitude indexes in the
          'lat' dimension.
        - Set up the mapping from years to lists of months and year indexes
          that we will use.
        - Load arrays for olrByYear and olrYearAvgsByLatitude from a file in
          the current directory if we have formerly saved them.
        - Compute and save the olrByYear and olrYearAvgsByLatitude arrays, if
          we cannot load them from a file.
        '''
        self.dset = Dataset(filename, 'r')
        self.latDim = self.dset.variables['lat']
        self.longDim = self.dset.variables['lon']
        self.timeDim = self.dset.variables['time']
        self.values = self.dset.variables['olr']

        #
        # Setup the units used by the num2Date netcdf routine.
        #
        self.units = self.timeDim.units
        
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

        fileNameBase = os.path.basename(filename)
        self.processedDataFname = fileNameBase + ".npz"
        print("Processed filename is %s\n" % self.processedDataFname)
        if not self.loadArrays():
            self.setupYearlyData()
            self.getLatitudeYearlyAvgs()
            self.saveArrays()
    
    def getLats(self):
        """
        Return an array containing all the latitude values.
        """
        return self.latDim[:]
    
    def getLongs(self):
        """
        Return an array containing all the longitude values.
        """        
        return self.longDim[:]
    
    def getTimes(self):
        """
        Return an array containing all the time values.
        """        
        return self.timeDim[:]
    
    def getValues(self):
        """
        Return the 3-dimensional array containing all the OLR values indexed
        by (time, lat, long).
        """
        return self.values[:]

    def setupYearMappings(self):
        """
        Setup mappings between years and months by:
        
        1. Convert each month number into a date, and store the month number in
        a dictionary indexed by the year of the date.

        2. For each year, store an index that starts at 0 for the earliest year,
        and increases by 1 for each subsequent year - assuming the years form
        a sequence.

        3. Store a month index for each month read, in a separate list to the
        actual month nubers.

        4. Store a mapping between each month number, and the year it is in.

        5. Store a mapping for each month number in hours since a start date,
        to the actual date represented by that month number.
        """
        times_dim = self.getTimes()
        self.startingYear = num2date(times_dim[0], self.units).year        
        num_months = len(times_dim)
        self.yearsToMonths = dict()
        self.monthsToYear  = dict() 
        self.hoursToDate   = dict()
        for i in range(0, num_months):
            monthNumber = times_dim[i]
            m_date = num2date(monthNumber, self.units)
            year = m_date.year
            if not (year in self.yearsToMonths.keys()):
                self.yearsToMonths[year] = dict()
                self.yearsToMonths[year][self.yi] = year - self.startingYear
                self.yearsToMonths[year][self.ml] = list()
                self.yearsToMonths[year][self.mi] = i
            self.yearsToMonths[year][self.ml].append(monthNumber)
            self.monthsToYear[monthNumber] = year
            self.hoursToDate[monthNumber] = m_date
    
    def setupYearlyData(self):
        """
        Create the olrByYear array as a 3 dimensional array indexed by
        (year, lat, long). This is done by storing an average yearly value for
        the olr for each (lat, long) coordinate, by averaging over the monthly
        averages recorded in the NOAA data.

        The yearly array calulated here is stored as the olrByYear instance
        variable, and is saved to disk in a compressed form.
        """
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
        """
        getYearAvgs is passed a 2 dimensional array indexed by latitude and
        longitude indices, that is to be filled in with the yearly averages
        calculated from the monthly averages read from the original NOAA data.
        """
        lats = self.getLats()
        longs = self.getLongs()
        for lat in range(0, oBY.shape[0]):
            for long in range(0, oBY.shape[1], 1):
                #
                # monthData is a 1 dimensional slice of the array containing
                # OLR averages for each month for the given (lat,long)
                # coordinate.
                #
                monthData = self.values[:,lat,long]
                oBY[lat,long] = np.average(monthData[monthIndex:(monthIndex + len(month_list))])
    
    def getLatitudeYearlyAvgs(self):
        """
        This method adds a further step to the averaging over years that we did
        previously, by averaging over every longitude for each latitude.

        What we end up with is a 2 dimensional array indexed by (lat, year)
        that contains the averaged OLR readings for that latitude and year. This
        array is stored as the olrYearAvgsByLatitude instance variable, and is
        saved to disk in a compressed form, along with the olrByYear array.
        """
        if self.olrByYear is None:
            self.setupYearlyData()

        # setup the new olrYearAvgsByLatitude array
        numYears = len(self.yearsToMonths.keys())
        numLats = self.latDim.shape[0]
        self.olrYearAvgsByLatitude = np.zeros((numLats, numYears), dtype=np.float64)
        for lat in range(0, numLats, 1):
            print("Getting Averages for Latitude = %d" % self.latDim[lat])
            for yr in range(0, numYears):
                self.olrYearAvgsByLatitude[lat,yr] = self.getLatAvg(lat, yr)

    def plotOlrForYearsAndLatitudes(self, year_list, lat_list):
        """
        Draw a plot of the averages of OLR for each year in year_list, at each
        latitude in lat_list. Draw a best fit line to display the trend.
        """
        years_sorted = sorted(year_list)
        years_sorted = [int(x) for x in years_sorted]
        if self.olrYearAvgsByLatitude is None:
            self.getLatitudeYearlyAvgs()

        starting_idx = self.latsToIndexes[lat_list[0]]
        ending_idx = self.latsToIndexes[lat_list[1]]
        if starting_idx > ending_idx:
            ei = starting_idx
            starting_idx = ending_idx
            ending_idx = ei

        olrs = list()
        for year in years_sorted:
            yr_idx = self.yearsToMonths[year][self.yi]
            new_avg = 0.0
            idx = starting_idx    
            while idx <= ending_idx:
                new_avg = self.olrYearAvgsByLatitude[idx, yr_idx] + new_avg
                idx = idx + 1
            new_avg = new_avg/float(ending_idx - starting_idx + 1)
            olrs.append(new_avg)

        years_modified = [int(x + 2000) if x < 1973 else int(x) for x in years_sorted]
        plt.plot(years_modified, olrs)
        plt.plot(years_modified, olrs, 'g^')
        plt.axis([years_modified[0]-1, years_modified[len(years_modified) - 1]+1, new_avg - 5, new_avg + 5])
        plt.title(r'Average Yearly OLR for Latitudes %d$^{\circ} \rightarrow$ %d$^{\circ}$ Degrees North' \
                  % (self.getLats()[starting_idx],self.getLats()[ending_idx]))
        plt.xlabel("Year")
        plt.ylabel(r'Average OLR $Wm^{-2}$')
        if len(years_modified) > 1:
            m = np.polyfit(years_modified, olrs, 1)
            yfits = np.polyval(m, years_modified)
            plt.plot(years_modified, yfits)
        plt.show()
        return
    
    def getLatAvg(self, lat, yr):
        """
        Calculate the average OLR over all longitudes for latitude lat, for
        year yr.
        """
        byLong = self.olrByYear[yr,lat,:]
        latAvg = np.average(byLong)
        return latAvg
    
    def printLatYearAvgs(self):
        """
        A debug function to printou out the data in the olrYearAvgsByLatitude
        array.
        """
        if self.olrYearAvgsByLatitude is None:
            self.getLatitudeYearlyAvgs()
        for lat in range(0, self.latDim.shape[0],1):
            for year in self.yearsToMonths.keys():
                yi = self.yearsToMonths[year][self.yi]
                avg = self.olrYearAvgsByLatitude[lat, yi]
                print("The average for latitude(%d) for year %d is: %f" % \
                      (self.getLats()[lat], year, avg))
        
    def getYearlyData(self):
        """
        Return the array containing yearly averages for all years, latitudes
        and longitudes.
        """
        return self.olrByYear

    def getYearlyDataByLatitude(self):
        """
        Return the array containing yearly averages for all years, latitudes
        over all longitudes.
        """        
        return self.olrYearAvgsByLatitude

    def getHoursToDate(self):
        """
        Return the mapping between hour values used in the netcdf file, and
        dates.
        """
        return self.hoursToDate
    
    def getYearsToMonths(self):
        """
        Return the yearsToMonths mapping.
        """
        return self.yearsToMonths
    
    def getMonthsToYear(self):
        """
        Return the months to years mapping.
        """
        return self.monthsToYear

    def saveArrays(self):
        """
        Save the calculated arrays of yearly averages to a file in the current
        directory. The name of the file is the same as the input file processed,
        with the added .npz extension.
        """
        if (self.olrByYear is not None) and                                 \
          (self.olrYearAvgsByLatitude is not None):
            np.savez                                                        \
            (                                                               \
                self.processedDataFname,                                    \
                olrByYear=self.olrByYear,                                   \
                olrYearAvgsByLatitude=self.olrYearAvgsByLatitude            \
            )
        
    def loadArrays(self):
        """
        """
        if not os.path.exists(self.processedDataFname):
            return False
        else:
            data = np.load(self.processedDataFname)
            if data is not None:
                self.olrByYear = data['olrByYear']
                self.olrYearAvgsByLatitude = data['olrYearAvgsByLatitude']
                data.close()
                return True
            else:
                return False
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage = sys.argv[0] + \
       " generates averages by year and by latitude from NOAA OLR data.")
    parser.add_argument                                                     \
    (                                                                       \
        '-s', '--src',                                                      \
        dest='src_file', action='store',                                    \
        default='./olr.mon.ltm.nc',                                         \
        help='The netcdf file containing monthly OLR that we want to process.' \
    )
    parser.add_argument                                                     \
    (                                                                       \
        '-p', '--plot',                                                     \
        action='store_true',                                                \
        help='Plot the data for the apropriate latitudes and years.'        \
    )

    parser.add_argument                                                     \
    (                                                                       \
        '-y', '--years', nargs=2,                                           \
        help="If -p is selected, the range of years to plot. The range includes both endpoints."\
    )
    parser.add_argument                                                     \
    (                                                                       \
        '-l', '--lats', nargs=2,                                            \
        help= "%s%s\n%s" %                                                  \
        (                                                                   \
            "If -p is selected, the range of latitudes to plot.",           \
            "The range is inclusive of the two latitudes given.",           \
            "The units of latitude are degrees North."                      \
        )                                                                   \
    )
    args = parser.parse_args()

    if not os.path.exists(args.src_file):
        print("Error: Cannot find the %s olr file." % args.src_file)
        sys.exit(-1)

    #
    # Setup the datset.
    #    
    dset = OlrGlobal(args.src_file)
    
    ytom = dset.getYearsToMonths()

    if args.plot:
        year_start = None
        year_end = None
        if args.years != None:
            year_start = min(int(args.years[0]), int(args.years[1]))
            year_end = max(int(args.years[0]), int(args.years[1]))
        else:
            year_start = ytom.keys()[0]
            year_end = ytom.keys()[len(ytom.keys()) - 1]

        lats_start = None
        lats_end = None

        if args.lats != None:
            lats_start = min(int(args.lats[0]), int(args.lats[1]))
            lats_end = max(int(args.lats[0]), int(args.lats[1]))
        else:
            lats_start = dset.getLats()[0]
            lats_end = dset.getLats()[len(dset.getLats()) - 1]

        
        dset.plotOlrForYearsAndLatitudes(range(year_start,year_end + 1), [lats_start, lats_end])
    
