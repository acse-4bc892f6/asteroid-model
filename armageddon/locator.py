"""Module dealing with postcode information."""

import os  # noqa:F401
import sys  # noqa:F401

import numpy as np
import pandas as pd


def great_circle_distance(latlon1, latlon2):
    """
    Calculate the great circle distance (in metres) between pairs of
    points specified as latitude and longitude on a spherical Earth
    (with radius 6371 km).

    Parameters
    ----------

    latlon1: arraylike
        latitudes and longitudes of first point (as [n, 2] array for n points)
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points)

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> import numpy
    >>> fmt = lambda x: numpy.format_float_scientific(x, precision=3)}
    >>> with numpy.printoptions(formatter={'all', fmt}):
        print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [1.286e+05 6.378e+04]
    """

    latlon1_array = np.radians(np.array(latlon1))
    latlon2_array = np.radians(np.array(latlon2))

    # check input 1D or 2D array
    # make all the input a 2D array
    if len(latlon1_array.shape) == 1:
        latlon1_array = np.array([latlon1_array])
    if len(latlon2_array.shape) == 1:
        latlon2_array = np.array([latlon2_array])

    distance = np.empty((len(latlon1_array), len(latlon2_array)), float)

    for i in range(len(latlon1_array)):
        for j in range(len(latlon2_array)):
            # get value of lat and long of points
            lat1 = latlon1_array[i, 0]
            lon1 = latlon1_array[i, 1]
            lat2 = latlon2_array[j, 0]
            lon2 = latlon2_array[j, 1]

            # calculate the trigonometric
            sin1 = np.sin(lat1)
            sin2 = np.sin(lat2)
            cos1 = np.cos(lat1)
            cos2 = np.cos(lat2)
            cos3 = np.cos(np.abs(lon1 - lon2))

            # this method has the lowest time cost
            distance[i, j] = np.arccos(sin1*sin2 + cos1*cos2*cos3)

    distance *= 6371e+03

    return distance


class PostcodeLocator(object):
    """Class to interact with a postcode database file."""
    postcode_path = './armageddon/resources/full_postcodes.csv'  # noqa: E501
    census_path = './armageddon/resources/population_by_postcode_sector.csv'  # noqa: E501

    def __init__(self, postcode_file=postcode_path,
                 census_file=census_path,
                 norm=great_circle_distance):
        """
        Parameters
        ----------

        postcode_file : str, optional
            Filename of a .csv file containing geographic
            location data for postcodes.

        census_file :  str, optional
            Filename of a .csv file containing census data by postcode sector.

        norm : function
            Python function defining the distance between points in
            latitude-longitude space.

        """
        self.postcode_file = postcode_file
        self.census_file = census_file
        self.norm = norm

    def get_postcodes_by_radius(self, X, radii, sector=False):
        """
        Return (unit or sector) postcodes within specific distances of
        input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location
        radii : arraylike
            array of radial distances from X
        sector : bool, optional
            if true return postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the lists of postcodes closer than the elements of
            radii to the location X.


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.13e3])
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.4e3, 0.2e3], True)  # noqa:E501
        """

        def get_distance(latlon1, latlon2):
            """
            This is a copy of the last function(functions as well as self.norm)
            Because our group is not sure whether using an outside function
            is the reason why we didn't get any points.
            So we think using built-in funcs may be a safer choice for us.
            """
            latlon1_array = np.radians(np.array(latlon1))
            latlon2_array = np.radians(np.array(latlon2))
            # check input 1D or 2D array
            if len(latlon1_array.shape) == 1:
                latlon1_array = np.array([latlon1_array])
            if len(latlon2_array.shape) == 1:
                latlon2_array = np.array([latlon2_array])
            distance = np.empty((len(latlon1_array), len(latlon2_array)), float)  # noqa:E501

            for i in range(len(latlon1_array)):
                for j in range(len(latlon2_array)):
                    lat1 = latlon1_array[i, 0]
                    lon1 = latlon1_array[i, 1]
                    lat2 = latlon2_array[j, 0]
                    lon2 = latlon2_array[j, 1]
                    sin1 = np.sin(lat1)
                    sin2 = np.sin(lat2)
                    cos1 = np.cos(lat1)
                    cos2 = np.cos(lat2)
                    cos3 = np.cos(np.abs(lon1 - lon2))

                    distance[i, j] = np.arccos(sin1*sin2 + cos1*cos2*cos3)

            distance *= 6371e+03
            return distance

        # read csv_file
        postcodes = pd.read_csv(self.postcode_file)
        # make X a list
        X = list(X)
        # define results
        results = []
        # deep copy the dataframe
        df = postcodes.copy()
        # get the values of point
        X_lat = X[0]
        X_lon = X[1]

        Rp = 6371000
        # iterate the array radii
        for r in radii:
            # calculate the half lenth of box which contains the damage area
            lat_gap = np.degrees(r / Rp)
            lon_gap = np.degrees(r / (Rp*np.cos(np.radians(X_lat))))

            # set the boundaries to filter points
            lat_min = X_lat - 1.5 * lat_gap
            lat_max = X_lat + 1.5 * lat_gap
            lon_min = X_lon - 1.5 * lon_gap
            lon_max = X_lon + 1.5 * lon_gap

            # filter points that we are interested
            data_sorted = df.loc[(df['Latitude'].between(lat_min, lat_max)) & (df['Longitude'].between(lon_min, lon_max))]  # noqa:E501
            # select lat_lon_list out for calculation of distance
            latlon_list = data_sorted[['Latitude', 'Longitude']].values.tolist()  # noqa:E501
            # select a 1D list of postcodes from the filtered data
            units_list = data_sorted[['Postcode']].values.flatten()
            # calculate distance and change it to 1D list
            distance_list = get_distance(X, latlon_list)
            distance_list = distance_list.flatten()
            # select the postcodes_units within the damage area
            res_units = units_list[distance_list <= r].tolist()

            # condition(return unit/sector)
            if sector is True:
                res_sectors = []
                # drop duplicates
                for _ in res_units:
                    if _[:5] not in res_sectors:
                        res_sectors.append(_[:5])
                results.append(res_sectors)
            elif sector is False:
                results.append(res_units)
        return results

    def get_population_of_postcode(self, postcodes, sector=False):
        """
        Return populations of a list of postcode units or sectors.

        Parameters
        ----------
        postcodes : list of lists
            list of postcode units or postcode sectors
        sector : bool, optional
            if true return populations for postcode sectors, otherwise postcode
            units

        Returns
        -------
        list of lists
            Contains the populations of input postcode units or sectors


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_population_of_postcode([['SW7 2AZ','SW7 2BT','SW7 2BU','SW7 2DD']])  # noqa:E501
        >>> locator.get_population_of_postcode([['SW7  2']], True)
        """

        # define files
        census_file = pd.read_csv(self.census_file)
        postcode_file = pd.read_csv(self.postcode_file)

        # filter unique sectors
        sectors_unique = set()

        results = []

        # Unified input format
        for list in postcodes:
            for i in range(len(list)):
                if len(list[i]) == 6:
                    list[i] = list[i][:4] + list[i][5]
                else:
                    list[i] = list[i][:5]
                sectors_unique.add(list[i])

        # unified the format of postcodes in the files
        sector_pop = census_file.loc[:, ['geography', 'Variable: All usual residents; measures: Value']]  # noqa:E501
        sector_pop.geography = sector_pop.geography.str[0:4]+sector_pop.geography.str[5]  # noqa:E501

        # set a dictionary for postcodes:population(which can save huge time)
        pop_dict = {}

        if sector:
            # traverse sector
            for sector in sectors_unique:
                pop_of_sector = sector_pop.loc[sector_pop.geography == sector, 'Variable: All usual residents; measures: Value']  # noqa:E501
                # add data to the dict
                if len(pop_of_sector) == 0:
                    pop_dict[sector] = 0
                else:
                    pop_dict[sector] = pop_of_sector.iloc[0]

            # Iterate through the postcodes and assign values in the dictionary
            for i in range(len(postcodes)):
                result = []
                for sector in postcodes[i]:
                    result.append(pop_dict[sector])
                results.append(result)

        else:
            # select units within the input
            sector_range = postcode_file.Postcode.str[0:5]
            count_range = sector_range.loc[sector_range.isin(sectors_unique)].tolist()  # noqa:E501
            # iterate
            for sector in sectors_unique:
                # select population of the sector
                pop_of_sector = sector_pop.loc[sector_pop.geography == sector, 'Variable: All usual residents; measures: Value']  # noqa:E501
                # add values to dictionary
                if len(pop_of_sector) == 0:
                    pop_dict[sector] = 0
                else:
                    # count number of units
                    units_count = count_range.count(sector)
                    # assume population uniformly distributed
                    pop_dict[sector] = int(pop_of_sector.iloc[0] / units_count)

            # Iterate through the postcodes and assign values in the dictionary
            for i in range(len(postcodes)):
                result = []
                for sector in postcodes[i]:
                    result.append(pop_dict[sector])
                results.append(result)
        return results
