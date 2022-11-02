import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from armageddon import locator
from armageddon import mapping


def damage_zones(outcome, lat, lon, bearing, pressures):
    """
    Calculate the latitude and longitude of the surface zero location and the
    list of airblast damage radii (m) for a given impact scenario.

    Parameters
    ----------

    outcome: Dict
        the outcome dictionary from an impact scenario
    lat: float
        latitude of the meteoroid entry point (degrees)
    lon: float
        longitude of the meteoroid entry point (degrees)
    bearing: float
        Bearing (azimuth) relative to north of meteoroid trajectory (degrees)
    pressures: float, arraylike
        List of threshold pressures to define airblast damage levels

    Returns
    -------

    blat: float
        latitude of the surface zero point (degrees)
    blon: float
        longitude of the surface zero point (degrees)
    damrad: arraylike, float
        List of distances specifying the blast radii
        for the input damage levels

    Examples
    --------

    # >>> import armageddon
    # >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3,
    #                'burst_distance': 90e3, 'burst_peak_dedz': 1e3,
    #                'outcome': 'Airburst'}
    # >>> armageddon.damage_zones(outcome, 52.79, -2.95, 135,
    # pressures=[1e3, 3.5e3, 27e3, 43e3])
    """
    # calculate the longitude and  the latitude
    Rp = 6.371e6
    r = outcome['burst_distance']
    pressure = np.copy(pressures)

    lat, lon, bearing = np.deg2rad([lat, lon, bearing])

    blat = np.arcsin(
        np.sin(lat) * np.cos(r / Rp) +
        np.cos(lat) * np.sin(r / Rp) * np.cos(bearing))
    blon = lon + np.arctan(
        np.sin(bearing) * np.sin(r / Rp) * np.cos(lat) /
        (np.cos(r / Rp) - np.sin(lat) * np.sin(blat)))
    blat, blon = np.rad2deg([blat, blon])

    def Pre_slover(r,
                   zb=outcome['burst_altitude'],
                   Ek=outcome['burst_energy']
                   ):  # outcome.burst_altitude，outcome.burst_energy
        func = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            func[i] = 3.14e11 * (((r[i]**2) + (zb**2)) /
                                 (Ek**(2 / 3)))**(-1.3) + 1.8e7 * (
                                     ((r[i]**2) +
                                      (zb**2)) / (Ek**(2 / 3)))**(-0.565) - p
        return func

    damrad = fsolve(Pre_slover, [5000.] * len(pressures))
    blat = float(blat)
    blon = float(blon)
    for index in range(len(damrad)):
        if damrad[index] < 1:
            damrad[index] = 0
    damrad = damrad.tolist()

    return blat, blon, damrad


fiducial_means = {
    'radius': 10,
    'angle': 20,
    'strength': 1e6,
    'density': 3000,
    'velocity': 19e3,
    'lat': 51.5,
    'lon': 1.5,
    'bearing': -45.
}
fiducial_stdevs = {
    'radius': 1,
    'angle': 1,
    'strength': 5e5,
    'density': 500,
    'velocity': 1e3,
    'lat': 0.025,
    'lon': 0.025,
    'bearing': 0.5
}


def impact_risk(planet,
                means=fiducial_means,
                stdevs=fiducial_stdevs,
                pressure=27.e3,
                nsamples=100,
                sector=True):
    """
    Perform an uncertainty analysis to calculate the risk for each affected
    UK postcode or postcode sector

    Parameters
    ----------
    planet: armageddon.Planet instance
        The Planet instance from which to solve the atmospheric entry

    means: dict
        A dictionary of mean input values for the uncertainty analysis. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    stdevs: dict
        A dictionary of standard deviations for each input value. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    pressure: float
        The pressure at which to calculate the damage zone for each impact

    nsamples: int
        The number of iterations to perform in the uncertainty analysis

    sector: logical, optional
        If True (default) calculate the risk for postcode sectors, otherwise
        calculate the risk for postcodes

    Returns
    -------
    risk: DataFrame
        A pandas DataFrame with columns for postcode (or postcode sector) and
        the associated risk. These should be called ``postcode`` or ``sector``,
        and ``risk``.
    """
    fiducial_param = {}
    blast_lat_list = []
    blast_lon_list = []
    damage_rad_list = []
    # generate parameters through Gaussian Distribution randomly
    for key in means:
        fiducial_param[key] = np.random.normal(means.get(key), stdevs.get(key),
                                               nsamples)
    postcode_all1 = []
    my_postcodelocator = locator.PostcodeLocator()  # instance
    for i in range(nsamples):
        # Solve the system of differential equations
        # for a given impact scenario
        result = planet.solve_atmospheric_entry(
            radius=fiducial_param['radius'][i],
            angle=fiducial_param['angle'][i],
            strength=fiducial_param['strength'][i],
            density=fiducial_param['density'][i],
            velocity=fiducial_param['velocity'][i])

        result = planet.calculate_energy(
            result)  # Calculate the kinetic energy
        outcome = planet.analyse_outcome(
            result)  # Calculate the impact and airburst stats

        # Calculate surface zero location and the list of airblast damage radii
        blast_lat, blast_lon, damage_rad = damage_zones(
            outcome,
            lat=fiducial_param['lat'][i],
            lon=fiducial_param['lon'][i],
            bearing=fiducial_param['bearing'][i],
            pressures=[pressure])
        # (unit or sector) postcodes within
        # specific distances of input location
        postcodes = my_postcodelocator.get_postcodes_by_radius(
            (blast_lat, blast_lon), radii=damage_rad, sector=sector)
        postcode_all1 += postcodes  # store total postcodes/sector of samples
        blast_lat_list.append(blast_lat)
        blast_lon_list.append(blast_lon)
        damage_rad_list.append(damage_rad)

        print(damage_rad)

    uncertainty_plot = mapping.plot_uncertain(blast_lat_list,
                                              blast_lon_list,
                                              damage_rad_list,
                                              map=None)
    uncertainty_plot.save("uncertainty_plot.html")

    postcode_all2 = []
    # change list of list to list
    for item in postcode_all1:
        postcode_all2 += item
    index = set(postcode_all2)
    for i in set(postcode_all2):
        postcode_all2.count(i)  # calculate the postcode/sector appearing times
    probability = {}  # key is postcode/sector, value is probability
    for i in index:
        probability[i] = postcode_all2.count(i) / nsamples
    # postcode/ postcode sectors that had been hit
    # at least once inside the n times simulation
    probability = {
        key: value
        for key, value in probability.items() if value != 0
    }
    risk = []
    sectors = []
    units = []
    for key in probability:
        if sector:
            # sector is true
            sectors += [key]
        else:
            units += [key]

        # calculate risk
        risk += [
            probability[key] * my_postcodelocator.get_population_of_postcode(
                [[key]], sector=sector)[-1][-1]
        ]

    if sector:
        return pd.DataFrame({'sector': sectors, 'risk': risk})
    else:
        return pd.DataFrame({'postcode': units, 'risk': risk})
