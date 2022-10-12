import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time


class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """
    def __init__(self,
                 atmos_func='exponential',
                 atmo_filename='armageddon/resources/AltitudeDensityTable.csv',
                 Cd=1.,
                 Ch=0.1,
                 Q=1e7,
                 Cl=1e-3,
                 alpha=0.3,
                 Rp=6371e3,
                 g=9.81,
                 H=8000.,
                 rho0=1.2):

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0
        self.atmos_filename = atmo_filename

        try:
            # set function to define atmoshperic density
            if atmos_func == 'exponential':
                self.rhoa = lambda z: rho0 * np.exp(-z / H)
            elif atmos_func == 'tabular':
                DF = pd.read_csv(
                    atmo_filename,
                    sep=' ',
                    skiprows=6,
                    names=['Altitude', 'Atmospheric_density', 'H'])
                altitude = DF['Altitude'].to_numpy(dtype=float)
                density = DF['Atmospheric_density'].to_numpy(dtype=float)
                self.rhoa = lambda z: density[(np.abs(altitude - z)).argmin()]
                # for z outside range of altitude in table
                # use density of the closest altitude
                self.rhoa = lambda z: density[(np.abs(altitude - z)).argmin()]
            elif atmos_func == 'constant':
                self.rhoa = lambda x: rho0
            else:
                raise NotImplementedError(
                    "atmos_func must be 'exponential', 'tabular' or 'constant'"
                )
        except NotImplementedError:
            print("atmos_func {} not implemented yet.".format(atmos_func))
            print("Falling back to constant density atmosphere for now")
            self.rhoa = lambda x: rho0

    def solve_atmospheric_entry(self,
                                radius,
                                velocity,
                                density,
                                strength,
                                angle,
                                init_altitude=100e3,
                                dt=0.005,
                                radians=False):
        """
        Solve the system of differential equations for a given impact scenario
        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters
        velocity : float
            The entery speed of the asteroid in meters/second
        density : float
            The density of the asteroid in kg/m^3
        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2
        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians
        init_altitude : float, optional
            Initial altitude in m
        dt : float, optional
            The output timestep, in s
        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the dataframe will have the same units as the
            input
        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
            Includes the following columns:
            'velocity', 'mass', 'angle', 'altitude',
            'distance', 'radius', 'time'
        """

        # Enter your code here to solve the differential equations

        # convert angle to radians for sin and cosine functions
        if not radians:
            angle = angle / 360 * (2 * np.pi)

        def RHS(t, state):
            """
            Obtain RHS of ODEs for time stepping
            Parameters
            ----------
            t : float
                Time step
            state : numpy.ndarray
                Current state of the system
            Returns
            -------
            f : numpy.ndarray
                Array of RHS of ODEs
            """
            # initialise array to store RHS of ODEs for each variable
            f = np.zeros_like(state)

            velocity = state[0]
            mass = state[1]
            angle = state[2]
            altitude = state[3]
            state[4]
            radius = state[5]

            f[0] = -(self.Cd * self.rhoa(altitude) * np.pi * radius**2 *
                     velocity**2) / (2 * mass) + self.g * np.sin(angle)
            f[1] = -(self.Ch * self.rhoa(altitude) * np.pi * radius**2 *
                     velocity**3) / (2 * self.Q)
            f[2] = (self.g * np.cos(angle)) / velocity - \
                   (self.Cl * self.rhoa(altitude) * np.pi * radius ** 2 *
                    velocity) / \
                   (2 * mass) - (velocity * np.cos(angle)) / \
                   (self.Rp + altitude)
            f[3] = -velocity * np.sin(angle)
            f[4] = (velocity * np.cos(angle)) / (1 + altitude / self.Rp)

            # airburst criteria
            if self.rhoa(altitude) * velocity**2 < strength:
                f[5] = 0
            else:
                f[5] = np.sqrt((7 * self.alpha * self.rhoa(altitude)) /
                               (2 * density)) * velocity

            return f

        def RK4(f, state0, dt=0.05, t0=0.):
            """
            Implement Runge-Kutta 4 (RK4) stage method
            Parameters
            ----------
            f : function
                RHS
            state0 : list or numpy.ndarray
                Initial state
            dt : float
                Time step size
            t0 : float
                Initial time step
            Returns
            -------
            state_all : numpy.ndarray
                State of system at each time step
            t_all : numpy.ndarray
                All time steps
            """
            # initialise array to store current state and time
            state = np.array(state0)
            t = np.array(t0)

            # initialise list to store state at each time step
            state_all = [state0]

            # initialise list to store all time steps
            t_all = [t0]

            # set limit to number ot iterations
            num_iteration = 0

            while state[3] > 0 and num_iteration < 1e4:
                k1 = dt * f(t, state)
                k2 = dt * f(t + 0.5 * dt, state + 0.5 * k1)
                k3 = dt * f(t + 0.5 * dt, state + 0.5 * k2)
                k4 = dt * f(t + dt, state + k3)
                state = state + (1. / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
                state_all.append(state)
                t = t + dt
                t_all.append(t)
                num_iteration += 1

            return np.array(state_all[:-1][:]), np.array(t_all[:-1][:])

        # initial state from function input
        state0 = [
            velocity, density * 4. / 3. * np.pi * radius**3, angle,
            init_altitude, 0., radius
        ]

        # use RK4 to solve system
        state_all, time_all = RK4(RHS, state0, dt)

        # obtain variables over time steps
        velocity = state_all[:, 0]
        mass = state_all[:, 1]
        angle = state_all[:, 2]
        if not radians:
            angle = angle / (2 * np.pi) * 360
        altitude = state_all[:, 3]
        distance = state_all[:, 4]
        radius = state_all[:, 5]
        # store variables in dataframe
        return pd.DataFrame({
            'velocity': velocity,
            'mass': mass,
            'angle': angle,
            'altitude': altitude,
            'distance': distance,
            'radius': radius,
            'time': time_all
        })

    def calculate_energy(self, result):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.
        Parameters
        ----------
        result : DataFrame
            A pandas dataframe with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time
        Returns : DataFrame
            Returns the dataframe with additional column ``dedz`` which is the
            kinetic energy lost per unit altitude
        """

        # Replace these lines with your code to add the dedz column to
        # the result DataFrame
        result = result.copy()

        # save mass, altitude and velocity as numpy arrays for operations
        mass = result['mass'].to_numpy(dtype=float)
        altitude = result['altitude'].to_numpy(dtype=float)
        velocity = result['velocity'].to_numpy(dtype=float)

        # initialise array to store change in KE
        # and change in KE per unit altitude
        de = np.empty_like(mass)
        dedz = np.empty_like(mass)
        de[0] = 0.
        dedz[0] = 0.

        KEs = 0.5 * mass * velocity**2
        PEs = mass * self.g * altitude

        de[1:] = np.abs((KEs[:-1] + PEs[:-1] - PEs[1:]) -
                        0.5 * mass[1:] * velocity[1:]**2)

        dedz[1:] = np.abs(de[1:] / ((altitude[:-1] - altitude[1:]) / 1000))

        # convert from Joules to kilo tonnes of TNT
        dedz = dedz / 4.184e12

        # insert dedz as column in data frame
        result.insert(len(result.columns), 'dedz', dedz)

        return result

    def analyse_outcome(self, result):
        """
        Inspect a pre-found solution to calculate the impact and airburst stats
        Parameters
        ----------
        result : DataFrame
            pandas dataframe with velocity, mass, angle, altitude, horizontal
            distance, radius and dedz as a function of time
        Returns
        -------
        outcome : Dict
            dictionary with details of the impact event,
            which should contain
            the key ``outcome``
            (which should contain one of the following strings:
            ``Airburst`` or ``Cratering``)
            , as well as the following keys:
            ``burst_peak_dedz``, ``burst_altitude``, ``burst_distance``, ``burst_energy``  # noqa: E501
        """

        # find maximum dedz value and index
        if len(result['dedz']) > 1:
            dedz_max = result['dedz'].max()
            index = result.loc[(result['dedz'] == dedz_max)].index[-1]
        else:
            dedz_max = result['dedz'].iloc[0]
            index = 0

        # find maximum dedz value and distance
        burst_altitude = result['altitude'].iloc[index]
        burst_distance = result['distance'].iloc[index]

        start_mass = result['mass'].iloc[0]
        start_velocity = result['velocity'].iloc[0]
        start_KE = 0.5 * start_mass * start_velocity**2

        if index < (len(result['dedz'])-1):
            status = 'Airburst'
            # find total KE lost from 100km altitude
            # to burst altitude equivalent to Ek
            burst_mass = result['mass'].iloc[index]
            burst_velocity = result['velocity'].iloc[index]
            burst_KE = 0.5 * burst_mass * burst_velocity**2
            total_KE_lost = abs(burst_KE - start_KE) / 4.184e12
            outcome = {
                'outcome': status,
                'burst_peak_dedz': dedz_max,
                'burst_altitude': burst_altitude,
                'burst_distance': burst_distance,
                'burst_energy': total_KE_lost
            }
        else:
            status = 'Cratering'

            altitude_zero = result.loc[result['altitude'] > 0]
            min_altitude = altitude_zero['altitude'].min()
            min_altitude_ind = altitude_zero.loc[altitude_zero['altitude'] ==
                                                 min_altitude].index[-1]
            min_altitude_mass = result['mass'].iloc[min_altitude_ind]
            min_altitude_velocity = result['velocity'].iloc[min_altitude_ind]

            residual_KE = 0.5 * min_altitude_mass * min_altitude_velocity**2

            total_KE_lost = abs(start_KE - residual_KE)

            # the larger out of residual KE and total KE lost is Ek
            Ek = max([total_KE_lost, residual_KE])
            Ek = Ek / 4.184e12

            outcome = {
                'outcome': status,
                'burst_peak_dedz': dedz_max,
                'burst_altitude': burst_altitude,
                'burst_distance': burst_distance,
                'burst_energy': Ek
            }

        return outcome

    def extension2(self):
        """
        Function to find asteroid parameters from the Chelyabinsk
        meteoroid energy deposition curve
        """
        # find best fitted parameter
        data = pd.read_csv(
            'armageddon/resources/ChelyabinskEnergyAltitude.csv')
        data.rename(columns={
            'Height (km)': 'altitude',
            'Energy Per Unit Length (kt Km^-1)': 'dedz'
        },
                    inplace=True)
        data_peak_dedz = data['dedz'].max()
        strength = 1e6
        radius = 20
        velocity = 1.92e4
        density = 3300
        angle = 18.3

        raduis_all = []
        strength_all = []
        altitude_all = []
        dedz_all = []
        err2_all = []
        raduis_all.append(radius)
        strength_all.append(strength)
        error = 50
        D = 0
        while error > 10:
            A = time.time()
            # reduce error2
            result = self.solve_atmospheric_entry(radius,
                                                  velocity,
                                                  density,
                                                  strength,
                                                  angle,
                                                  init_altitude=43e3,
                                                  dt=0.05,
                                                  radians=False)
            result = self.calculate_energy(result)
            altitude_en = result['altitude'] / 1000
            dedz_en = result['dedz']
            altitude_all.append(altitude_en)
            dedz_all.append(dedz_en)
            outcome = self.analyse_outcome(result)
            error = np.abs(outcome['burst_peak_dedz'] - data_peak_dedz)
            error = error
            err2_all.append(error)
            # err_all.append(error)
            # print('error:')
            # print(error)
            if outcome['burst_peak_dedz'] > data_peak_dedz:
                if error > 15:
                    radius -= 0.6
                else:
                    radius -= 0.4
                    strength -= 200
            else:
                if error > 30:
                    radius += 0.6
                else:
                    radius += 0.4
                    strength += 200
            raduis_all.append(radius)
            strength_all.append(strength)
            B = time.time()
            C = B - A
            D += C
            # print(D)
            if D >= 50:
                break
        min_error = np.min(err2_all)
        min_error_index = err2_all.index(min_error)

        final_radius = raduis_all[min_error_index]
        final_strength = strength_all[min_error_index]
        print('the approximate radius is :')  # final best fitted parameter
        print(final_radius)
        print('the approximate strength is:')  # final best fitted parameter
        print(final_strength)
        final_altitude = altitude_all[min_error_index]
        final_dedz = dedz_all[min_error_index]
        plt.plot(final_altitude,
                 final_dedz,
                 'r-',
                 label='fitted energy deposition curve')
        plt.plot(data['altitude'],
                 data['dedz'],
                 'b-',
                 label='exact energy deposition curve')
        plt.xlabel('altitude(m)')
        plt.ylabel('Energy Per Unit Length (kt Km^-1)')
        plt.legend(loc='best')
        plt.grid(True)
        plt.show()
