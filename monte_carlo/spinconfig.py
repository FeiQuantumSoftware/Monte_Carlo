"""
SpinConfig class and methods for initialization, manipulation and analyzing property
"""

import numpy as np


# def spin configuration class
class SpinConfig:
    def __init__(self, N_length=10):
        """Create a class of 1-d Ising model, with length of the spinlist as N_length.

        Parameters
        ----------
        N_length: integer , optional
            Length of the list.

        Returns
        -------
        SpinConfig : class
            A class of spinlist with length N_length. The total possible spin configurations number is iMax = 2**N_length.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> myspin.N_length
        8
        >>> myspin.iMax
        256
        """
        self.N_length = N_length
        self.iMax = 2**self.N_length
        self.spinlist = []

    # spinlist initialization
    def init_input_decimal(self, decimal_input):
        """Initialize spin configuration with a decimal input.

        Parameters
        ----------
        decimal_Input : integar
            The decimal value of a binary list.

        Returns
        -------
        self.spinlist : list
            A spin list represented by a binary list. '0' represents: spin down. '1' represents: spin up.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> myspin.init_input_decimal(10)
        [0, 0, 0, 0, 1, 0, 1, 0]
        """
        if decimal_input < self.iMax:
            binary_list = [int(element) for element in list(bin(decimal_input)[2:])]
            while len(binary_list) < self.N_length:
                binary_list = [0] + binary_list
            self.spinlist = binary_list
        else:
            raise ValueError(f"input decimal({decimal_input}) exceeds the biggest spinconfig 2**N={self.iMax}. ")

        return self.spinlist

    def init_rand_spinlist(self):
        """Initialize the spin configuration with random integer in range[0,iMax).

        Parameters
        ----------

        Returns
        -------
        self.spinlist : list
            A binary list of N_length.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> myspin.init_rand_spinlist()
        [0, 1, 1, 0, 1, 0, 1, 0]
        """
        return self.init_input_decimal(np.random.randint(0, self.iMax - 1))

    # spinlist manipulation:
    def random_flip(self):
        """Random flip spin on a random site for a given spinlist.

        Parameters
        ----------

        Returns
        -------
        self.spinlist : list
            A binary spinlist with the spin on a random site flipped.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> myspin.init_rand_spinlist()
        [0, 1, 1, 0, 1, 0, 1, 0]
        >>> myspin.random_flip()
        [0, 1, 0, 0, 1, 0, 1, 0]
        """
        random_site = np.random.randint(0, self.N_length - 1)
        if self.spinlist[random_site] == 0:
            self.spinlist[random_site] = 1
        else:
            self.spinlist[random_site] = 0

        return self.spinlist

    def input_str(self, str_input):
        """Translate the string of the input string in '+' and '-' into a binary list.

        Parameters
        ----------
        str_input: string
            A spin list represented in string with '+': spin up, and '-': spin down.

        Returns
        -------
        spinlist2 : list
            A binary spin list: '0' represents spin down, and '1' represents spin up.

        Examples
        --------
        >>> myspin = SpinConfig()
        >>> myspin.input_str("++-+---+--+")
        [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1]
        """
        binary_list2 = list()
        for element in str_input:
            if element == "+":
                binary_list2.append(1)
            elif element == "-":
                binary_list2.append(0)
            else:
                raise TypeError("input_str: input should be a string of - and +.")
                break

        return binary_list2

    # spinlist properties:
    def magnetization(self):
        """Calculate the magnetization of the spinlist.

        Parameters
        ----------

        Returns
        -------
        magnet: integer
            magnetization of the spinlist

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> mySpin.init_input_decimal(10)
        >>> mySpin.magnetization()
        -4
        """
        self.magnet = 2 * self.spinlist.count(1) - self.N_length

        return self.magnet

    def hamiltonian(self, J=-2, u=1.1):
        """Calculate the energy of the given spinlist.

        Parameters
        ----------
        J: float, optional
            Coupling parameter, default J=-2 .
        u: float, optional
            External field strength, default u=1.1 .

        Returns
        -------
        energy : float
            Total energy from external field and the coupling between the nearest neighbors.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> mySpin.init_input_decimal(10)
        >>> mySpin.hamiltonian()
        -4.4
        """

        self.u = u
        self.J = J
        self.energy = 0.0
        # energy from external field H = Sum_i(u*S[i])
        self.energy = self.u * (2 * self.spinlist.count(1) - self.N_length)

        # energy from coupling between nearest spins.
        newList = self.spinlist[1:]
        newList.append(self.spinlist[0])
        for spinx, spiny in zip(self.spinlist, newList):
            if spinx == spiny:
                self.energy += -self.J * 1
            else:
                self.energy += -self.J * (-1)

        return self.energy

    # Observable
    def observable_theory(self, T=10, J=-2, u=1.1):
        """Calculate oberservables of 1-d Ising model with N_length theoretically
         under temperature T, wtih external field parameter u and coupling parameter J.

        Parameters
        ----------
        T : float, optional
            Temperature
        J: float, optional
            Coupling parameter, default J=-2 .
        u: float, optional
            External field strength, default u=1.1 .

        Returns
        -------
        E, m, C, ms : set
            Expectation of energy, average magnetism, heat capacibility, magnetic susceptbility.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> mySpin.observable_theory()
        (-3.6772068591549063,
         -0.5894627003462397,
         0.32593415709340545,
         0.5351140013397603)
        """

        self.J = J
        self.u = u

        # Sum up obserable of all possible spin configurations
        Zsum = 0.0
        E_theory = 0.0
        EE_theory = 0.0
        m_theory = 0.0
        mm_theory = 0.0

        for i_list in range(self.iMax):
            self.spinlist = self.init_input_decimal(i_list)
            self.magnetization()
            self.hamiltonian(self.J, self.u)
            Zi = np.exp(-self.energy / T)

            Zsum += Zi
            E_theory += Zi * self.energy
            EE_theory += Zi * self.energy**2
            m_theory += Zi * self.magnet
            mm_theory += Zi * self.magnet**2

        # Normalize over Zsum
        self.E_theory = E_theory / Zsum
        EE_theory = EE_theory / Zsum
        self.m_theory = m_theory / Zsum
        mm_theory = mm_theory / Zsum

        # get capacity and magnetic susceptibility
        self.C_theory = (EE_theory - self.E_theory**2) / (T * T)
        self.ms_theory = (mm_theory - self.m_theory**2) / (T)

        return self.E_theory, self.m_theory, self.C_theory, self.ms_theory

    def observable_metropolis_sampling(self, T=10, sample_size_M=10000, u=1.1, J=-2):
        """
        Simulated averaged energy, magnetization, heat Capacity and magnetic susceptbility
         of 1-d Monte Carlo of sample_size_M under temperature T.

        Parameters
        ----------
        T : float, optional
            Temperature
        sample_size_M : int
            sample size of the metroplis sample, default sample_size_M = 10000
        u: float, optional
            External field strength, default u=1.1
        J: float, optional
            Coupling parameter, default J=-2

        Returns
        -------
        E, m, C, ms : set
            Average energy, average magnetism, heat capacibility, magnetic susceptbility.

        Examples
        --------
        >>> myspin = SpinConfig(8)
        >>> mySpin.observable_metropolis_sampling()
        (-1.2713199999999993, -0.2052, 0.36647321457601023, 0.6763492959999999)
        """

        self.J = J
        self.u = u

        # initalize the 1st sample
        self.init_rand_spinlist()

        spin_metro_sample = self.spinlist
        E_metro_sample = self.hamiltonian()
        m_metro_sample = self.magnetization()

        E_metro_sum = E_metro_sample
        m_metro_sum = m_metro_sample
        EE_metro_sum = E_metro_sample**2
        mm_metro_sum = m_metro_sample**2

        # generate the rest M-1 metroplis samples by random flip one spin with decision:
        j = 1
        while j < sample_size_M:
            self.random_flip()
            dE = self.hamiltonian() - E_metro_sample

            # decision
            if dE < 0 or np.random.random() < np.exp(-dE / T):
                j += 1
                spin_metro_sample = self.spinlist
                E_metro_sample = self.energy
                self.magnetization()
                m_metro_sample = self.magnet

                E_metro_sum += E_metro_sample
                EE_metro_sum += E_metro_sample**2
                m_metro_sum += m_metro_sample
                mm_metro_sum += m_metro_sample**2

            else:
                self.spinlist = spin_metro_sample

        # average to get the simulated observables
        self.E_metropolis = E_metro_sum / sample_size_M
        self.m_metropolis = m_metro_sum / sample_size_M
        self.C_metropolis = (EE_metro_sum / sample_size_M - self.E_metropolis**2) / (
            T * T
        )
        self.ms_metropolis = (mm_metro_sum / sample_size_M - self.m_metropolis**2) / (
            T
        )

        return (
            self.E_metropolis,
            self.m_metropolis,
            self.C_metropolis,
            self.ms_metropolis,
        )
