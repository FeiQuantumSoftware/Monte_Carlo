Getting Started
===============

This page details how to get started with monte_carlo. monte_carlo is a package 
which was developed for the MolSSI Bes t Practices workshop.

Installation 
------------
To install monte_carlo, you will need an environment with the following packages:

* Python 3.9
* NumPy 

Once you have these packages installed, you can install molecool in the same environment using 
:: 

    pip install -e .


SpinConfig
----------
Once installed, you can use the package. This example shows how to initialize a spinlist.
::

    import monte_carlo

    example_spinlist = molecool.SpinConfig(8)
    example_spinlist.init_decimal_input(10) = [0, 0, 0, 0, 1, 0, 1, 0]


----------
::
    class SpinConfig:
    def __init__(self, N_length=10):
        
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



Observable Theory
-----------------
::
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

result
------
::
    import matplotlib.pyplot as plt
    import numpy as np

    mySpin = SpinConfig(8)
    Tlist=np.linspace(0.1,10,num=100)

    Elist = list()
    mlist = list()
    Clist = list()
    mslist = list()

    for temperature in Tlist:
        E, m, C, ms = mySpin.observable_theory(temperature)
        Elist.append(E) 
        mlist.append(m) 
        Clist.append(C) 
        mslist.append(ms) 
    
    plt.figure(num = 0, dpi = 120)
    plt.plot(Tlist, Elist,label="<E>")
    plt.plot(Tlist, mlist,label="<m>")
    plt.plot(Tlist, Clist,label="C")
    plt.plot(Tlist, mslist,label="ms")
    plt.legend()
    plt.xlabel("T")


Metroplis Sampling
------------------




