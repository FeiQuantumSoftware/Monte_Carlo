"""Provide the primary functions."""


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())


def zen(with_attribution=True):
    quote = """Beautiful is better than ugly.
    Explicit is better than implicit.
    Simple is better than complex.
    Complex is better than complicated.
    Flat is better than nested.
    Sparse is better than dense.
    Readability counts.
    Special cases aren't special enough to break the rules.
    Although practicality beats purity.
    Errors should never pass silently.
    Unless explicitly silenced.
    In the face of ambiguity, refuse the temptation to guess.
    There should be one-- and preferably only one --obvious way to do it.
    Although that way may not be obvious at first unless you're Dutch.
    Now is better than never.
    Although never is often better than *right* now.
    If the implementation is hard to explain, it's a bad idea.
    If the implementation is easy to explain, it may be a good idea.
    Namespaces are one honking great idea -- let's do more of those!"""

    if with_attribution:
        quote += "\n\tTim Peters"

    return quote



#monte carlo package functions:
import matplotlib.pyplot as plt
import numpy as np
import random
from math import exp


#def spin configuration class
class SpinConfig():
    
    def __init__(self,N_length=10):
        """N = total spin_number"""
        self.N_length = N_length
        self.iMax = 2**self.N_length
        self.spinlist=[]



    #spinlist initialization
    def init_input_decimal(self, decimal_input):
  
        if decimal_input < self.iMax:
            binary_list = [int(element) for element in list(bin(decimal_input)[2:])]
            while len(binary_list) < self.N_length:
                binary_list = [0] + binary_list
            self.spinlist = binary_list
        else:
            print("Error:input decimal is bigger than possible spinconfig can present")
        
        return self.spinlist

    
    
    def init_rand_spinlist(self):
        
        return self.init_input_decimal(random.randint(0,self.iMax-1))



    #spinlist manipulation: 
    def random_flip(self):
        
        random_site=random.randint(0,self.N_length-1)
        if self.spinlist[random_site] == 0:
            self.spinlist[random_site]=1
        else:
            self.spinlist[random_site]=0

        return self.spinlist


    
    def input_str(self, str_input):

        binary_list2=list()
        for element in str_input:
            if element =="+":
                binary_list2.append(1)
            elif element == "-":
                binary_list2.append(0)
            else:
                print("Error:input should be '+' or '-' string'")
                break
        
        return binary_list2
    


    # spinlist properties: 
    def magnetization(self):
        
        self.magnet = 2 * self.spinlist.count(1)-self.N_length
        
        return self.magnet
        
#         self.magnet = 0   
#         for eachspin in self.spinlist:
#             if eachspin == 1:
#                 self.magnet += 1
#             elif eachspin == 0:
#                 self.magnet += -1
#             else:
#                 pass
                
#         return self.magnet



    def hamiltonian(self,J=-2,u=1.1):
        
        self.u = u
        self.J = J
        self.energy = 0.0
        
        self.energy = self.u * (2* self.spinlist.count(1) - self.N_length)

        newList = self.spinlist[1:]
        newList.append(self.spinlist[0])
        for spinx, spiny in zip(self.spinlist, newList):
            if spinx==spiny:
                self.energy += -self.J * 1
            elif spinx!=spiny:
                self.energy += -self.J * (-1)
            else:
                print("Error: hamiltonian of the spinlist has error")
                
        return self.energy



    #Observable
    def observable_theory(self,T=10,J=-2,u=1.1):
        
        self.J = J
        self.u = u
        
        Zsum = 0.0
        E_theory = 0.0
        EE_theory = 0.0
        m_theory = 0.0
        mm_theory = 0.0
    
        for i_list in range(self.iMax):
            self.spinlist = self.init_input_decimal(i_list)
            self.magnetization()
            self.hamiltonian(self.J,self.u)
            Zi = exp(-self.energy/T)
            
            Zsum +=Zi 
            E_theory += Zi * self.energy
            EE_theory += Zi * self.energy**2
            m_theory += Zi * self.magnet
            mm_theory += Zi * self.magnet**2
        
        """Normalize over Zsum"""
        self.E_theory = E_theory/Zsum
        EE_theory = EE_theory/Zsum
        self.m_theory = m_theory/Zsum
        mm_theory = mm_theory/Zsum
        #get capacity and magnetic susceptibility
        self.C_theory = (EE_theory - self.E_theory**2)/(T*T)
        self.ms_theory = (mm_theory - self.m_theory**2)/(T)
        
        return self.E_theory, self.m_theory, self.C_theory, self.ms_theory


    
    def observable_metropolis_sampling(self,T=10,sample_size_M=10000,u=1.1,J=-2):
        
                
        self.J = J
        self.u = u
        
        #initalize
        self.init_rand_spinlist()
        
        spin_metro_sample = self.spinlist
        E_metro_sample = self.hamiltonian()
        m_metro_sample = self.magnetization()

        E_metro_sum = E_metro_sample
        m_metro_sum = m_metro_sample
        EE_metro_sum = E_metro_sample**2
        mm_metro_sum = m_metro_sample**2
        
        #generate the rest M-1 random_flipped spinlist samples:
        j=1
        while j < sample_size_M:
            self.random_flip()
            dE = self.hamiltonian() - E_metro_sample
            
            #decision
            if dE < 0 or random.random() < exp(-dE/T):
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
                
        #average
        self.E_metropolis = E_metro_sum/sample_size_M
        self.m_metropolis = m_metro_sum/sample_size_M
        self.C_metropolis = (EE_metro_sum/sample_size_M - self.E_metropolis**2)/(T*T)
        self.ms_metropolis = (mm_metro_sum/sample_size_M - self.m_metropolis**2)/(T)
        
        return self.E_metropolis, self.m_metropolis, self.C_metropolis, self.ms_metropolis


        
        