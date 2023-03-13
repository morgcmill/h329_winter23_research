#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from scipy.optimize import curve_fit
import numpy as np



class Galaxy_Scales:
    """ A way to calculate the scale length and scale height of a galaxy given certain pybody Profile arrays
    
    Attributes:
    radius: radial array of the stars' positions in the galaxy in kpc
    mass in radial direction: the mass of each star corresponding with the radial position above
    z-direction: array of the positions in kpc of each star in the z-direction
    mass in z direction: every stellar mass corresponding to the star in the z positions above
    """
    
    def __init__(self, radius, mass_rad, z_dir, mass_z_dir):
        """Return a galactic scale object"""
        self.radius = radius
        self.mass_rad = mass_rad
        self.z_dir = z_dir
        self.mass_z_dir = mass_z_dir
        
    def scale_length(self):
        """Return scale length of given object"""
        
        def linear(my_x, my_m, my_b):
            y = (my_m*my_x) + my_b
            return y
        
        ## We want to define a condition to isolate desired radial section and make arrays specifically for 5-20 kpc
        condition = np.logical_and(self.radius > 5, self.radius < 20)
        
        ## Here are the arrays for 5-20 kpc. We will use this in the fitting
        new_radius = np.extract(condition, self.radius)
        new_mass_rad = np.extract(condition, self.mass_rad)
        
        ## I am using curve_fit function and defining the fit parameters as fitpars_lin and the error as error_lin
        fitpars_lin, error_lin = curve_fit(linear, new_radius, np.log(new_mass_rad))
            
            ## Need to np.log the fit for the mass because those are the graph's scaling that we are fitting to
            ## fitpars_lin returns our fit variables as (m, b) or (slope, additive constant)
        
        ## scale length (r_0) = -1 / slope
        r_0 = -1 / fitpars_lin[0]
        
        return r_0
    
    def scale_height_parab(self):
        """Return scale length of given object - METHOD: using parabolic fit
        
        Parabolic fit is the typical best fit since we took the log of a gaussian when we did our z vs stellar mass plots.
        So the log of a Gaussian should give us a parabolic graph to fit to. With the data it does not always fit perfect."""
        
        ## z_var is the z-axis array
        ## my_p2 represents the term in front of z^2 in the polynomial, 
        ## and my_const is the additive constant
        
        def quad_cen(z_var, my_p2, my_const):
            y = (my_p2 *(z_var**2)) + my_const
            return y
        
        ## fitting to the z and mass in the z direction to the parabolic function
        fitpars_quad, error_quad = curve_fit(quad_cen, self.z_dir, np.log(self.mass_z_dir)) 
            
            ## Need to np.log the fit for the mass because those are the graph's scaling that we are fitting to 
        
        ## Using p2 variable (defined above) to get scale height
        z_0_parab = np.sqrt(-1/fitpars_quad[0])
        
        return z_0_parab
    
    def scale_height_sech2(self):
        """Return scale length of given object - METHOD: using sech^2 fit
        
        The sech^2 fit is another option I wanted to provide because, for h329 at least, the curve typically fits better
        and we can still obtain the scale height approximation."""
        
        ## my_z is the z direction position array
        ## my_z0 is the fitted variable for the curve
        ## my_amp is the additive constant
        
        def sech2(my_z, my_z0, my_amp):
            y = ((np.cosh(my_z*my_z0))**(-2))*my_amp
            return y
        
        ## fitting the z and mass in the z direction to the sech^2 function
        fitpars_sech2, error_sech2 = curve_fit(sech2, self.z_dir, self.mass_z_dir) ## we don't use np.log in this one due to def
        
        ## Using z0 variable to calculate the scale height. Using absolute value since symmetric over the y-axis to avoid neg.
        z_0_sech2 = np.abs(-1/(fitpars_sech2[0]))
        
        return z_0_sech2
    

