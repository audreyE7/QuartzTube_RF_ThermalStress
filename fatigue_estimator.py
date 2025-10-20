"""
Low-cycle thermal fatigue helper using Coffin–Manson and Morrow mean stress correction.
"""
import numpy as np

def coffin_manson_life(delta_eps_p_over_2, eps_f_prime=0.5, c=0.6):
    """Solve 2Nf from Δε_p/2 = εf' * (2Nf)^(-c)"""
    return (eps_f_prime / (delta_eps_p_over_2))**(1.0/c)

def morrow_corrected_strain_range(delta_sigma, E, sigma_m=0.0, sigma_f_prime=900e6, b=0.09):
    """
    Simplified Morrow correction: adjust elastic term with mean stress.
    Return Δε_e and Δε_p/2.
    """
    delta_eps_e = delta_sigma / E
    sigma_f_eff = max(1.0, sigma_f_prime - sigma_m)
    m = 0.25
    delta_eps_p_over_2 = 0.002 * (delta_sigma / sigma_f_eff)**(1.0/m)
    return delta_eps_e, delta_eps_p_over_2
