#!/bin/env python

import unyt
import unyt.dimensions

cgs_unit = {
    "U_I" : "A",
    "U_L" : "cm",
    "U_M" : "g",
    "U_T" : "K",
    "U_t" : "s",
}

def units_from_attributes(dset):
    cgs_factor = dset.attrs["Conversion factor to CGS (not including cosmological corrections)"][0]
    U_I = dset.attrs["U_I exponent"][0]
    U_L = dset.attrs["U_L exponent"][0]
    U_M = dset.attrs["U_M exponent"][0]
    U_T = dset.attrs["U_T exponent"][0]
    U_t = dset.attrs["U_t exponent"][0]
    return cgs_factor * (unyt.A**U_I) * (unyt.cm**U_L) * (unyt.g**U_M) * (unyt.K**U_T) * (unyt.s**U_t) 

def cgs_expression(exponents):

    if all([e==0.0 for e in exponents.values()]):
        return "[ - ]"

    expression = ""
    for base_unit, power in exponents.items():
        if power == 0:
            pass
        elif power == 1.0:
            expression += base_unit + " "
        elif power % 1.0 == 0.0:
            expression += ("%s^%d " % (base_unit, power))
        else:
            expression += ("%s^%7.4f " % (base_unit, power))

    expression += "[ "
    for base_unit, power in exponents.items():
        if power == 0:
            pass
        elif power == 1.0:
            expression += cgs_unit[base_unit] + " "
        elif power % 1.0 == 0.0:
            expression += ("%s^%d " % (cgs_unit[base_unit], power))
        else:
            expression += ("%s^%7.4f " % (cgs_unit[base_unit], power))
    expression += "]"
    
    return expression

def write_unit_attributes(dset, unit):

    # Get conversion to CGS
    cgs_factor = unit.get_conversion_factor(unit.get_cgs_equivalent())[0]

    # Find powers of mass, length, time etc
    powers = {}
    coeff, dims = unit.dimensions.as_coeff_mul()
    for dim in dims:
        base, exp = dim.as_base_exp()
        exp = float(exp)
        if base is unyt.dimensions.length:
            powers["U_L"] = exp
        elif base is unyt.dimensions.mass:
            powers["U_M"] = exp
        elif base is unyt.dimensions.current_mks:
            powers["U_I"] = exp
        elif base is unyt.dimensions.temperature:
            powers["U_T"] = exp
        elif base is unyt.dimensions.time:
            powers["U_t"] = exp

    # Add in any missing dimensions
    for dim in "LMITt":
        name = "U_"+dim
        if name not in powers:
            powers[name] = 0.0

    # Find expression for conversion to CGS
    cgs_expr = cgs_expression(powers)

    # Write attributes to the dataset
    dset.attrs["Expression for physical CGS units"] = cgs_expr
    dset.attrs["Conversion factor to CGS (not including cosmological corrections)"] = [cgs_factor,]
    dset.attrs["U_I exponent"] = [powers["U_I"],]
    dset.attrs["U_L exponent"] = [powers["U_L"],]
    dset.attrs["U_M exponent"] = [powers["U_M"],]
    dset.attrs["U_T exponent"] = [powers["U_T"],]
    dset.attrs["U_t exponent"] = [powers["U_t"],]
