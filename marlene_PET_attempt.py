import numpy as np
from scipy.optimize import root
from scipy.optimize import root_scalar

# -------------------- CONSTANTS -----------------------------------------------------------------
blood_density = 1050  # kg/m^3
blood_spec_heatcapacity = 3600  # J/(kg*K)
skin_emissivity = 0.98
sigma = 5.67e-8  # Stefan-Boltzmann constant (W/(m^2*K^4))
standard_pressure = 1013.25  # Standard atmospheric pressure (kPa)
eta_running = 0.2  # Mechanical efficiency from Höppe (1993) (no unit)
r = 2260e3  # Vaporization heat of water (J/kg)
air_spec_heatcapacity = 1010  # spec heat capacity for air, J/kgK
f_eff = 0.7  # effective surface for radiative exchanges (no unit)
epsilon_p = 0.98  # emission coeff for the human body, Bernard et al (2013)
m = 1 / (0.79 * 10e7)  # permeance coeff of the skin for water vapor (kg / (s·Pa))


# ---------------------------- SIDE FUNCTIONS -----------------------------------------------------------------

def get_saturation_vapor_pressure(temp):
    """Calculates the saturation vapor pressure."""
    SVP = 6.1078 * (10 ** (7.5 * temp / (temp + 237.3)))
    return SVP


def get_RTM(met_rate):
    """Calculates the mass of air expired per second (RTM)."""
    RTM = met_rate * 0.000144 / 1000 * 1.2 # PLACEHOLDER
    return RTM


def get_bloodflow(Tc, Tsk):
    """Find the blood flow from core to skin."""
    blood_flow = (6.3 + 75 * (Tc - 36.6)) / (1 + 0.5 * (34.0 - Tsk))
    blood_flow = 0.002
    return blood_flow


def get_sweatrate(Tsk, Tc, sex, body_area):
    """
    Find the sweat rate in kg/sm^2.
    """
    tbody = 0.1 * Tsk + 0.9 * Tc  # body temperature, FRÅN FREDRIK

    if sex == 0:
        # sweat_rate = 8.47 * 10e-5 * ((0.1 * Tsk + 0.9 * Tc) - 36.6)  # for men, HÖPPE
        sweat_rate = 304.94 * (tbody - 36.6) * body_area / 3600000  # sweating rate kg/s, FREDRIK

    if sex == 1:
        sweat_rate = 0.7 * 8.47 * 10e-5 * ((0.1 * Tsk + 0.9 * Tc) - 36.6)  # for women

    return sweat_rate


# ----------------------- STEP 1: FIND TSK AND TC FOR RUNNING (NAKED) -----------------------------------------------

def running_energy_balance_naked(vars, Ta, Tmrt, RH, v, met_rate, body_mass, height, sex):
    Tsk, Tc = vars  # x[0] = Tsk, x[1] = Tc

    A = 0.203 * body_mass ** 0.425 * height ** 0.725  # Body surface area, Du-Bois formula (m^2)
    vps = get_saturation_vapor_pressure (Ta) # kPa
    VP = (RH * vps / 100) # ambient vapor pressure, kPa
    T_ex = 0.47 * Ta + 21.0 # temp expired air, °C
    SVP_Tsk = get_saturation_vapor_pressure(Tsk)
    SVP_Tex = get_saturation_vapor_pressure(T_ex)
    RTM = get_RTM(met_rate)
    h_c = 1.16*(10.45-v+10*v**(1/2)) # Heat transfer coeff dep on wind speed, Engineers toolbox, W/(K*m^2)
    #h_c = 2.75 + 6.5 * v ** 0.67 # FREDRIK ANVÄNDER DENNA
    h_e = 0.633 * h_c / (standard_pressure * air_spec_heatcapacity) # FREDRIK ANVÄNDER DENNA I E_Sw


    # FYSIOLOGISKA BERÄKNINGAR
    blood_flow = get_bloodflow(Tc, Tsk)
    sweat_rate = get_sweatrate(Tsk, Tc, sex, A)

    # VÄRMETRANSPORT - Simplified for naked body
    Fcs = blood_flow * blood_density * blood_spec_heatcapacity * (Tc - Tsk)  # core → skin

    E_res = RTM * air_spec_heatcapacity * (Ta - T_ex)
    E_rel = RTM * r * (SVP_Tex - VP) / standard_pressure
    E_re = E_res + E_rel
    R = A * f_eff * epsilon_p * sigma * ((Tmrt) ** 4 - (Tsk) ** 4) # Strålning
    C = A * h_c * (Ta - Tsk) # Konvektion
    E_d = m * r * (SVP_Tsk - VP)

    # Begränsad svettavdunstning
    E_sw_phy = sweat_rate * r * A # Fredrik har negativ SW
    E_sw_pot = h_e * A *  r * (0.622 / standard_pressure) * (SVP_Tsk - VP) # Något med enheterna här???
    E_Sw = min(E_sw_phy, E_sw_pot)
    print(f"  E_sw_phy: {E_sw_phy:.2f} W, E_sw_pot: {E_sw_pot:.2f} W")

    M = met_rate
    W = eta_running * M
    H = M - W
    S = 0  # Steady state

    # SYSTEMET: Vi vill att summan av alla flöden = 0 (balans)

    eq1 = Fcs - (R + C + E_d + E_Sw)  # Hudbalans: in (blod) = ut (rad + conv + evap)
    eq2 = H - (Fcs + E_re)  # Kropp: in (metaboliskt) = ut (blod + respiration)

    # --- DEBUGGING PRINTS ---
    print("\n--- Heat Flow Values (Debugging - Naked) ---")
    print(f"  Tsk: {Tsk:.2f} °C, Tc: {Tc:.2f} °C")
    print(f"  Fcs (Core->Skin): {Fcs:.2f} W")
    print(f"  R (Radiation Loss/Gain): {R:.2f} W")
    print(f"  C (Convection Loss/Gain): {C:.2f} W")
    print(f"  E_d (Diffusion Evaporation): {E_d:.2f} W")
    print(f"  E_Sw (Sweat Evaporation): {E_Sw:.2f} W")
    print(f"  E_re (Respiration): {E_re:.2f} W")
    print(f"  H (Internal Heat Production): {H:.2f} W")
    print(f"  Eq1 Residual (Skin Balance): {eq1:.2f} W")
    print(f"  Eq2 Residual (Core Balance): {eq2:.2f} W")
    print("--- End Heat Flow Values - Naked ---")

    return [eq1, eq2]

def solve_Temps_naked(Ta, Tmrt, RH, v, met_rate, body_mass, height, sex):
    initial_guess = np.array([34.0, 39])  # [Tsk, Tc]

    solution = root(running_energy_balance_naked, initial_guess,
                    args=(Ta, Tmrt, RH, v, met_rate, body_mass, height, sex),
                    method='lm')

    if solution.success:
        Tsk, Tc = solution.x
        print(f"Steg 1 löst: Tsk = {Tsk:.2f} °C, Tc = {Tc:.2f} °C (Naked)")
        return Tsk, Tc
    else:
        raise RuntimeError("Optimering misslyckades: " + solution.message)

# ----------------------- STEP 2: FIND PET ----------------------------------------------------------------------

def pet_energy_balance(Ta_ref, Tsk, Tc, body_mass, height, sex): # EJ KORRIGERAD FÖR KLÄDER
    #Icl = 0.9 * 0.155  # ref clothing
    met_rate = 80  # ref MET
    A = 0.203 * body_mass ** 0.425 * height ** 0.725  # m²
    v = 0.1
    RH = 50
    Tmrt = Ta_ref  # standard för PET

    vps = get_saturation_vapor_pressure(Ta_ref)
    VP = RH * vps / 100
    T_ex = 0.47 * Ta_ref + 21.0
    SVP_Tsk = get_saturation_vapor_pressure(Tsk)
    SVP_Tex = get_saturation_vapor_pressure(T_ex)
    RTM = get_RTM(met_rate)
    h_c = 2.75 + 6.5 * v ** 0.67
    #Tcl = (Ta_ref + Tmrt + Tsk) / 3

    blood_flow = get_bloodflow(Tc, Tsk)
    sweat_rate = get_sweatrate(Tsk, Tc, sex, A)

    Fcs = blood_flow * blood_density * blood_spec_heatcapacity * (Tc - Tsk)
    Fsc = (Tsk - Tcl) / Icl
    E_res = RTM * air_spec_heatcapacity * (Ta_ref - T_ex)
    E_rel = RTM * r * (SVP_Tex - VP) / standard_pressure
    E_re = E_res + E_rel
    R = A * fcl_ref * f_eff * epsilon_p * sigma * ((Tmrt + 273.15) ** 4 - (Tsk + 273.15) ** 4)
    C = A * fcl_ref * h_c * (Ta_ref - Tsk)
    E_d = m * r * (SVP_Tsk - VP)
    E_Sw = A * sweat_rate * r

    M = met_rate
    W = 0 # eta_running = 0 for PET reference environment!
    H = M - W
    S = 0

    eq1 = Fcs - (R + C + E_d + E_Sw)
    eq2 = H - (Fcs + E_re)

    return eq1 + eq2  # Summerad obalans i hela kroppen

def calculate_PET(Tsk, Tc, body_mass, height, sex):
    def func(Ta_ref):
        return pet_energy_balance(Ta_ref, Tsk, Tc, body_mass, height, sex)

    sol = root_scalar(func, bracket=[-50, 80], method='brentq')

    if sol.converged:
        PET = sol.root
        print(f"PET: {PET:.2f} °C")
        return PET
    else:
        raise RuntimeError("PET-lösning konvergerade inte.")

if __name__ == "__main__":
    # Steg 1: Verkliga förhållanden (utan kläder)
    Ta = 25
    Tmrt = 40
    RH = 50
    v = 4
    met_rate = 600
    body_mass = 70 # kg
    height = 1.75 # m
    sex = 0

    Tsk_naked, Tc_naked = solve_Temps_naked(Ta, Tmrt, RH, v, met_rate, body_mass, height, sex)
    #PET = calculate_PET(Tsk, Tc, body_mass, height, sex)
