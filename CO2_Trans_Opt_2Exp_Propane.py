
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from Dearman.TDN import PSI, TDN, mlog, mt, sort


""" 
SYSTEM DEFINITION:

1- OPTIMIZED PISTON EXPANDER AT 2 PRESSURE LEVEL
2- TWO-STAGE INTERCOOLED VCC
3- EXPANSION AT 2 STAGES
4- PROPANE AS REFRIGERANT 
    
"""

"""________________________________________________NOTICE:___________________________________________________________"""
"""____________________________________________All dimensions are in SI unit_____________________________________________"""

"""_____________________________________________________________________________________________________________________
_____________________________________________INITIAL CONDITION__________________________________________________________
_____________________________________________________________________________________________________________________"""
R = ['R404A', 'R134A', 'n-Propane', 'CO2']
C = ['Nitrogen', 'Air', 'Methane', 'Hydrogen']
HX = ['INCOMP::MEG[0.6]', 'INCOMP::AEG[0.6]']

"""-----------INITIAL GUESS OF PRESSURE AND WORK----------"""
pig = [100e5, 40e5, 2e5]

list = []
ss = []
s_sum = []
m_dot = []
t_boil = []
s_tot = []
eta_2_law = np.zeros(len(C))

for cryo in range(len(C)):

    """____FLUIDS____ : """
    r = R[2]  # Refrigerant in the refrigeration cycle
    c = C[cryo]  # Cryogenic fluid used for cooling (Nitrogen in this case)
    hx = HX[0]

    """   __ Model 1 :  # CoolProp mixture Water/ASHRAE, Ethylene Glycol 60% (volume based)
          __ Model 2 :  # Aake Melinder Properties of Secondary Working fluids published 2010 IIR"""

    print('"""____________________________' + str(C[cryo]) + '_____________________________________"""')
    print()

    """_________ Reference Condition (Ambient used in 2nd law) __________ : """
    t0 = 293.15
    p0 = 101325
    rgas = 8.31446261815324  # Gas constant J⋅K−1⋅mol−1

    """___Temperature [K]___ : """

    dt_approach = 3  # Temp difference at heat exchanger
    dt_approach_air = 6  # Min Temp difference at aircooled heat exchanger
    t_evap = 243.15  # Evaporator temperature (T_Min) - 30 C
    t_reefer = t_evap + dt_approach
    t_amb = t0
    t_mt = 268.15

    t_st = PSI('T', 'P', 1.01323e5, 'Q', 0, c)  # Cryogen storage tank temperature(Ambient pressure and quality of 0)
    t_chx = t_evap        # Cryogen after heat exchange in CHX
    t_exp_ig = t_amb # Piston expander average temperature

    dt_subcool_ig = -5  # DT of sub-cooling in (k) (with minus sign)
    dt_superheat_ig = 5    # DT of superheat in (k)

    dt_chx = t_chx - t_st
    """------------------------------------"""
    """----------- Initial Guess ----------"""
    t_cond = t_amb + dt_approach  # Condenser temperature (T_Max)  + 25 C

    """------------------------------------"""
    """___Pressure [Pa]___  : """

    p_evap = PSI('P', 'Q', 1, 'T', t_evap, r)  # Evaporator pressure (P_Min)

    # if t_cond< PSI('Tcrit', r):
    #
    p_cond = PSI('P', 'Q', 1, 'T', t_cond, r)  # Condenser pressure (P_Max)
        # p_cond = PSI('P', 'Q', 1, 'T', t_cond, r)  # Condenser pressure (P_Max)
    #    #
    # else:
    #
    #     p_mt = PSI('P', 'Q', 1, 'T', t_mt, r)
    #
    #     """-------OBJECTIVE FUNCTION TO BE MINIMIZED------"""
    #
    #     def objective(x):
    #
    #         return - ((TDN(x[0],0,t_mt,0,0,r).h).h-TDN(x[0],0,t_mt,0,0,r).h)/compwork)
    #
    #     """--------CONSTRAINTS THAT SHOULD BE HELD TRUE -----------"""
    #
    #     def cont0(x):  # INEQ  # f(x)>= 0 inequality always means >=
    #         return
    #
    #     """___characteristics___ : """
    #     """----------- Initial Guess ----------"""
    #     zero = [f_ig]
    #     """------RANGE OF VARIABLE BONDS --------"""
    #     bon0 = (f_ig, 100)
    #     bon = [bon0]
    #     """ CONSTRANTS EQ OR INEQ SPESIFIED """
    #     con0 = {'type': 'ineq', 'fun': cont0}  # INEQ  always means >= 0
    #     con = [con0]
    #     """------SOLUTION FOR OPTIMIZATION PROBLEM-----"""
    #     sol = minimize(objective, zero, method='SLSQP', bounds=bon, constraints=con)
    #     print('x[0]: ', 'Frequency of piston expander :      ', sol.x[0])



    p_low = p_evap  # p_low of VCC (pa)
    p_high = p_cond  # p_high of VCC (pa)

    p_amb = p0        # Ambient pressure
    p_st = p_amb * 2  # Cryogen pressure at storage tank
    p_exp_max = 4e6   # Max pressure at the piston expander

    p_lim_exp = 4.0e6  # Mechanical limiting pressure for expader

    p_exp_exh = p_amb * 1.5  # Exhaust pressure at the piston expander
    p_exp_min = p_exp_exh    # Min pressure at the piston expander

    dp_suc = 2e4   # Pressure drop at the suction side of compressor (pa)
    dp_dis = 2e4   # Pressure drop at the discharge side of compressor (pa)
    dp_evap = 7e4  # Pressure drop at the Evaporator in (pa)
    dp_cond = 2e4  # Pressure drop at the condenser in (pa)
    dp_liq = 1e4   # Pressure drop at the liquid line (pa)

    """------------------------------------"""

    """___Number of Data Points___ : """

    num_data = 9  # Number of data points for Refrigerant, 4 Component, 2 dp, 2 dt, 1 Molier

    num_data_C = 6  # Number of data points for Cryogen open cycle

    num_data_hx = 6  # Number of data points for HX closed cycle

    num_data_exp_i = 6  # Number of data points Expander engine per seconds

    num_data_exp_r = 100  # Number of data points Expander engine in each RPM

    """------------------------------------"""

    """___Efficiencies___ : """
    eta_hx = 0.95   # Efficiency of heat exchanger
    eta_shaft = 0.97  # Expander to Compressor
    eta_exp_ig = 0.65  # Efficiency of expander
    eta_isen = 0.65  # Isentropic efficiency of the compressor
    eta_aux = 0.95  # Auxillary work required ratio to total work
    eta_pump = 0.9  # Efficiency of heat exchanger
    c_disch = 0.75   #  Discharge coeficient at cryogen inlet of expander

    """------------------------------------"""

    """----------- Cooling required ----------"""
    cool_req = 15e3  # Cooling capacity required from VCC amd cryogenic HX  (W)

    """------------------------------------"""

    """----------- Initial Guess ----------"""
    ig_mdot_c = (500 * 10 ** (-3) / 0.00123) / 6 / 3600  # Initial guess of N2 M_dot entering piston expander in each revolution
    rpm_ig = 1200              # RPM of the
    q_out_ig = cool_req / 3  # Heat exchanged in condensor gas cooler befrore exchanging heat with HX fluid

    intake_bore_ratio_ig = 4
    exhaust_bore_ratio_ig = 1.7

    """------------------------------------"""

    """_____________________________________________________________________________________________________________________
    _____________________________________________SYSTEM MODELING____________________________________________________________
    _____________________________________________________________________________________________________________________"""

    q_out = q_out_ig
    toll_dq = 0.0001
    dq = toll_dq + 1

    """____To converge on the amount of heat needed to be discipated in gas cooler before getting to HX fluid_____"""

    m_dot_hx = .2  # Initial guess Just to do the first loop without error
    m_dot_r = m_dot_hx / 10
    m_dot_c = m_dot_hx / 20

    while dq >= toll_dq:

        """_______________________________________CRYOGENIC FLUID TDN STATES_________________________________________"""

        CTS = [TDN for i in range(num_data_C)]
        CTS_C = [TDN for i in range(num_data_C)]

        CTS[0] = TDN(p_st, 0, t_st, 0, 0, c)  # Starting from cryogenic tank
        CTS_C[0] = CTS[0].__copy__()
        tank_out = CTS[0]

        CTS[1] = TDN(p_exp_max, 0, CTS_C[0].dt(dt_chx), 0, 0, c)  # Heat exchange of cryogen with cooling compartment
        CTS_C[1] = CTS[1].__copy__()  # Flow expanded so to reach around 40 bar (Max piston can hold)
        exp_c_in = CTS[1]

        CTS[2] = TDN(pig[0], 0, t_exp_ig, 0, 0, c)
        CTS_C[2] = CTS[2].__copy__()

        CTS[3] = TDN(pig[1], 0, t_exp_ig, 0, 0, c)
        CTS_C[3] = CTS[3].__copy__()
        exp_c_out = CTS[3]

        CTS[4] = TDN(p_exp_min, 0, CTS_C[2].t, 0, 0,
                     c)  # Gas in near ambient pressure and temperature leaves piston exp
        CTS_C[4] = CTS[4].__copy__()
        exp_exh_c = CTS[4]

        CTS[5] = TDN(p_amb, 0, CTS_C[2].t, 0, 0, c)  # Gas in near ambient pressure and temperature leaves piston exp
        CTS_C[5] = CTS[5].__copy__()
        exp_amb_c = CTS[5]

        """----------------------------------------------"""

        C0 = TDN(p_amb, 0, t_amb, 0, 0, c)


        """-------------------------------------------------------------------------------------------------------"""

        """__________________________________________EXPANDER MODEL________________________________________________"""

        """ASSUMPTION: HX fluid is neglected in volume and pressure(Liquid vs N2 in gas form) and its pressence is only as 
        heat sourse to keep the process ISOTHERMAL"""

        """_____Ideal Isothermal Expansion_____"""
        """
        ASSUMPTINS : 
        1-Cryogen assumed to be ideal fluid 
        2-No effect of high RPM on open and closed valves
        3-Staedy states and no dynamics of two phase heat exchange involved
        4-Volume of HX neglected compared to Cryogen
        """

        ETSI = [TDN for i in
                range(num_data_exp_i)]  # Thermodynamic Instance List(Range can be varied on TDN states needed)
        ETS_C = [TDN for i in range(num_data_exp_i)]  # Thermodynamic Instance List copy

        ETSI[0] = TDN(p_exp_max, 0, CTS[1].t, 0, 0, c)  # 0   -Top Dead Centre-HX closed-Cryogen Open

        """ __Cryogen in__ """

        ETSI[1] = TDN(p_exp_max, 0, CTS[2].t, 0, 0, c)  # 60  -All closed

        """ __Power stroke_____Isothermal process___ """

        ETSI[2] = TDN(p_exp_exh, 0, CTS[2].t, 0, 0, c)  # 180 -Bottom Dead Centre Exhaust open

        """ __Return stroke__ """

        ETSI[3] = TDN(p_exp_min, 0, 0, 0, ETSI[2].d, c)  # exhaust still open and pressure drops to lowest point

        ETSI[4] = TDN(p_exp_min, 0, 0, 0, ETSI[0].d, c)  # 300 -HX open-Cryogen closed

        """ __HX fluid in__ """

        ETSI[5] = TDN(p_exp_max, 0, CTS[1].t, 0, 0, c)  # 360 -Top Dead Centre All closed

        w_exp = eta_exp_ig * (ETSI[1].p / ETSI[1].d * np.log((ETSI[1].d / ETSI[2].d)))  # ISOTHERMAL work extraction in expander

        """----------------------------------------------------------------------------------------------------------"""

        """Work from expander"""
        w_up = (pig[0] / TDN(pig[0], 0, exp_c_in.t, 0, 0, c).d) * np.log(
            (TDN(pig[0], 0, exp_c_in.t, 0, 0, c).d / TDN(pig[1], 0, t_exp_ig, 0, 0, c).d))
        w_down = (pig[1] / TDN(pig[1], 0, t_exp_ig, 0, 0, c).d) * np.log(
            (TDN(pig[1], 0, t_exp_ig, 0, 0, c).d / TDN(pig[2], 0, t_exp_ig, 0, 0, c).d))

        """______________________________VCC REFRIGERANT TDN STATES__________________________________________________"""

        """ VCC states for the R404A refrigerant """

        RTS = [TDN for i in
               range(num_data)]  # Thermodynamic Instance List (Range can be varied based on TDN states needed)
        RTS_C = [TDN for i in range(num_data)]  # Thermodynamic Instance List copy

        RTS[0] = TDN(p_low, PSI('H', 'Q', 1, 'P', p_low, r), 0, 0, 0, r)  # starting state from P low and quality 1
        RTS_C[0] = RTS[0].__copy__()

        RTS[1] = TDN(RTS_C[0].p, 0, RTS_C[0].dt(dt_superheat_ig), 0, 0,
                     r)  # Pressure drop in suction side entrance of compressor
        RTS_C[1] = RTS[1].__copy__()

        RTS[2] = TDN(RTS_C[1].dp(dp_suc), RTS_C[1].h, 0, 0, 0,
                     r)  # Pressure drop in suction side entrance of compressor
        RTS_C[2] = RTS[2].__copy__()
        comp_in = RTS[2]  # TDN properties at compressor entrance m3/kg

        s3_isen = TDN(p_high, 0, 0, RTS_C[2].s, 0, r)  # Isentropic compressor
        s3_isenc = s3_isen.__copy__()

        eta_isen_calc = 0.96-0.00046*rpm_ig+9.4e-8*rpm_ig**2+ 0.07*p_high/p_low-0.0018*(p_high/p_low)**2
        print(eta_isen_calc,eta_isen)

        RTS[3] = TDN(p_high, RTS_C[2].comp(eta_isen, s3_isenc.h), 0, 0, 0,
                     r)  # Isentropic efficiency of compressor correction
        RTS_C[3] = RTS[3].__copy__()
        comp_out = RTS[3]

        RTS[4] = TDN(RTS_C[3].dp(dp_dis), RTS_C[3].h, 0, 0, 0,
                     r)  # Pressure drop in discharge side entrance of compressor
        RTS_C[4] = RTS[4].__copy__()
        cond_in = RTS[4]

        RTS[5] = TDN(RTS[4].p, PSI('H', 'Q', 0, 'P', RTS_C[4].dp(dp_cond), r), 0, 0, 0,
                     r)  # H calculated by quality from PSI
        RTS_C[5] = RTS[5].__copy__()
        cond_out = RTS[5]

        RTS[6] = TDN(RTS[5].p, 0, RTS_C[5].dt(dt_subcool_ig), 0, 0, r)  # at constant pressure it is sub-cooled
        RTS_C[6] = RTS[6].__copy__()
        Molier_3 = RTS[6]

        RTS[7] = TDN(p_low + dp_evap, RTS[6].h, 0, 0, 0, r)  # Back calculated from p_low point with dp of evaporator
        RTS_C[7] = RTS[7].__copy__()
        evap_in = RTS[7]

        RTS[8] = TDN(RTS_C[7].dp(dp_evap), RTS_C[0].h, 0, 0, 0, r)
        evap_out = RTS[8]

        """----------------------------------------------"""

        R0 = TDN(p_amb, 0, t_amb, 0, 0, r)

        """----------------------------------------------------------------------------------------------------------"""

        COP_1 = (evap_out.h-evap_in.h)/(comp_out.h-comp_in.h)
        num_data_r_2r_2c_int = 17

        """__________________________________VCC       2 STAGE  with 2 Rec ____________________________"""

        """--------------RECIEVER ASSUMING STAEDY STATE---------- """
        """ VCC states for refrigerant """

        """ØØØØØØØØØØØØØØØØ     Equation to be solved for  finding mid pressure    ØØØØØØØØØØØØØØØØ"""

        p_mid = p_mt

        RTS2R2I = [TDN for i in range(num_data_r_2r_2c_int)]
        RTS2R2I_C = [TDN for i in range(num_data_r_2r_2c_int)]

        RTS2R2I[0] = TDN(p_high, PSI('H', 'Q', 0, 'P', p_high, r), 0, 0, 0, r)

        RTS2R2I[0] = TDN(p_high-dp_cond, 0, RTS2R2I[0].dt(dt_subcool_ig), 0, 0, r)

        RTS2R2I[1] = TDN(p_mid, RTS2R2I[0].h, 0, 0, 0, r)

        q_2r_2c_ht_i = RTS2R2I[1].q()

        RTS2R2I[2] = TDN(RTS2R2I[1].p, PSI('H', 'Q', 1, 'P', RTS2R2I[1].p, r), 0, 0, 0, r)  # gas from reciever outlet

        RTS2R2I[3] = TDN(RTS2R2I[1].p, PSI('H', 'Q', 0, 'P', RTS2R2I[1].p, r), 0, 0, 0,
                         r)  # Liquid from reciever outlet

        RTS2R2I[4] = TDN(p_low + dp_evap, RTS2R2I[3].h, 0, 0, 0, r)
        RTS2R2I_C[4] = RTS2R2I[4].__copy__()

        q_2r_2c_lt_i = RTS2R2I_C[4].q()

        RTS2R2I[5] = TDN(RTS2R2I[4].p, PSI('H', 'Q', 1, 'P', RTS2R2I[4].p , r), 0, 0, 0,
                         r)  # gas from reciever outlet

        RTS2R2I[6] = TDN(RTS2R2I[4].p, PSI('H', 'Q', 0, 'P', RTS2R2I[4].p, r), 0, 0, 0,
                         r)  # Liquid from reciever outlet
        RTS2R2I_C[6] = RTS2R2I[6].__copy__()

        RTS2R2I[7] = TDN(RTS2R2I_C[6].dp(dp_evap), PSI('H', 'Q', 1, 'P', RTS2R2I[6].p - dp_evap, r), 0, 0, 0, r)
        RTS2R2I_C[7] = RTS2R2I[7].__copy__()

        RTS2R2I[8] = TDN(RTS2R2I[7].p, 0, RTS2R2I_C[7].dt(dt_superheat_ig), 0, 0, r)
        RTS2R2I_C[8] = RTS2R2I[8].__copy__()

        """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
        RTS2R2I[9] = TDN((RTS2R2I[5].p * q_2r_2c_lt_i + RTS2R2I[8].p * (1 - q_2r_2c_lt_i)), (RTS2R2I[5].h * q_2r_2c_lt_i + RTS2R2I[8].h * (1 - q_2r_2c_lt_i)), 0, 0, 0, r)
        RTS2R2I_C[9] = RTS2R2I[9].__copy__()
        """-------------------------------"""

        RTS2R2I[10] = TDN(RTS2R2I_C[9].dp(dp_suc), RTS2R2I[9].h, 0, 0, 0, r)
        RTS2R2I_C[10] = RTS2R2I[10].__copy__()

        s3_isen_2r_2c_lp_i = TDN(p_mid, 0, 0, RTS2R2I[10].s, 0, r)  # Isentropic compressor
        s3_isenc_2r_2c_lp_i = s3_isen_2r_2c_lp_i.__copy__()

        RTS2R2I[11] = TDN(p_mid, RTS2R2I_C[10].comp(eta_isen, s3_isenc_2r_2c_lp_i.h), 0, 0, 0, r)
        RTS2R2I_C[11] = RTS2R2I[11].__copy__()

        RTS2R2I[12] = TDN(p_mid, 0, t_amb + dt_approach_air, 0, 0, r)
        RTS2R2I_C[12] = RTS2R2I[12].__copy__()

        """Energy balance on superheated gas from HT compressor and gas directly coming from reciever"""
        RTS2R2I[13] = TDN(RTS2R2I[2].p, (RTS2R2I[2].h * q_2r_2c_ht_i + RTS2R2I[12].h * (1 - q_2r_2c_ht_i)), 0, 0, 0, r)
        RTS2R2I_C[13] = RTS2R2I[13].__copy__()
        """-------------------------------"""

        s3_isen_2r_2c_hp_i = TDN(p_high, 0, 0, RTS2R2I[13].s, 0, r)  # Isentropic compressor
        s3_isenc_2r_2c_hp_i = s3_isen_2r_2c_hp_i.__copy__()

        RTS2R2I[14] = TDN(p_high, RTS2R2I_C[13].comp(eta_isen, s3_isenc_2r_2c_hp_i.h), 0, 0, 0, r)
        RTS2R2I_C[14] = RTS2R2I[14].__copy__()

        RTS2R2I[15] = TDN(RTS2R2I_C[14].dp(dp_dis), RTS2R2I_C[14].h, 0, 0, 0,
                          r)  # Pressure drop in discharge side entrance of compressor
        RTS2R2I_C[15] = RTS2R2I[15].__copy__()

        RTS2R2I[16] = TDN(RTS2R2I_C[15].dp(dp_cond), PSI('H', 'Q', 0, 'P', RTS2R2I[14].p - dp_cond, r), 0, 0, 0,
                          r)  # H calculated by quality
        RTS2R2I_C[16] = RTS2R2I[16].__copy__()

        cop_2r_2c_i = (RTS2R2I[7].h - RTS2R2I[6].h) * (1 - q_2r_2c_ht_i) * (1 - q_2r_2c_lt_i) / ((RTS2R2I[11].h - RTS2R2I[10].h) * (1 - q_2r_2c_lt_i) + (RTS2R2I[14].h - RTS2R2I[13].h))

        """---------------------------------------------------------------------------------------------------------"""

        """________________________________________HX FLUID STATES___________________________________________________"""

        HXTS = [TDN for i in range(num_data_hx)]
        HXTS_C = [TDN for i in range(num_data_hx)]

        HXTS[0] = TDN(p_amb, 0, CTS[3].t, 0, 0, hx)  # starting from piston expander output
        HXTS_C[0] = HXTS[0].__copy__()
        hx_in_r = HXTS[0]

        HXTS[1] = TDN(p_amb, HXTS[0].h + (m_dot_r / m_dot_hx) * (cond_in.h - cond_out.h), 0, 0, 0,
                      hx)  # HX fluid exchanges heat with in the condensor
        HXTS_C[1] = HXTS[1].__copy__()
        hx_out_r = HXTS[1]

        HXTS[2] = TDN(p_amb, 0, t_amb + dt_approach_air, 0, 0, hx)  # HX fluid exchanges heat with ambient(gas cooler)
        HXTS_C[2] = HXTS[2].__copy__()
        hx_out_amb = HXTS[2]

        HXTS[3] = TDN(p_exp_max, 0, HXTS_C[2].t, 0, 0, hx)  # HXF passes through pump to get high P
        HXTS_C[3] = HXTS[3].__copy__()  # Enough to enter piston
        hx_in_exp = HXTS[3]

        q_exp_c = exp_c_out.h - exp_c_in.h

        # print(HXTS[3].h - (w_up + w_down + q_exp_c)*m_dot_c/m_dot_hx)
        # HXTS[4] = TDN(exp_c_out.p, HXTS[3].h - (w_up + w_down + q_exp_c)*m_dot_c/m_dot_hx, 0, 0, 0, hx)  # HXF exchanges heat with piston expander

        HXTS[4] = TDN(exp_c_out.p, 0, t_amb, 0, 0,hx)  # HXF exchanges heat with piston expander
        HXTS_C[4] = HXTS[4].__copy__()
        hx_mix_exp = HXTS[4]

        HXTS[5] = TDN(CTS_C[3].p, 0, CTS[3].t, 0, 0, hx)  # HXF exchanges heat with piston expander
        HXTS_C[5] = HXTS[5].__copy__()
        hx_out_exp = HXTS[5]

        """----------------------------------------------"""

        HX0 = TDN(p_amb, 0, t_amb, 0, 0, hx)

        """----------------------------------------------------------------------------------------------------------"""

        """___________________________________________________________________________________________________________
        _____________________________________________M_Dot calculation__________________________________________________
        ___________________________________________________________________________________________________________"""

        """_______Per kg of Cryogen_________"""

        w_exp = eta_exp_ig * (ETSI[1].p / ETSI[1].d * np.log((ETSI[1].d / ETSI[2].d)))  # ISOTHERMAL work extraction in expander
        q_hx_c = exp_c_in.h - tank_out.h  # Heat exchanged in cryogenic HX
        q_exp_c = exp_c_out.h - exp_c_in.h  # Heat exchanged in isothermal expnsion for cryogen

        """_______Per kg of Refrigerant_________"""
        w_comp = ((RTS2R2I[11].h - RTS2R2I[10].h) * (1 - q_2r_2c_lt_i) + (RTS2R2I[14].h - RTS2R2I[13].h))  # Compressor work required (J/Kg refrierant)
        q_cond_r =(RTS2R2I[15].h - RTS2R2I[16].h)  # Heat discipated at the condensor in (J/Kg refrierant)
        q_evap_r = (RTS2R2I[7].h - RTS2R2I[6].h) * (1 - q_2r_2c_ht_i) * (1 - q_2r_2c_lt_i)  # Heat absorbed from Cooling effect  (J/Kg refrierant)

        """_______Per kg of HX fluid_________"""
        q_hx_r = hx_out_r.h - hx_in_r.h  # Heat exchanged for HX fluid at condensor
        q_hx_exp = hx_in_exp.h - hx_in_r.h  # Heat exchanged for HX fluid at expander
        q_hx_amb = hx_out_r.h - hx_out_amb.h

        """----------------------------------------------------------------------------------------------------------"""

        """_________________________________________Matrix solution_________________________________________________"""

        a = np.array([[w_up+w_down + q_exp_c, 0, -q_hx_exp],

                      [0, q_cond_r, -q_hx_exp],

                      [(w_up+w_down)*eta_exp_ig * eta_shaft * eta_aux, -w_comp, 0]])

        b = np.array([0, q_out, 0])

        try:
            x = np.linalg.solve(a, b)

        except np.linalg.LinAlgError:
            x = np.linalg.lstsq(a, b)[0]

        m_dot_c = x[0]  # Mass flow rate for cryogen in [kg/s]
        m_dot_r = x[1]  # Mass flow rate for Refrigrant in [kg/s]
        m_dot_hx = x[2]  # Mass flow rate for HX fluid in [kg/s]

        """----------------------------------------------------------------------------------------------------------"""

        """__________________________________________Cooling effect__________________________________________________"""

        cool_prod = (q_hx_c * m_dot_c + q_evap_r * m_dot_r)  # Total cooling effect
        q_out_cal = q_out * (cool_req / cool_prod)  # Heat discipated to ambient by gas cooler
        dq = abs(q_out_cal - q_out)  # Check for convergance
        q_out = q_out_cal  # Updated value of q_out

    """--------------------------------------------------------------------------------------------------------------"""

    """_____________________________________________________________________________________________________________________
    _____________________________________________Performance calculation____________________________________________________
    _____________________________________________________________________________________________________________________"""

    """Creating list of m_dots and boiling temp and Total entropy variation from start to end for various cryogens"""
    for xx in x:
        m_dot.append(xx)

    t_boil.append(t_st)

    s_tot.append((CTS[4].s-CTS[0].s) * m_dot_c)

    cop = (q_evap_r * m_dot_r) / (w_comp * m_dot_r)  # Coeficient of performance
    v_dot_suc = m_dot_r / comp_in.d  # Volume flow rate at the suction side of compressor  (M3/s)
    p_ratio = comp_out.p / comp_in.p  # Pressure ratio of outlet over inlet of compressor
    cool_kg_c = q_hx_c + (w_up+w_down)*eta_exp_ig * eta_aux * eta_shaft * cop_2r_2c_i  # Cooling per kg of cryogen

    w_vcc_sup = w_exp * eta_aux * eta_shaft * m_dot_c  # Work supplied to VCC
    w_vcc_req = w_comp * m_dot_r  # Work required by VCC

    """--------------------------------------------------------------------------------------------------------------"""

    """_____________________________________________________________________________________________________________________
    _____________________________________________Expander States____________________________________________________
    _____________________________________________________________________________________________________________________"""
    start = time.time()
    """___Expander engine characteristics___ : """
    """----------- Initial Guess ----------"""
    """DESIGN OF EXPANDER:
    1- The lower the valve opening degree the higher the expansion ratio and the higher the work extraction
    (There is mechanical limit for that and huge throttling loss(though can be countered by letting cryogen enter at higher
    pressure than initially planned based on structural limits of expander)
    2- The higher the ratio of connecting rod over shaft radius the higher the expansion ratio and the higher the extracted 
    work this trend remains until some point when it start to decrease .
    3- Bore diameter has no effect on the expansion ratio but controls the absolute value of p_max P_min and 
    V_max V_min as a result changing the value of extracted work  
    """
    f_ig = 20  # frequency

    l_con_ig = 0.1  # length of connecting rod to shaft
    r_shaft_ig = 0.01  # Shaft radius

    l_dead_ig = 0.001  # length of dead space top
    d_bore_ig = 0.06  # Bore diameter

    """------------------------------------"""

    """r_shaft_ig should be smaller than l_con_ig"""

    rpm = f_ig * 60

    l_str_ig = r_shaft_ig * 2  # length of expander chamber

    v_dead = (l_dead_ig * mt.pi * d_bore_ig ** 2) / 4  # Volume of dead space at the top
    v_tot = (l_str_ig * mt.pi * d_bore_ig ** 2) / 4  # Total volume

    theta_exp_l = []  # Angle of shaft
    z_exp_l = []  # Position of piston in expander
    v_exp_l = []  # Volume of expander
    p_exp_l = []  # Pressure withing piston expander

    """Degree at which the cryogen valve opens"""
    theta_0_ig = 0

    """Cryogen valve opening degree.It must be after zero not to exert negative power on piston,0 is the best choice"""
    theta_c_open_ig = 20

    """Degree at which the cryogen valve opens"""
    theta_exh_open_ig = 180

    """HX fluid valve opening degree"""
    theta_hx_open_ig = 330

    """Assumption : Total volume of HX at top dead state is taken as th volume of Top dead state """

    """ For the whole expantion process the volume of HX taken constant equal to that of top dead space """

    # v_dead_top = (m_dot_hx / hx_in_exp.d) / f_ig
    # v_dead_bot_ideal = (m_dot_c / exp_exh_c.d) / f_ig
    # cos_open = mt.cos((theta_c_open_ig * mt.pi) / 180)
    # v_ratio_open = mt.cos((theta_c_open_ig * mt.pi) / 180)
    # v_c_o_ideal = (m_dot_c / exp_c_out.d * v_ratio_open) / f_ig
    # p_power = exp_c_out.p / v_ratio_open
    # p_exh = p_power * (1 - cos_open) / 2
    #
    # v_pis_c_o = r_shaft_ig * (d_bore_ig / 2) ** 2 * mt.pi * (1 - cos_open)
    # v_pis_bot = r_shaft_ig * (d_bore_ig / 2) ** 2 * mt.pi * 2

    dead_top = 360  # Dead top angle
    dead_bot = 180  # Dead bottom angle

    for th in range(0, 360):

        z_c = r_shaft_ig * (1 - mt.cos(th / 180 * mt.pi)) + l_con_ig * (1 - mt.sqrt(1 - (r_shaft_ig / l_con_ig) ** 2 *
            mt.sin(th / 180 * mt.pi) ** 2))  # Piston distance from the dead top position
        v_c = (z_c + l_dead_ig) * mt.pi * d_bore_ig ** 2 / 4  # Total volume of piston expander based on rot degree of shaft

        """ _1_ Cryogen in__ """
        """ _2_ Power stroke_____Isothermal process___ """
        """ _3_ Return stroke__ """
        """ _4_ HX fluid in__ """
        """ _5_ Top dead center__"""

        if 0 <= th < theta_c_open_ig:
            p_exp_l.append(p_exp_max)

        elif theta_c_open_ig <= th < dead_bot:
            p_exp_l.append((ETSI[1].p / ETSI[1].d) * m_dot_c / (v_c * f_ig))

        elif dead_bot <= th < theta_hx_open_ig:
            p_exp_l.append(p_exp_exh)

        elif theta_hx_open_ig <= th < dead_top - 1:
            p_exp_l.append(p_exp_exh)

        else:
            p_exp_l.append(p_exp_max)

        z_exp_l.append(z_c)
        v_exp_l.append(v_c)
        theta_exp_l.append(th)

    p1 = p_exp_l[theta_c_open_ig]  # Highest pressure in the piston expander
    p2 = p_exp_l[dead_bot - 1]  # Lowest pressure in the piston expander
    d1 = TDN(p_exp_l[theta_c_open_ig], 0, t_exp_ig, 0, 0, c).d  # Highest density in the piston expander
    d2 = p2 / p1 * d1  # Lowest density in the piston expander

    w_exp_real = m_dot_c * (p1 / d1) * np.log(d1 / d2)  # Work produced based on pressures that achieved in expander
    dead_min = 1 / HXTS[3].d * f_ig * m_dot_hx  # Min dead volume required based on HX Fluid volume

    # ISOTHERMAL work extraction in expander
    RPM = (1 / d2) * m_dot_c / max(v_exp_l) * 60  # Revolution per second for piston expander


    """ Intake stroke   """
    d_bore_intake = d_bore_ig/intake_bore_ratio_ig
    gam_int = exp_c_in.cp() / exp_c_in.cv()
    p_star = exp_c_in.p * (2 / (gam_int + 1)) ** ((gam_int) / (gam_int - 1))
    m_dot_lim = c_disch * (d_bore_intake**2)/4*mt.pi * mt.sqrt(gam_int * exp_c_in.d * exp_c_in.p * ((2 / (gam_int + 1) ** ((gam_int + 1) / (gam_int - 1)))))
    w_exp_exact_ig = (TDN(p1, 0, t_amb, 0, 0, c).u() - TDN(p2, 0, t_amb, 0, 0, c).u() - TDN(p1, 0, t_amb, 0, 0, c).t * (TDN(p1, 0, t_amb, 0, 0, c).s - TDN(p2, 0, t_amb, 0, 0, c).s)) * m_dot_c

    """--------------------------------------------------------------------------------------------------------------"""

    """_____________________________________________________________________________________________________________________
      _____________________________________________EXPANDER OPTIMIZATION____________________________________________________
      _____________________________________________________________________________________________________________________"""

    """-------OBJECTIVE FUNCTION TO BE MINIMIZED------"""
    def objective(x):

        v_c_min = (x[5] * (1 - mt.cos(x[1] / 180 * mt.pi)) + x[4] * (1 - mt.sqrt(1 - (x[5] / x[4]) ** 2 *mt.sin(x[1] / 180 * mt.pi) ** 2))+ x[3]) * mt.pi * x[2] ** 2 / 4

        v_c_max = (x[5] * (1 - mt.cos(theta_exh_open_ig / 180 * mt.pi)) + x[4] * (1 - mt.sqrt(1 - (x[5] / x[4]) ** 2 *mt.sin(theta_exh_open_ig / 180 * mt.pi) ** 2))+ x[3]) * mt.pi * x[2] ** 2 / 4
        # print(- v_c_min * mt.log(v_c_max / v_c_min)*TDN(0, 0, t_amb, 0, 1/v_c_min/x[0] * m_dot_c, c).p)
        return - v_c_min * mt.log(v_c_max / v_c_min)*TDN(0, 0, t_amb, 0, 1/v_c_min/x[0] * m_dot_c, c).p

    """--------CONSTRAINTS THAT SHOULD BE HELD TRUE -----------"""
    def cont0(x):  # INEQ  # f(x)>= 0 inequality always means >=

        # print(x[4] - x[5] * 2)
        return x[4] - x[5] * 2

    def cont1(x):  # EQ

        v_c_min = (x[5] * (1 - mt.cos(x[1] / 180 * mt.pi)) + x[4] * (1 - mt.sqrt(1 - (x[5] / x[4]) ** 2 * mt.sin(x[1] / 180 * mt.pi) ** 2)) + x[3]) * mt.pi * x[2] ** 2 / 4

        # print(v_c_min - m_dot_c / exp_c_out.d / x[0])
        return v_c_min - m_dot_c / exp_c_out.d / x[0]

    def cont2(x):  # INEQ

        """ Intake stroke   """
        d_bore_intake = x[2] / x[6]
        gam = exp_c_in.cp() / exp_c_in.cv()
        m_dot_lim = c_disch * (d_bore_intake ** 2) / 4 * mt.pi * mt.sqrt(gam * exp_c_in.d * exp_c_in.p * ((2 / (gam + 1) ** ((gam + 1) / (gam - 1)))))

        # print(m_dot_lim * x[1] / 360 - m_dot_c)
        return m_dot_lim * x[1] / 360 - m_dot_c

    def cont3(x):  # INEQ

        v_c_min = (x[5] * (1 - mt.cos(x[1] / 180 * mt.pi)) + x[4] * (1 - mt.sqrt(1 - (x[5] / x[4]) ** 2 * mt.sin(x[1] / 180 * mt.pi) ** 2)) + x[3]) * mt.pi * x[2] ** 2 / 4

        # print(p_lim_exp - TDN(0, 0, t_amb, 0, 1/v_c_min/x[0] * m_dot_c, c).p)
        return p_lim_exp - TDN(0, 0, t_amb, 0, 1/v_c_min/x[0] * m_dot_c, c).p

    """___Expander engine characteristics___ : """
    """----------- Initial Guess ----------"""
    """ Zero = [0 f    1 theta       2 d_bore   3 l_dv    4 l_cr    5 r_sh    6  intake_bore_ratio ] """
    zero = [f_ig, theta_c_open_ig, d_bore_ig, l_dead_ig, l_con_ig, r_shaft_ig, intake_bore_ratio_ig]

    """------------------------------------"""
    print(zero)
    print('Initial design objective function value   :', objective(zero))

    """------RANGE OF VARIABLE BONDS --------"""
    bon0 = (f_ig, 100)
    bon1 = (theta_c_open_ig, 35)
    bon2 = (1e-2, 1e-1)
    bon3 = (1e-3, 1e-2)
    bon4 = (1e-2, l_con_ig)
    bon5 = (1e-2, 1e-1)
    bon6 = (3, 10)

    bon = [bon0, bon1, bon2, bon3, bon4, bon5, bon6]

    """ CONSTRANTS EQ OR INEQ SPESIFIED """

    con0 = {'type': 'ineq', 'fun': cont0}  # INEQ  always means >= 0
    con1 = {'type': 'eq', 'fun': cont1}
    con2 = {'type': 'ineq', 'fun': cont2}
    con3 = {'type': 'ineq', 'fun': cont3}

    con = [con0, con1, con2, con3]

    """------SOLUTION FOR OPTIMIZATION PROBLEM-----"""
    sol = minimize(objective, zero, method='SLSQP', bounds=bon, constraints=con)
    print('Optimized design objective function value :', objective([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4], sol.x[5], sol.x[6]]))
    print(sol)
    end = time.time()
    print("CODE RUNTIME for expander optimization is  : " + str(end - start))
    print('--- Results of the Optimization Problem ---')
    print('x[0]: ', 'Frequency of piston expander :      ', sol.x[0])
    print('x[1]: ', 'Opening valve degree for cryogen :  ', sol.x[1])
    print('x[2]: ', 'Bore diameter :                     ', sol.x[2])
    print('x[3]: ', 'Dead volume on top :                ', sol.x[3], '+', m_dot_hx/HXTS[2].d/mt.pi / sol.x[0] /  (sol.x[2]**2)*4, ' For HXF')
    print('x[4]: ', 'Length of connecting rod :          ', sol.x[4])
    print('x[5]: ', 'Shaft diameter :                    ', sol.x[5])
    print('x[6]: ', 'Intake valve to bore radius ratio : ', sol.x[6])
    print('-------------------------------------------')

    """------OPTIMIZED VALUES REPLACED WITH INITIAL GUESS------"""
    f = sol.x[0]
    theta_c_open = sol.x[1]
    theta_c_open = int(round(sol.x[1]))
    d_bore = sol.x[2]
    l_dead = sol.x[3]
    l_con = sol.x[4]
    r_shaft = sol.x[5]
    intake_bore_ratio = sol.x[6]


    #
    # """ TO REMOVE THE  EFFECT OF PISTON EXPANDER"""
    #
    # f = f_ig
    # theta_c_open = theta_c_open_ig
    # d_bore = d_bore_ig
    # l_dead = l_dead_ig
    # l_con = l_con_ig
    # r_shaft = r_shaft_ig
    # intake_bore_ratio = intake_bore_ratio_ig
    # theta_c_open = int(round(theta_c_open))
    #
    # """REMOVE LATER"""
    #




    p_star = exp_c_in.p * (2 / (gam_int + 1)) ** ((gam_int) / (gam_int - 1))
    l_str = r_shaft * 2  # length of expander chamber

    v_dead = (l_dead * mt.pi * d_bore ** 2) / 4  # Volume of dead space at the top
    v_tot = (l_str_ig * mt.pi * d_bore ** 2) / 4  # Total volume

    theta_exp_l = []  # Angle of shaft
    z_exp_l = []  # Position of piston in expander
    v_exp_l = []  # Volume of expander
    p_exp_l = []  # Pressure withing piston expander

    """ For the whole expantion process the volume of HX taken constant equal to that of top dead space """

    v_dead_top = (m_dot_hx / hx_in_exp.d) / f_ig
    v_dead_bot_ideal = (m_dot_c / exp_exh_c.d) / f_ig
    cos_open = mt.cos((theta_c_open * mt.pi) / 180)
    v_ratio_open = mt.cos((theta_c_open * mt.pi) / 180)
    v_c_o_ideal = (m_dot_c / exp_c_out.d * v_ratio_open) / f_ig
    p_power = exp_c_out.p / v_ratio_open
    p_exh = p_power * (1 - cos_open) / (2)

    v_pis_c_o = r_shaft * (d_bore / 2) ** 2 * mt.pi * (1 - cos_open)
    v_pis_bot = (r_shaft*2 + m_dot_hx / HXTS[2].d / mt.pi / f / (d_bore ** 2) * 4 + l_dead) * (d_bore / 2) ** 2 * mt.pi

    dead_top = 360  # Dead top angle
    dead_bot = 180  # Dead bottom angle

    for th in range(0, 360):

        z_c = r_shaft * (1 - mt.cos(th / 180 * mt.pi)) + l_con * (1 - mt.sqrt(1 - (r_shaft / l_con) ** 2 *
                mt.sin(th / 180 * mt.pi) ** 2))  # Piston distance from the dead top position
        v_c = (z_c + l_dead) * mt.pi * d_bore ** 2 / 4  # Total volume of piston expander based on rot degree of shaft

        """ _1_ Cryogen in__ """
        """ _2_ Power stroke_____Isothermal process___ """
        """ _3_ Return stroke__ """
        """ _4_ HX fluid in__ """
        """ _5_ Top dead center__"""

        if 0 <= th < theta_c_open:
            p_exp_l.append(p_exp_max)

        elif theta_c_open <= th < dead_bot:
            p_exp_l.append((ETSI[1].p / ETSI[1].d) * m_dot_c / (v_c * f_ig))

        elif dead_bot <= th < theta_hx_open_ig:
            p_exp_l.append(p_exp_exh)

        elif theta_hx_open_ig <= th < dead_top - 1:
            p_exp_l.append(p_exp_exh)

        else:
            p_exp_l.append(p_exp_max)

        z_exp_l.append(z_c)
        v_exp_l.append(v_c)
        theta_exp_l.append(th)

    p1 = p_exp_l[theta_c_open]  # Highest pressure in the piston expander
    p2 = p_exp_l[dead_bot - 1]  # Lowest pressure in the piston expander
    CTS[1] = TDN(pig[0], 0, exp_c_in.t, 0, 0,c)
    CTS[3] = TDN(p1, 0, CTS_C[2].t, 0, 0,c)
    CTS[4] = TDN(p2, 0, CTS_C[2].t, 0, 0,c)

    d1 = TDN(p_exp_l[theta_c_open], 0, t_exp_ig, 0, 0, c).d  # Highest density in the piston expander
    d2 = p2 / p1 * d1  # Lowest density in the piston expander

    intake = TDN(p1, 0, t_amb, 0, 0, c)
    exhaust = TDN(p2, 0, t_amb, 0, 0, c)

    w_exp_real = m_dot_c * (p1 / d1) * np.log(d1 / d2)  # Work produced based on pressures that achieved in expander
    w_exp_real_TDN = (intake.u() - exhaust.u() - intake.t * (intake.s - exhaust.s)) * m_dot_c

    dead_min = 1 / HXTS[3].d * f_ig * m_dot_hx  # Min dead volume required based on HX Fluid volume

    # ISOTHERMAL work extraction in expander
    RPM = (1 / d2) * m_dot_c / max(v_exp_l) * 60  # Revolution per second for piston expander
    RPM_f = rpm

    """--------------------------------------------------------------------------------------------------------------"""
    """--------------------------------------------------------------------------------------------------------------"""
    """--------------------------------------------------------------------------------------------------------------"""

    """ Exhaust stroke   """
    d_bore_exh = d_bore_ig / exhaust_bore_ratio_ig
    gam_exh = exhaust.cp()/exhaust.cv()
    p_star_exh = exhaust.p * (2/(gam_exh+1))**((gam_exh)/(gam_exh-1))

    Y_in = 1-(0.351 + 0.256 *(1/intake_bore_ratio)**4 + 0.93 *(1/intake_bore_ratio)**8)*(1-(p_amb/p1)**(1/gam_int))
    Y_out = 1-(0.351 + 0.256 *(1/exhaust_bore_ratio_ig)**4 + 0.93 *(1/exhaust_bore_ratio_ig)**8)*(1-(p_amb/p2)**(1/gam_int))

    m_dot_lim_exh = c_disch * (d_bore_exh**2)/4*mt.pi * mt.sqrt(gam_exh*exhaust.d*exhaust.p*((2/(gam_exh+1)**((gam_exh+1)/(gam_exh-1)))))
    p_drop_exh = 0.5 * (HXTS[4].d+exhaust.d)*(1-(1/exhaust_bore_ratio_ig)**4)*((m_dot_hx/HXTS[4].d)/(((180 - theta_c_open) / 360)*0.75* d_bore_exh**2 *mt.pi/4))**2 +0.5 * (exhaust.d)*(1-(1/exhaust_bore_ratio_ig)**4)*((m_dot_c/exhaust.d)/(((180 - theta_c_open) / 360)*0.75*Y_out* d_bore_exh**2 *mt.pi/4))**2

    """ Intake stroke   """
    d_bore_intake = sol.x[2] / sol.x[6]
    gam_int = exp_c_in.cp() / exp_c_in.cv()
    m_dot_lim = c_disch * (d_bore_intake ** 2) / 4 * mt.pi * mt.sqrt(gam_int * exp_c_in.d * exp_c_in.p * ((2 / (gam_int + 1) ** ((gam_int + 1) / (gam_int - 1)))))
    p_drop_int = 0.5 * (intake.d) * (1 - (1 / intake_bore_ratio) ** 4) * ((m_dot_c / intake.d) / (((theta_c_open) / 360) * 0.75 *Y_in * d_bore_intake ** 2 * mt.pi / 4)) ** 2

    """_______________________________________________________________________________________________________________
    _____________________________________________2nd Law Analysis____________________________________________________
    _______________________________________________________________________________________________________________"""

    """_____________________________________________CRYOGEN____________________________________________________"""
    """Separation of Cryogen from HX fluid is not considered here, since two fluids are inimicible and at different phases """

    """Cryogenic heat exchanger-----Assuming constant CP"""
    s_q_c_hx = q_hx_c * (1 / mlog(tank_out.t, exp_c_in.t) - 1 / t_evap) * m_dot_c
    s_tot_c_hx = - (tank_out.s - exp_c_in.s) * m_dot_c

    """Ass : for simpilicity the expansion process is devided into two different process:
    1- Heat exchange of Cryogen with HX fluid up until reaching the equilibrium TEMP
    2- Isothermal expansion of cryogen in the expander(HX fluid keeps the temp constant)
    """

    """Cryogenic HX fluid mixture in expander-----Assuming cons,tant CP"""
    s_q_exp_c = q_exp_c * (1 / mlog(exp_c_in.t, exp_c_out.t) - 1 / mlog(hx_in_exp.t, hx_mix_exp.t)) * m_dot_c
    s_tot_exp_c = (exp_c_out.s - exp_c_in.s) * m_dot_c

    """Cryogenic HX fluid mixture expansion in expander-----Assuming constant temperature"""
    """Ass :
    Cryogen behaves idealy in the range p_max to p_exh
    No heat loss occurs in the expander
    """

    s_exp_c = (exp_c_out.z() + exp_exh_c.z()) / 2 * (1 - eta_exp_ig) * rgas / exp_c_in.m() * mt.log(
        p_ratio) * m_dot_c  # Isothermal expander entropy
    s_tot_exh_c = (exp_exh_c.s - exp_c_out.s) * m_dot_c

    """Expansion from lowest p reached within expander to ambient"""
    s_exp_exh = (exp_amb_c.z() + exp_exh_c.z()) / 2 * rgas / exp_c_in.m() * mt.log(
        p_exp_l[dead_bot - 1] / p_amb) * m_dot_c  # Isothermal expander entropy
    s_tot_exp_exh = (exp_amb_c.s - exp_exh_c.s) * m_dot_c

    """_____________________________________________HX FLUID_______________________________________________________"""
    """HX Fluid exchanging heat in condensor"""
    """ The heat exchange entropy is already calculated in refrigeration condensor"""
    s_tot_hx_cond = (hx_out_r.s - hx_in_r.s) * m_dot_hx

    """HX Fluid exchanging heat within the air-cooled heat exchanger"""
    s_hx_cond = q_out * m_dot_hx * (1 / t_amb - 1 / mlog(hx_out_r.t, hx_out_amb.t))
    s_tot_hx_cond = (hx_out_r.s - hx_out_amb.s) * m_dot_hx

    """HX Fluid pumped to reach the pressure high enough to enter expader"""
    s_w_pump_hx = m_dot_hx / hx_in_exp.d * (hx_in_exp.p - hx_out_amb.p) / hx_in_exp.t * ((1 - eta_pump) / (eta_pump))
    s_tot_pump_hx = -m_dot_hx * (hx_in_exp.s - hx_out_amb.s)

    """HX Fluid exchanging heat cryogen in expander"""
    """This consist of 2 phases of heat exchange :
    1-Heat transfer for getting cryogen from t_evap to expander temperature
    2-Heat transfer to keep expander at isothermal condition during expansion

    ASS: All heat transfer occurs with entropy gen as if all happened within first process
    """
    s_hx_c_exp = m_dot_hx * (hx_in_exp.h - hx_out_exp.h) * (1 / mlog(exp_exh_c.t, exp_c_in.t) - 1 / mlog(hx_in_exp.t, hx_out_exp.t))
    s_tot_hx_c_exp = m_dot_hx * (hx_in_exp.s - hx_out_exp.s)
    """------------------------------------------------------------------------------------------------------------"""

    """______________________________________ Refrigeration cycle _____________________________________________"""

    """Pressure drops"""
    s_dp_suc = m_dot_r * rgas / R0.m() * comp_in.z() * mt.log(1 + dp_suc / comp_in.p)
    s_dp_dis = m_dot_r * rgas / R0.m() * comp_out.z() * mt.log(1 + dp_dis / cond_in.p)
    s_dp_evap = m_dot_r * rgas / R0.m() * comp_in.z() * mt.log(1 + dp_evap / evap_out.p)
    s_dp_cond = m_dot_r * rgas / R0.m() * (cond_in.z() + cond_out.z()) / 2 * mt.log(1 + dp_cond / cond_out.p)

    """Compressor """
    s_comp_r = m_dot_r * rgas / comp_in.m() * (comp_in.z() + comp_out.z()) / 2 * ((1 - eta_isen) / eta_isen) * mt.log(
        p_ratio)

    """Condensor heat exchange entropy gen"""
    s_cond_r = q_cond_r * m_dot_r * (1 / mlog(hx_out_r.t, hx_in_r.t) - 1 / mlog(cond_in.t, cond_out.t))

    """TXV isentalpic throttling"""
    s_trot_r = m_dot_r * rgas / cond_out.m() * mt.log(evap_in.p / evap_out.p)

    """Evaporator heat exchange entropy"""
    s_evap_r = q_evap_r * m_dot_r * (1 / mlog(evap_in.t, evap_out.t) - 1 / t_reefer)

    """--------------------------------------------------------------------------------------------------------------"""
    """__Sorted List of entropy generations : __"""
    s = [s_q_c_hx, s_q_exp_c, s_hx_cond, s_exp_exh, s_w_pump_hx, s_dp_suc, s_dp_dis, s_comp_r, s_dp_evap,
         s_dp_cond, s_cond_r, s_trot_r, s_evap_r]

    for value in s:
        ss.append(value)
    s_list = [['s_q_c_hx', s_q_c_hx], ['s_q_exp_c', s_q_exp_c], ['s_exp_exh', s_exp_exh], ['s_hx_cond', s_hx_cond],
              ['s_w_pump_hx', s_w_pump_hx], ['s_dp_suc', s_dp_suc], ['s_dp_dis', s_dp_dis], ['s_dp_evap', s_dp_evap],
              ['s_dp_cond', s_dp_cond], ['s_comp_r', s_comp_r], ['s_cond_r', s_cond_r], ['s_trot_r', s_trot_r],
              ['s_evap_r', s_evap_r]]
    sort_s = sort(s_list)
    six = ['s_q_c_hx', 's_exp_exh', 's_comp_r', 's_cond_r', 's_w_pump_hx', 's_q_exp_c']

    for i in range(len(six)):
        for j in range(len(s_list)):
            a = s_list[j]
            if a[0] == six[i]:
                list.append(a[1])

    s_sum.append(sum([s_q_c_hx, s_q_exp_c, s_exp_exh, s_hx_cond, s_w_pump_hx, s_dp_suc, s_dp_dis, s_dp_evap, s_dp_cond, s_comp_r,s_cond_r, s_trot_r, s_evap_r]))

    """Work based on initial and final condition internal energy and entropy::"""
    exp_in = TDN(p1, 0, t_exp_ig, 0, 0, c)
    exp_out = TDN(p2, 0, t_exp_ig, 0, 0, c)

    w_exp_max = (intake.u()-exp_amb_c.u()-intake.t*(intake.s-exp_amb_c.s))*m_dot_c

    w_exp_max_ex = (intake.ex()-exp_amb_c.ex())*m_dot_c

    w_exp_exact =(exp_in.u()-exp_out.u()-exp_in.t*(exp_in.s-exp_out.s))*m_dot_c
    w_up = (pig[0] / TDN(pig[0], 0,t_exp_ig, 0, 0, c).d) * np.log(
        (TDN(pig[0], 0, t_exp_ig, 0, 0, c).d / TDN(p1, 0, t_exp_ig, 0, 0, c).d))
    w_down = w_exp_exact/m_dot_c

    #
    #
    #
    #
    # """ TO MAKE THE 2ND EXPANDER EFFECT TO ZERO"""
    # w_up = 0
    # """REMOVE LATER"""
    #
    #
    #

    cool_kg_c = q_hx_c + (w_up+w_down)*eta_exp_ig * eta_aux * eta_shaft * cop_2r_2c_i  # Cooling per kg of cryogen

    #
    #
    #
    # """ TO MAKE THE VCC single stage"""
    # cool_kg_c = q_hx_c + (w_up + w_down) * eta_exp_ig * eta_aux * eta_shaft * COP_1  # Cooling per kg of cryogen
    # """REMOVE LATER"""
    #
    #
    #

    """--------------------------------------------------------------------------------------------------------------"""
    """EXERGY input to system"""
    ex_hx = (HXTS[3].ex() - HXTS[5].ex())*m_dot_hx
    ex_c = (CTS[0].ex()- CTS[4].ex())*m_dot_c
    """_________________________________________________________________________________________________________________
    ________________________________________________PRINT_______________________________________________________________
    _________________________________________________________________________________________________________________"""

    # print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Refrigerant'],
    #       "\n", TDN.prop(test), "\n", TDN.prop(test))

    if p1 > p_lim_exp+5e5:
        print("!!!! WARNING!!!! :Wrong selection of open cryogen volume(not enough space for Cryogen)")

    if m_dot_c > m_dot_lim/360*theta_c_open:
        print("!!!! WARNING!!!! :Wrong selection of intake diameter or valve opening degree(choked flow limits m_dot")

    if m_dot_lim_exh * ((theta_hx_open_ig-theta_exh_open_ig) / 360) < m_dot_c:
        print("!!!! WARNING!!!! :Wrong selection of exhaust diameter or valve opening degree(choked flow limits m_dot")
        print(m_dot_lim_exh * ((180 - theta_c_open) / 360), ' < ', m_dot_c)

    print("Lowest pressure in the expander       : " + str(p2/p_amb)+ "[atm]")
    print("Highest pressure in the expander      : " + str(p1/p_amb)+ "[atm]")
    print("Expansion pressure ratio of expander  : " + str(p1 / p2))
    print("Pressure drop during intake stroke    : " + str(p_drop_int / p_amb) + "[atm]")
    print("Pressure drop during exhaust stroke   : " + str(p_drop_exh/p_amb) + "[atm]")
    print("Total volume of piston expander       : " + str(v_pis_bot) + "[m^3]")
    print("Mass flow rate for Cryogen in [kg/s]      : " + str(m_dot_c) + ' of '+ c)
    print("Mass flow rate for Refrigerant in [kg/s]  : " + str(m_dot_r) + ' of '+ r)
    print("Mass flow rate for HX fluid in [kg/s]     : " + str(m_dot_hx) + ' of '+ hx)

    print("The limiting cryogen at choking condition  " + str(m_dot_lim / 360 * 30))

    print("Cooling effect (of VCC) per Kg of refrigerant  : " + str(q_evap_r))
    print("TOTAL COOLING effect per Kg of Cryogen         : " + str(cool_kg_c))
    print("TOTAL COOLING output                           : " + str(cool_prod))
    print("COP 2comp 2rec int is :  ____ " + str(cop_2r_2c_i)+ "____")

    print("Volume flow rate at the suction of comp    :  " + str(v_dot_suc))
    print("Pressure ratio in compressor               :  " + str(p_ratio))

    print("Heat rejected to ambient at air-cooled HX  : " + str(q_hx_amb * m_dot_hx))
    print("Heat transfered to piston expander by HXF  : " + str((HXTS[3].h - HXTS[4].h) * m_dot_hx))
    print("Heat exchanged at condenser     :            " + str(q_cond_r*m_dot_r))
    print("Heat expander :                              " + str((exp_c_out.h-exp_c_in.h)*m_dot_c))

    print("RPM   : " + str(RPM) + "      &   Frequency is : " + str(f))

    print("Work required for VCC  :   " + str(w_vcc_req))
    print("W_exp_real based on U-S, Optimal design : " + str((w_up+w_down)*m_dot_c*eta_exp_ig * eta_shaft * eta_aux))
    print()
    print("Mass flow rate per hour of cryogen : " + str(m_dot_c*3600))
    print()
    print("Max work that can be extracted from the flow : " + str(w_exp_max), w_exp_max_ex)

    print("Z_in exp  : " + str(CTS[2].z())+ "      Z_out exp : " + str(CTS[3].z()))

    print("Sum of entropies : " + str(sum(
        [s_q_c_hx, s_q_exp_c, s_exp_exh, s_hx_cond, s_w_pump_hx, s_hx_c_exp, s_dp_suc, s_dp_dis, s_dp_evap, s_dp_cond,
         s_comp_r, s_cond_r, s_trot_r, s_evap_r])))

    print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Refrigerant'])
    print('-------------------2 stage + 2 receiver + 1 intercooler VCC------------')

    i = 0
    for til2r2i in RTS2R2I:
        if i < 10:
            print(i, ' , ', TDN.prop(til2r2i))
        else:
            print(i, ', ', TDN.prop(til2r2i))
        i += 1
    i = 0
    # print('-------------------Cryogen TDN states------------')
    # print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Cryogenic Fluid'])
    for cts in CTS:
        print(i, ' , ', TDN.prop(cts))
        i += 1
    i = 0

    # print('-------------------HXF TDN states------------')
    # print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'HX Fluid'])
    for hxt in HXTS:
        print(i, ' , ', TDN.propl(hxt))
        i += 1
    i = 0

    """--------------------------------------------------------------------------------------------------------------"""
    print()
    print('"""-----------------------------------------------------------------------------------------------"""')
    print()

    """_________________________________________________________________________________________________________________
    ____________________________________________________PLOTS___________________________________________________________
    _________________________________________________________________________________________________________________"""

    """___________________________Thermodynamic states plot____________________________________________________ : """

    h_plot = [] # First plots the lines and then 'ro' spesifies red dots as the TDN states
    p_plot = []

    h_plot_2r_2c_i = []
    p_plot_2r_2c_i = []

    for i in range(len(RTS)):
        h_plot.append(RTS[i].h)
        p_plot.append(RTS[i].p)

    for i in range(len(RTS2R2I)):
        h_plot_2r_2c_i.append(RTS2R2I[i].h)
        p_plot_2r_2c_i.append(RTS2R2I[i].p)

    h_plot_2r_2c_i.append(RTS2R2I[0].h)
    p_plot_2r_2c_i.append(RTS2R2I[0].p)

    plt.figure('Vapor Compressionn cycle TDN states of refrigerant' )
    h_p = plt.semilogy(h_plot, p_plot)
    h_p_2r_2c_i = plt.semilogy(h_plot_2r_2c_i, p_plot_2r_2c_i,'gold')
    h_p_ro = plt.semilogy(h_plot, p_plot, 'ro')
    h_p_ro_2r_2c_i = plt.semilogy(h_plot_2r_2c_i, p_plot_2r_2c_i, 'cd')
    plt.title('Vapor Compressionn cycle TDN states of refrigerant R404A')
    plt.xlabel('Entalpy [j/kg]')
    plt.ylabel('Log Pressure, Log(P) [Pa]')

    """--------------------------------------------------------------------------------------------------------------"""

    """    Liquid & Vapor saturation plot"""

    sat_liq_h = []
    sat_vap_h = []
    sat__p = []

    """______________Two-Phase plot from T_min to P_Critical for Refrigerant__________________"""
    for p in range(int(round(PSI('P', 'Q', 0, 'T', PSI(r, 'Tmin'), r), -4)+1e4), int(round(PSI(r, 'pcrit'), -4)-1e4), 1000):

        sat_liq_h.append(PSI('H', 'Q', 0, 'P', p, r))  # At vapor quality of 0 the sat_liquid line is obtained
        sat_vap_h.append(PSI('H', 'Q', 1, 'P', p, r))  # At vapor quality of 1 the sat_vapor line is obtained
        sat__p.append(p)

    plot_liq_sat = plt.semilogy(sat_liq_h, sat__p)
    plot_vap_sat = plt.semilogy(sat_vap_h, sat__p)
    plt.grid()

    """CRYOGENIC FLUID"""
    """   Thermodynamic states plot"""
    h_plot_c = []
    h_plot_ro_c = []  # First plots the lines and then 'ro' specifies red dots as the TDN states
    p_plot_c = []
    p_plot_ro_c = []

    for i in range(len(CTS)):
        h_plot_c.append(CTS[i].h)
        h_plot_ro_c.append(CTS[i].h)
        p_plot_c.append(CTS[i].p)
        p_plot_ro_c.append(CTS[i].p)

    plt.figure('Thermodynamic states plot of cryogenic flow of  '+ str(C[cryo]))

    h_p_c = plt.semilogy(h_plot_c, p_plot_c)
    h_p_ro_c = plt.semilogy(h_plot_ro_c, p_plot_ro_c, 'ro')
    plt.title('TDN states of cryogen ' + c)
    plt.xlabel('Entalpy [j/kg]')
    plt.ylabel('Log Pressure, Log(P) [Pa]')

    """--------------------------------------------------------------------------------------------------------------"""
    leg = C[0], C[1], C[2], C[3]
    """   Liquid & Vapor saturation plot"""

    sat_liq_h_c = []
    sat_vap_h_c = []
    sat__p_c = []

    """______________Two-Phase plot from T_min to P_Critical for Cryogen__________________"""
    for p in range(int(round(PSI('P', 'Q', 0,'T', PSI(c,'Tmin'), c), -4)+1e4), int(round(PSI(c, 'pcrit'), -4)-1e4), 1000):
        sat_liq_h_c.append(PSI('H', 'Q', 0, 'P', p, c))  # At vapor quality of 0 the sat_liquid line is obtained
        sat_vap_h_c.append(PSI('H', 'Q', 1, 'P', p, c))  # At vapor quality of 1 the sat_vapor line is obtained
        sat__p_c.append(p)

    plot_liq_sat_c = plt.semilogy(sat_liq_h_c, sat__p_c)
    plot_vap_sat_c = plt.semilogy(sat_vap_h_c, sat__p_c)
    plt.grid()

    """ HX FLUID """
    """  Thermodynamic states plot"""
    h_plot_hx = []
    h_plot_ro_hx = []  # First plots the lines and then 'ro' spesifies red dots as the TDN states

    t_plot_hx = []
    t_plot_ro_hx = []

    for i in range(len(HXTS)):
        h_plot_hx.append(HXTS[i].h)
        t_plot_hx.append(HXTS[i].t+ cryo * 0.2)
        h_plot_ro_hx.append(HXTS[i].h)
        t_plot_ro_hx.append(HXTS[i].t+ cryo * 0.2)

    plt.figure(' Thermodynamic states plot of HX fluid for variouss cyogens')
    h_p_hx = plt.plot(h_plot_hx, t_plot_hx)
    h_p_ro_hx = plt.plot(h_plot_ro_hx, t_plot_ro_hx, 'ro')
    plt.title('TDN states of HX fluid  Ethylene Glycol 60% ')
    plt.xlabel('Entalpy j/kg]')
    plt.ylabel('Temperature,  [K]')

    plt.grid()  # Grid on


    """___________________________....  TS PLOT ....________________________________________"""

    """CRYOGENIC FLUID"""
    """   Thermodynamic states plot"""
    s_plot_c = []
    s_plot_ro_c = []  # First plots the lines and then 'ro' specifies red dots as the TDN states
    t_plot_c = []
    t_plot_ro_c = []

    for i in range(len(CTS)):
        s_plot_c.append(CTS[i].s)
        s_plot_ro_c.append(CTS[i].s)
        t_plot_c.append(CTS[i].t)
        t_plot_ro_c.append(CTS[i].t)

    plt.figure('Thermodynamic states plot of cryogenic flow  ' + str(C[cryo]))

    s_p_c = plt.plot(s_plot_c, t_plot_c)
    s_p_ro_c = plt.plot(s_plot_ro_c, t_plot_ro_c, 'ro')
    plt.title('TDN states of cryogen ' + c)
    plt.xlabel('Entropy [j/kg/k]')
    plt.ylabel('Temperature, T [K]')


    leg = C[0], C[1], C[2], C[3]
    """   Liquid & Vapor saturation plot"""

    sat_liq_s_c = []
    sat_vap_s_c = []
    sat__t_c = []

    """______________Two-Phase plot from T_min to P_Critical for Cryogen__________________"""
    for t in range(int(PSI(c, 'Tmin')+2), int(round(PSI(c, 'Tcrit')))):
        sat_liq_s_c.append(PSI('S', 'Q', 0, 'T', t, c))  # At vapor quality of 0 the sat_liquid line is obtained
        sat_vap_s_c.append(PSI('S', 'Q', 1, 'T', t, c))  # At vapor quality of 1 the sat_vapor line is obtained
        sat__t_c.append(t)

    sat_liq_s_c.append(PSI('S', 'Q', 0, 'T', PSI(c, 'Tcrit'), c))  # At vapor quality of 0 the sat_liquid line is obtained
    sat_vap_s_c.append(PSI('S', 'Q', 1, 'T', PSI(c, 'Tcrit'), c))  # At vapor quality of 1 the sat_vapor line is obtained
    sat__t_c.append(PSI(c, 'Tcrit'))


    plot_liq_sat_c = plt.plot(sat_liq_s_c, sat__t_c)
    plot_vap_sat_c = plt.plot(sat_vap_s_c, sat__t_c)
    plt.grid()


    """--------------------------------------------------------------------------------------------------------------"""

    """ Cryogen PV plot """

    v_plot_cc = []
    v_plot_ro_cc = []  # First plots the lines and then 'ro' spesifies red dots as the TDN states
    p_plot_cc = []
    p_plot_ro_cc = []

    """    Expander PV diagram"""
    leg = C[0], C[1], C[2], C[3]
    plt.figure('Expander PV diagram')
    plt.semilogy(v_exp_l, p_exp_l)
    plt.title('Expander  P-V plot for various shaft angle ')
    plt.xlabel('Volume [m3/kg]')
    plt.ylabel('Log Pressure, Log(P)  [pa]')
    plt.legend(leg)

    """--------------------------------------------------------------------------------------------------------------"""

"""Entropy generations for Nitrogen"""

plt.figure('Entropy generations for Nitrogen')

labels = 'Cryogenic HX', 'Expander heat', 'HXF heat', 'Expander exhaust', 'HXF pump', 'Dp suction side', 'Dp discharge side', 'Compressor', 'Dp evaporator', 'Dp condenser', 'Condenser', 'Trottling ', 'Evaporator'
sizes = ss[0:int(len(ss) / 4)]
colors = ['orange', 'brown', 'lightcoral', 'coral', 'indianred', 'silver', 'rosybrown', 'lightcoral', 'darksalmon', 'salmon', 'grey', 'rosybrown','orangered', 'gold']
explode = (0.15, 0.3, 0, 0, 0, 0.3, 0.15, 0, 0, 0.14, 0, 0, 0)  # explode 1st slice
plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)

plt.axis('equal')
plt.legend(bbox_to_anchor=(0.1, 0.3))

"""--------------------------------------------------------------------------------------------------------------"""

"""2ND LAW EFFICIENCY"""
for a in range(len(s_sum)):
    eta_2_law[a] = 1 - s_sum[a] / s_tot[a]

print(eta_2_law)

"""--------------------------------------------------------------------------------------------------------------"""

"""   Top 5 Entropy generations for various cryogens  """
Cryogen0 = (s_sum[0], list[0], list[1], list[2], list[3], list[4], list[5])
Cryogen1 = (s_sum[1], list[6], list[7], list[8], list[9], list[10], list[11])
Cryogen2 = (s_sum[2], list[12], list[13], list[14], list[15], list[16], list[17])
Cryogen3 = (s_sum[3], list[18], list[19], list[20], list[21], list[22], list[23])

ind = np.arange(len(Cryogen0))  # the x locations for the groups
width = 0.4  # the width of the bars

fig, ax = plt.subplots()
rects0 = ax.bar(ind - width / 2, Cryogen0, width / 2, label=C[0])
rects1 = ax.bar(ind - width / 4, Cryogen1, width / 2, label=C[1])
rects2 = ax.bar(ind + width / 4, Cryogen2, width / 2, label=C[2])
rects3 = ax.bar(ind + width / 2, Cryogen3, width / 2, label=C[3])

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Entropy generations [J/kg/K/s]')
ax.set_xlabel('Components causing highest losses')
ax.set_title('Entropy generation rates for various cryogen used(Losses) in [J/kg/K/s]')
ax.set_xticks(ind)
ax.set_xticklabels(('Total losses', 'Cryogenic HX', 'Expander exhaust', 'Compressor', 'Condensor', 'Pumps', 'Expander heat'))
ax.legend()

"""--------------------------------------------------------------------------------------------------------------"""

"""   Mass flow rate in [kg / s] for various cryogens  """
Cryogen0 = (m_dot[0], m_dot[1], m_dot[2])
Cryogen1 = (m_dot[3], m_dot[4], m_dot[5])
Cryogen2 = (m_dot[6], m_dot[7], m_dot[8])
Cryogen3 = (m_dot[9], m_dot[10], m_dot[11])

ind = np.arange(len(Cryogen0))  # the x locations for the groups
width = 0.5  # the width of the bars

fig, ax = plt.subplots()
re0 = ax.bar(ind - width / 2, Cryogen0, width / 2, label=C[0])
re1 = ax.bar(ind - width / 4, Cryogen1, width / 2, label=C[1])
re2 = ax.bar(ind + width / 4, Cryogen2, width / 2, label=C[2])
re3 = ax.bar(ind + width / 2, Cryogen3, width / 2, label=C[3])

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Mass flow rate in [kg / s]')
ax.set_xlabel('Mass flow rate for various flows')
ax.set_title('Mass flow rate in [kg / s] for ')
ax.set_xticks(ind)
ax.set_xticklabels(('Cryogenic', 'Refrigerant', 'HX fluid'))
ax.legend()

"""--------------------------------------------------------------------------------------------------------------"""

"""   Boiling temperature [k] for various cryogens  """

width = 0.4  # the width of the bars

fig, ax = plt.subplots()
rec0 = ax.bar(1, t_boil[0], width * 2, label=C[0])
rec1 = ax.bar(2, t_boil[1], width * 2, label=C[1])
rec2 = ax.bar(3, t_boil[2], width * 2, label=C[2])
rec3 = ax.bar(4, t_boil[3], width * 2, label=C[3])

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Boiling temperature [k] for various cryogens')
ax.set_xlabel('Various cryogens')
ax.set_title('Boiling temperature [k] for various cryogens ')
ax.set_xticks((1, 2, 3, 4))
ax.set_xticklabels((C[0], C[1], C[2], C[3]))
ax.legend()

"""--------------------------------------------------------------------------------------------------------------"""

"""   2nd law efficiency for various cryogens  """

width = 0.4  # the width of the bars

fig, ax = plt.subplots()
recer0 = ax.bar(1, eta_2_law[0], width * 2, label=C[0])
recer1 = ax.bar(2, eta_2_law[1], width * 2, label=C[1])
recer2 = ax.bar(3, eta_2_law[2], width * 2, label=C[2])
recer3 = ax.bar(4, eta_2_law[3], width * 2, label=C[3])

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('2nd law efficiency [-] for various cryogens')
ax.set_xlabel('Various cryogens')
ax.set_title('2nd law efficiency for various cryogens ')
ax.set_xticks((1, 2, 3, 4))
ax.set_xticklabels((C[0], C[1], C[2], C[3]))
ax.legend()

"""--------------------------------------------------------------------------------------------------------------"""

"""         Entropy comparison for all cryogens         """
fig, axs = plt.subplots(2, 2)

labels = 'Cryogenic HX', 'Expander heat', 'HXF heat', 'Expander exhaust', 'HXF pump', 'Dp suction side', 'Dp discharge side', 'Compressor', 'Dp evaporator', 'Dp condenser', 'Condenser', 'Trottling ', 'Evaporator'
sizes0 = ss[0:int(len(ss) / 4)]
sizes1 = ss[int(len(ss) / 4):int(len(ss) / 4) * 2]
sizes2 = ss[int(len(ss) / 4) * 2:int(len(ss) / 4) * 3]
sizes3 = ss[int(len(ss) / 4) * 3:int(len(ss) / 4) * 4]

explode = (0, 0.15, 0, 0, 0, 0.25, 0.15, 0, 0, 0.35, 0, 0.15, 0.2)  # explode 1st slice
axs[0, 0].pie(sizes0, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)
axs[0, 0].set_title(C[0])
axs[0, 1].pie(sizes1, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)
axs[0, 1].set_title(C[1])
axs[1, 0].pie(sizes2, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)
axs[1, 0].set_title(C[2])
axs[1, 1].pie(sizes3, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)
axs[1, 1].set_title(C[3])
plt.legend(bbox_to_anchor=(0.5, 0.7), bbox_transform=ax.transAxes)

"""--------------------------------------------------------------------------------------------------------------"""

"""         Entropy comparison for all cryogens         """

fig, axs = plt.subplots(2, 2)

labels = 'Cryogenic HX', 'Expander heat', 'HXF heat', 'Expander exhaust', 'HXF pump', 'Dp suction side', 'Dp discharge side', 'Compressor', 'Dp evaporator', 'Dp condenser', 'Condenser', 'Trottling ', 'Evaporator', '2nd Law Efficiency'
sizes0 = ss[0:int(len(ss) / 4)]
sizes1 = ss[int(len(ss) / 4):int(len(ss) / 4) * 2]
sizes2 = ss[int(len(ss) / 4) * 2:int(len(ss) / 4) * 3]
sizes3 = ss[int(len(ss) / 4) * 3:int(len(ss) / 4) * 4]

sizes0 += [s_tot[0] - s_sum[0]]
sizes1 += [s_tot[1] - s_sum[1]]
sizes2 += [s_tot[2] - s_sum[2]]
sizes3 += [s_tot[3] - s_sum[3]]

explode = (0, 0.15, 0, 0, 0, .25, 0.5, 0, 0, 0.75, 0, 0, 0.15, 0.2)   # explode 1st slice
colors = ['orange', 'lightsalmon', 'lightcoral', 'coral', 'indianred', 'silver', 'rosybrown', 'lightcoral', 'darksalmon', 'salmon', 'grey', 'rosybrown','orangered', 'gold']

axs[0, 0].pie(sizes0, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=0)
axs[0, 0].set_title(C[0])
axs[0, 1].pie(sizes1, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=0)
axs[0, 1].set_title(C[1])
axs[1, 0].pie(sizes2, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=0)
axs[1, 0].set_title(C[2])
axs[1, 1].pie(sizes3, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=0)
axs[1, 1].set_title(C[3])

plt.legend(bbox_to_anchor=(0.5, 0.7), bbox_transform=ax.transAxes)

"""--------------------------------------------------------------------------------------------------------------"""

"""NITROGEN """
explode = (0, 0.15, 0, 0, 0, 0.35, 0.25, 0, 0, 0.35, 0, 0, 0.15, 0.2)  # explode 1st slice

plt.figure(C[0])
plt.pie(sizes0, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=0)
plt.legend(bbox_to_anchor=(0.5, 0.7), bbox_transform=ax.transAxes)

"""--------------------------------------------------------------------------------------------------------------"""

"""_________________________________________________PLOT SHOW________________________________________________________"""

""" PLT.SHOW() :: shows multiple figures at the same time, if the figures are intended to be shown as one closes there
should be plt.show after each individual plot, Any comment after plt.show will be ran as soon as the figure is closed"""

plt.show()

"""--------------------------------------------------------------------------------------------------------------"""

