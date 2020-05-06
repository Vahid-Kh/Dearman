import numpy as np
import matplotlib.pyplot as plt
from Dearman.TDN import PSI, TDN, mlog, mt, sort

"""________________________________________________NOTICE:___________________________________________________________"""
"""________SOME PARAMETERS HAVE JUST BEEN ASSUMED AT THE MOMENT WHICH LATER SHOULD BE COMPUTED ITERATIVELY_______ """
"""____________________________________________All dimensions in SI unit_____________________________________________"""

"""_____________________________________________________________________________________________________________________
_____________________________________________INITIAL CONDITION__________________________________________________________
_____________________________________________________________________________________________________________________"""


# R = ['R404A', 'R22', 'n-Propane']

# R = ['R404A']
# R = ['R22']
R = ['n-Propane']
C = ['Nitrogen', 'Methane', 'Ethane', 'Hydrogen']
HX = ['INCOMP::MEG[0.6]', 'INCOMP::AEG[0.6]']

l ,  ss,  qq,  ww = [], [], [], []
s_tot, s_sum, m_dot = [], [], []
cop_l, cop_r_l, cop_r_2_l, cop_2r_2c_l, cop_2r_2c_i_l, cop_2r_2c_i_m_l = [], [], [], [], [], []

t_amb_in = 273.15
t_l = []
m_l = []
f_l = []
cop_o_l = []
num_temp = 38
cool_kg_c_l = []
eta_2_law = []
m_l_r = []
eta_l_r = []
cop_o_l_r = []
w_out_exp = []
for refrig in range(len(R)):
    print('"""____________________________' + str(R[refrig]) + '_____________________________________"""')
    print()
    for temp in range(num_temp):
        t_amb = t_amb_in + temp
        t_l.append(temp)
        """____FLUIDS____ : """
        r = R[refrig]  # Refrigerant in the refrigeration cycle
        c = C[0]  # Cryogenic fluid used for cooling (Nitrogen in this case)
        """   __ Model 1 :  # CoolProp mixture Water/ASHRAE, Ethylene Glycol 60% (volume based)
              __ Model 2 :  # Aake Melinder Properties of Secondary Working fluids published 2010 IIR"""

        hx = HX[0]

        """_________ Reference Condition __________ : """
        t0 = 293.15
        p0 = 1.01325e5
        rgas = 8.31446261815324             #Gas constant J⋅K−1⋅mol−1

        """___Temperature [K]___ : """

        dt_approach = 4  # Min Temp difference at watercooled heat exchanger
        dt_approach_air = 8  # Min Temp difference at aircooled heat exchanger
        t_evap = 243.15  # Evaporator temperature (T_Min) - 30 C
        t_reefer = t_evap + dt_approach

        t_st = PSI('T', 'P', p0, 'Q', 0, c)  # Cryogen storage tank temperature(Ambient pressure and quality of 0)
        t_chx = t_evap  # Cryogen after heat exchange in CHX
        t_exp_ig = t_amb  # Piston expander average temperature

        dt_subcool_ig = - 5  # DT of sub-cooling in (k) (with minus sign)
        dt_superheat_ig = 8  # DT of superheat in (k)

        dt_chx = t_chx - t_st
        """------------------------------------"""
        """----------- Initial Guess ----------"""
        t_cond = t_amb+dt_approach  # Condenser temperature (T_Max)  + 25 C

        """------------------------------------"""
        """___Pressure [Pa]___  : """

        p_evap = PSI('P', 'Q', 1, 'T', t_evap, r)  # Evaporator pressure (P_Min)

        p_cond = PSI('P', 'Q', 1, 'T', t_cond, r)  # Condenser pressure (P_Max)


        p_low = p_evap  # p_low of VCC (pa)
        p_high = p_cond  # p_high of VCC (pa)

        p_amb = p0  # Ambient pressure
        p_st = p_amb * 2  # Cryogen pressure at storage tank
        p_exp_max = 4.5e6  # Max pressure at the piston expander
        p_lim_exp = 4.8e6  # Mechanical limiting pressure for expader
        p_exp_exh = p_amb * 1.5  # Exhaust pressure at the piston expander
        p_exp_min = p_exp_exh  # Min pressure at the piston expander

        dp_suc = 2e4  # Pressure drop at the suction side of compressor (pa)
        dp_dis = 2e4  # Pressure drop at the discharge side of compressor (pa)
        dp_evap = 5e4  # Pressure drop at the Evaporator in (pa)
        dp_cond = 2e4  # Pressure drop at the condenser in (pa)
        dp_liq = 1e4  # Pressure drop at the liquid line (pa)
        """------------------------------------"""

        """___Number of Data Points___ : """

        num_data_r = 8  # Number of data points for Refrigerant, 4 Component, 2 dp, 2 dt, 1

        num_data_C = 5  # Number of data points for Cryogen open cycle

        num_data_hx = 6  # Number of data points for HX closed cycle

        num_data_exp_i = 6  # Number of data points Expander engine per seconds

        num_data_exp_r = 100  # Number of data points Expander engine in each RPM
        """------------------------------------"""

        """___Efficiencies___ : """
        eta_hx = 0.95      # Efficiency of heat exchanger
        eta_shaft = 0.97   # Expander to Compressor
        eta_exp_ig = 0.7   # Efficiency of expander
        eta_isen = 0.65    # Isentropic efficiency of the compressor
        eta_aux = 0.95     # Auxillary work required ratio to total work
        eta_pump = 0.9     # Efficiency of heat exchanger
        c_disch = 0.75  # Discharge coeficient at cryogen inlet of expander
        """------------------------------------"""

        """----------- Cooling required ----------"""
        cool_req = 15e3  # Cooling capacity required from VCC amd cryogenic HX  (W)
        """------------------------------------"""

        """----------- Initial Guess ----------"""
        ig_mdot_c = (500 * 10 ** (-3) / 0.00123) / 6 / 3600  # Initial guess of N2 M_dot entering piston expander in each revolution

        q_out_ig = cool_req / 3  # Heat exchanged in condensor gas cooler befrore exchanging heat with HX fluid

        intake_bore_ratio_ig = 4
        exhaust_bore_ratio_ig = 2.5
        """------------------------------------"""

        """_____________________________________________________________________________________________________________
        _____________________________________________SYSTEM MODELING____________________________________________________
        _____________________________________________________________________________________________________________"""

        q_out = q_out_ig
        toll_dq = 0.001
        dq = toll_dq + 1

        """____To converge on the amount of heat needed to be discipated in gas cooler before getting to HX fluid_____"""

        m_dot_hx = 1                    # Initial guess Just to do the first loop without error
        m_dot_r = m_dot_hx/2

        while dq >= toll_dq:

            """___________________________________________CRYOGENIC FLUID TDN STATES____________________________________"""

            CTS = [TDN for i in range(num_data_C)]
            CTS_C = [TDN for i in range(num_data_C)]

            CTS[0] = TDN(p_st, 0, t_st, 0, 0, c)  # Starting from cryogenic tank
            CTS_C[0] = CTS[0].__copy__()
            tank_out = CTS[0]

            CTS[1] = TDN(p_exp_max, 0, CTS_C[0].dt(dt_chx), 0, 0, c)  # Heat exchange of cryogen with cooling compartment
            CTS_C[1] = CTS[1].__copy__()  # Flow expanded so to reach around 40 bar (Max piston can hold)
            exp_c_in = CTS[1]

            if TDN(0, 0, t_exp_ig, 0, CTS[1].d, c).p > p_lim_exp:

                CTS[2] = TDN(p_lim_exp, 0, t_exp_ig, 0, 0, c)
            else:
                CTS[2] = TDN(0, 0, t_exp_ig, 0, CTS[1].d, c)  # Tempereture of cryogen rises to near ambient at near isochoric process

            CTS_C[2] = CTS[2].__copy__()
            exp_c_out = CTS[2]

            CTS[3] = TDN(p_exp_min, 0, CTS_C[2].t, 0, 0, c)  # Gas in near ambient pressure and temperature leaves piston exp
            CTS_C[3] = CTS[3].__copy__()
            exp_exh_c = CTS[3]

            CTS[4] = TDN(p_amb, 0, CTS_C[3].t, 0, 0, c)  # Gas in near ambient pressure and temperature leaves piston exp
            CTS_C[4] = CTS[4].__copy__()
            exp_amb_c = CTS[4]

            """----------------------------------------------"""

            C0 = TDN(p_amb, 0, t_amb , 0, 0, c)

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

            ETSI = [TDN for i in range(num_data_exp_i)]  # Thermodynamic Instance List(Range can be varied on TDN states needed)
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

            w_exp = eta_exp_ig * (
                        ETSI[1].p / ETSI[1].d * np.log((ETSI[1].d / ETSI[2].d)))  # ISOTHERMAL work extraction in expander

            """---------------------------------------------------------------------------------------------------------"""

            """__________________________________VCC REFRIGERANT TDN STATES_____________________________________________"""

            """
            ASSUMPTINS :
             
            1- Staedy states and no dynamics  involved
            2- The refrigeration cycle always run in sub-critical region
            """

            """ VCC states for the R404A refrigerant """

            RTS = [TDN for i in range(num_data_r)]  # Thermodynamic Instance List (Range can be varied based on TDN states needed)
            RTS_C = [TDN for i in range(num_data_r)]  # Thermodynamic Instance List copy

            RTS[0] = TDN(p_low, PSI('H', 'Q', 1, 'P', p_low, r), 0, 0, 0,r)  # starting state from P low and quality 1
            RTS_C[0] = RTS[0].__copy__()

            RTS[1] = TDN(RTS_C[0].p, 0, RTS_C[0].dt(dt_superheat_ig), 0, 0,
                         r)  # Pressure drop in suction side entrance of compressor
            RTS_C[1] = RTS[1].__copy__()

            RTS[2] = TDN(RTS_C[1].dp(dp_suc), RTS_C[1].h, 0, 0, 0, r)  # Pressure drop in suction side entrance of compressor
            RTS_C[2] = RTS[2].__copy__()
            comp_in = RTS[2]  # TDN properties at compressor entrance m3/kg

            s3_isen = TDN(p_high, 0, 0, RTS_C[2].s, 0, r)  # Isentropic compressor
            s3_isenc = s3_isen.__copy__()

            RTS[3] = TDN(p_high, RTS_C[2].comp(eta_isen, s3_isenc.h), 0, 0, 0,
                         r)  # Isentropic efficiency of compressor correction
            RTS_C[3] = RTS[3].__copy__()
            comp_out = RTS[3]

            RTS[4] = TDN(RTS_C[3].dp(dp_dis), RTS_C[3].h, 0, 0, 0, r)  # Pressure drop in discharge side entrance of compressor
            RTS_C[4] = RTS[4].__copy__()
            cond_in = RTS[4]

            RTS[5] = TDN(RTS[4].p, PSI('H', 'Q', 0, 'P', RTS_C[4].dp(dp_cond), r), 0, 0, 0,r)  # H calculated by quality from PSI
            RTS_C[5] = RTS[5].__copy__()
            cond_out = RTS[5]

            RTS[6] = TDN(RTS[5].p, 0, RTS_C[5].dt(dt_subcool_ig), 0, 0, r)  # at constant pressure it is sub-cooled
            RTS_C[6] = RTS[6].__copy__()
            Molier_3 = RTS[6]

            RTS[7] = TDN(p_low + dp_evap, RTS[6].h, 0, 0, 0, r)  # Back calculated from p_low point with dp of evaporator
            RTS_C[7] = RTS[7].__copy__()
            evap_in = RTS[7]

            """-----------------------------------------------------------------------------------------------------"""

            num_data_r_r = 10  # Number of data points for Refrigerant, 4 Component, 2 dp, 2 dt, 1 ; 1stage with reciever

            """__________________________________VCC with RECIEVER  REFRIGERANT TDN STATES____________________________"""
            """--------------RECIEVER ASSUMING STAEDY STATE---------- """
            """ VCC states for refrigerant """
            RTSR = [TDN for i in range(num_data_r_r)]
            RTSR_C = [TDN for i in range(num_data_r_r)]

            rec_liq = TDN(RTS[7].p, PSI('H', 'Q', 0, 'P', RTS[7].p, r), 0, 0, 0, r)

            RTSR[0] = rec_liq
            evap_in_rec = RTSR[0]
            RTSR_C[0] = RTSR[0].__copy__()
            evap_in_rec = RTSR[0]
            q_req = RTS[7].q()

            RTSR[1] = TDN(RTSR_C[0].dp(dp_evap), PSI('H', 'Q', 1, 'P', RTSR[0].p - dp_evap, r), 0, 0, 0, r)
            RTSR_C[1] = RTSR[1].__copy__()
            evap_out_rec = RTSR[1]
            evap_out = evap_out_rec

            RTSR[2] = TDN(RTSR_C[1].p, 0, RTSR_C[1].dt(dt_superheat_ig), 0, 0,
                          r)  # Pressure drop in suction side entrance of compressor
            RTSR_C[2] = RTSR[2].__copy__()
            rec_gas = TDN(RTSR[2].p, PSI('H', 'Q', 1, 'P', RTSR[2].p, r), 0, 0, 0, r)

            """Energy balance on superheated gas from evaporator and gas directly coming from reciever"""
            RTSR[3] = TDN(RTSR[2].p, (rec_gas.h * q_req + RTSR_C[2].h * (1 - q_req)), 0, 0, 0, r)
            RTSR_C[3] = RTSR[3].__copy__()

            RTSR[4] = TDN(RTSR_C[3].dp(dp_suc), RTSR_C[3].h, 0, 0, 0,
                          r)  # Pressure drop in suction side entrance of compressor
            RTSR_C[4] = RTSR[4].__copy__()
            comp_in_rec = RTSR[4]  # TDN properties at compressor entrance m3/kg

            s3_isen_rec = TDN(p_high, 0, 0, RTSR_C[4].s, 0, r)  # Isentropic compressor
            s3_isenc_rec = s3_isen_rec.__copy__()

            RTSR[5] = TDN(p_high, RTSR_C[4].comp(eta_isen, s3_isenc_rec.h), 0, 0, 0,
                          r)  # Isentropic efficiency of compressor correction
            RTSR_C[5] = RTSR[5].__copy__()
            comp_out_rec = RTSR[5]

            RTSR[6] = TDN(RTSR_C[5].dp(dp_dis), RTSR_C[5].h, 0, 0, 0,
                          r)  # Pressure drop in discharge side entrance of compressor
            RTSR_C[6] = RTSR[6].__copy__()
            cond_in_rec = RTSR[6]

            RTSR[7] = TDN(RTSR[6].p, PSI('H', 'Q', 0, 'P', RTSR_C[6].dp(dp_cond), r), 0, 0, 0,
                          r)  # H calculated by quality from PSI
            RTSR_C[7] = RTSR[7].__copy__()
            cond_out_rec = RTSR[7]

            RTSR[8] = TDN(RTSR[7].p, 0, RTSR_C[7].dt(dt_subcool_ig), 0, 0, r)  # at constant pressure it is sub-cooled
            RTSR_C[8] = RTSR[8].__copy__()
            Molier_3 = RTSR[8]

            RTSR[9] = TDN(p_low + dp_evap, RTSR[8].h, 0, 0, 0, r)  # Back calculated from p_low point with dp of evaporator
            RTSR_C[9] = RTSR[9].__copy__()


            R0 = TDN(p_amb, 0, t_amb , 0, 0, r)

            """---------------------------------------------------------------------------------------------------------"""

            num_data_r_r_2 = 13

            """__________________________________VCC       2 STAGE  with RECIEVER ____________________________"""

            """--------------RECIEVER ASSUMING STAEDY STATE---------- """
            """ VCC states for refrigerant """

            p_mid = (p_low + p_high) / 2

            RTSR2 = [TDN for i in range(num_data_r_r_2)]
            RTSR2_C = [TDN for i in range(num_data_r_r_2)]

            RTSR2[0] = TDN(p_high, PSI('H', 'Q', 0, 'P', p_high, r), 0, 0, 0, r)

            RTSR2[0] = TDN(p_high, 0, RTSR2[0].dt(dt_subcool_ig), 0, 0, r)

            RTSR2[1] = TDN(p_mid, RTSR2[0].h, 0, 0, 0, r)

            q_req_2 = RTSR2[1].q()

            RTSR2[2] = TDN(RTSR2[1].p, PSI('H', 'Q', 1, 'P', RTSR2[1].p, r), 0, 0, 0, r)  # gas from reciever outlet

            RTSR2[3] = TDN(RTSR2[1].p, PSI('H', 'Q', 0, 'P', RTSR2[1].p, r), 0, 0, 0, r)  # Liquid from reciever outlet

            RTSR2[4] = TDN(p_low+dp_evap, RTSR2[3].h, 0, 0, 0, r)
            RTSR2_C[4] = RTSR2[4].__copy__()

            RTSR2[5] = TDN(RTSR2_C[4].dp(dp_evap), PSI('H', 'Q', 1, 'P', RTSR2[4].p - dp_evap, r), 0, 0, 0, r)
            RTSR2_C[5] = RTSR2[5].__copy__()

            RTSR2[6] = TDN(RTSR2[5].p, 0, RTSR2_C[5].dt(dt_superheat_ig), 0, 0, r)
            RTSR2_C[6] = RTSR2[6].__copy__()

            RTSR2[7] = TDN(RTSR2_C[6].dp(dp_suc), RTSR2[6].h, 0, 0, 0, r)
            RTSR2_C[7] = RTSR2[7].__copy__()

            s3_isen_rec_2_lp = TDN(p_mid, 0, 0, RTSR2[7].s, 0, r)  # Isentropic compressor
            s3_isenc_rec_2_lp = s3_isen_rec_2_lp.__copy__()

            RTSR2[8] = TDN(p_mid, RTSR2_C[7].comp(eta_isen, s3_isenc_rec_2_lp.h), 0, 0, 0, r)
            RTSR2_C[8] = RTSR2[8].__copy__()

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
            RTSR2[9] = TDN(RTSR2[2].p, (RTSR2[2].h * q_req_2 + RTSR2[8].h * (1 - q_req_2)), 0, 0, 0, r)
            RTSR2_C[9] = RTSR2[9].__copy__()
            """-------------------------------"""

            s3_isen_rec_2_hp = TDN(p_high, 0, 0, RTSR2[9].s, 0, r)  # Isentropic compressor
            s3_isenc_rec_2_hp = s3_isen_rec_2_hp.__copy__()

            RTSR2[10] = TDN(p_high, RTSR2_C[9].comp(eta_isen, s3_isenc_rec_2_hp.h), 0, 0, 0, r)
            RTSR2_C[10] = RTSR2[10].__copy__()

            RTSR2[11] = TDN(RTSR2_C[10].dp(dp_dis), RTSR2_C[10].h, 0, 0, 0, r)  # Pressure drop in discharge side entrance of compressor
            RTSR2_C[11] = RTSR2[11].__copy__()

            RTSR2[12] = TDN(RTSR2_C[11].dp(dp_cond), PSI('H', 'Q', 0, 'P', RTSR2[11].p - dp_cond, r), 0, 0, 0,r)  # H calculated by quality
            RTSR2_C[12] = RTSR2[12].__copy__()

            cop_r_2 = (RTSR2[6].h-RTSR2[4].h)*(1-q_req_2)/((RTSR2[8].h-RTSR2[7].h)*(1-q_req_2)+(RTSR2[10].h-RTSR2[9].h))

            """---------------------------------------------------------------------------------------------------------"""

            num_data_r_2r_2c = 16

            """__________________________________VCC 2 STAGE  2 Rec 1 Intercooler ____________________________"""

            """--------------RECIEVER ASSUMING STAEDY STATE---------- """
            """ VCC states for refrigerant """

            p_mid = (p_low + p_high) / 2

            RTS2R2 = [TDN for i in range(num_data_r_2r_2c)]
            RTS2R2_C = [TDN for i in range(num_data_r_2r_2c)]

            RTS2R2[0] = TDN(p_high, PSI('H', 'Q', 0, 'P', p_high, r), 0, 0, 0, r)

            RTS2R2[0] = TDN(p_high, 0, RTS2R2[0].dt(dt_subcool_ig), 0, 0, r)

            RTS2R2[1] = TDN(p_mid, RTS2R2[0].h, 0, 0, 0, r)

            q_2r_2c_ht = RTS2R2[1].q()

            RTS2R2[2] = TDN(RTS2R2[1].p, PSI('H', 'Q', 1, 'P', RTS2R2[1].p, r), 0, 0, 0, r)  # gas from reciever outlet

            RTS2R2[3] = TDN(RTS2R2[1].p, PSI('H', 'Q', 0, 'P', RTS2R2[1].p, r), 0, 0, 0, r)  # Liquid from reciever outlet

            RTS2R2[4] = TDN(p_low + dp_evap, RTS2R2[3].h, 0, 0, 0, r)
            RTS2R2_C[4] = RTS2R2[4].__copy__()

            q_2r_2c_lt = RTS2R2_C[4].q()

            RTS2R2[5] = TDN(RTS2R2[4].p, PSI('H', 'Q', 1, 'P', RTS2R2[4].p-dp_evap, r), 0, 0, 0, r)  # gas from reciever outlet

            RTS2R2[6] = TDN(RTS2R2[4].p, PSI('H', 'Q', 0, 'P', RTS2R2[4].p, r), 0, 0, 0, r)  # Liquid from reciever outlet
            RTS2R2_C[6] = RTS2R2[6].__copy__()

            RTS2R2[7] = TDN(RTS2R2_C[6].dp(dp_evap), PSI('H', 'Q', 1, 'P', RTS2R2[6].p - dp_evap, r), 0, 0, 0, r)
            RTS2R2_C[7] = RTS2R2[7].__copy__()

            RTS2R2[8] = TDN(RTS2R2[7].p, 0, RTS2R2_C[7].dt(dt_superheat_ig), 0, 0, r)
            RTS2R2_C[8] = RTS2R2[8].__copy__()

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
            RTS2R2[9] = TDN(RTS2R2[8].p, (RTS2R2[5].h * q_2r_2c_lt + RTS2R2[8].h * (1 - q_2r_2c_lt)), 0, 0, 0, r)
            RTS2R2_C[9] = RTS2R2[9].__copy__()
            """-------------------------------"""

            RTS2R2[10] = TDN(RTS2R2_C[9].dp(dp_suc), RTS2R2[9].h, 0, 0, 0, r)
            RTS2R2_C[10] = RTS2R2[10].__copy__()

            s3_isen_2r_2c_lp = TDN(p_mid, 0, 0, RTS2R2[10].s, 0, r)  # Isentropic compressor
            s3_isenc_2r_2c_lp = s3_isen_2r_2c_lp.__copy__()

            RTS2R2[11] = TDN(p_mid, RTS2R2_C[10].comp(eta_isen, s3_isenc_2r_2c_lp.h), 0, 0, 0, r)
            RTS2R2_C[11] = RTS2R2[11].__copy__()

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
            RTS2R2[12] = TDN(RTS2R2[2].p, (RTS2R2[2].h * q_2r_2c_ht + RTS2R2[11].h * (1 - q_2r_2c_ht)), 0, 0, 0, r)
            RTS2R2_C[12] = RTS2R2[12].__copy__()
            """-------------------------------"""

            s3_isen_2r_2c_hp = TDN(p_high, 0, 0, RTS2R2[12].s, 0, r)  # Isentropic compressor
            s3_isenc_2r_2c_hp = s3_isen_2r_2c_hp.__copy__()

            RTS2R2[13] = TDN(p_high, RTS2R2_C[12].comp(eta_isen, s3_isenc_2r_2c_hp.h), 0, 0, 0, r)
            RTS2R2_C[13] = RTS2R2[13].__copy__()

            RTS2R2[14] = TDN(RTS2R2_C[13].dp(dp_dis), RTS2R2_C[13].h, 0, 0, 0,
                             r)  # Pressure drop in discharge side entrance of compressor
            RTS2R2_C[14] = RTS2R2[14].__copy__()

            RTS2R2[15] = TDN(RTS2R2_C[14].dp(dp_cond), PSI('H', 'Q', 0, 'P', RTS2R2[14].p - dp_cond, r), 0, 0, 0,
                             r)  # H calculated by quality
            RTS2R2_C[15] = RTS2R2[15].__copy__()

            cop_2r_2c = (RTS2R2[8].h - RTS2R2[6].h) * (1 - q_2r_2c_ht)* (1 - q_2r_2c_lt) / ((RTS2R2[11].h - RTS2R2[10].h) * (1 - q_2r_2c_lt) + (RTS2R2[13].h - RTS2R2[12].h))

            """---------------------------------------------------------------------------------------------------------"""

            num_data_r_2r_2c_int = 17

            """__________________________________VCC       2 STAGE  with 2 Rec ____________________________"""

            """--------------RECIEVER ASSUMING STAEDY STATE---------- """
            """ VCC states for refrigerant """

            p_mid = (p_low + p_high) / 2

            RTS2R2I = [TDN for i in range(num_data_r_2r_2c_int)]
            RTS2R2I_C = [TDN for i in range(num_data_r_2r_2c_int)]

            RTS2R2I[0] = TDN(p_high, PSI('H', 'Q', 0, 'P', p_high, r), 0, 0, 0, r)

            RTS2R2I[0] = TDN(p_high, 0, RTS2R2I[0].dt(dt_subcool_ig), 0, 0, r)

            RTS2R2I[1] = TDN(p_mid, RTS2R2I[0].h, 0, 0, 0, r)

            q_2r_2c_ht_i = RTS2R2I[1].q()

            RTS2R2I[2] = TDN(RTS2R2I[1].p, PSI('H', 'Q', 1, 'P', RTS2R2I[1].p, r), 0, 0, 0, r)  # gas from reciever outlet

            RTS2R2I[3] = TDN(RTS2R2I[1].p, PSI('H', 'Q', 0, 'P', RTS2R2I[1].p, r), 0, 0, 0,
                             r)  # Liquid from reciever outlet

            RTS2R2I[4] = TDN(p_low + dp_evap, RTS2R2I[3].h, 0, 0, 0, r)
            RTS2R2I_C[4] = RTS2R2I[4].__copy__()

            q_2r_2c_lt_i = RTS2R2I_C[4].q()

            RTS2R2I[5] = TDN(RTS2R2I[4].p, PSI('H', 'Q', 1, 'P', RTS2R2I[4].p - dp_evap, r), 0, 0, 0,
                             r)  # gas from reciever outlet

            RTS2R2I[6] = TDN(RTS2R2I[4].p, PSI('H', 'Q', 0, 'P', RTS2R2I[4].p, r), 0, 0, 0,
                             r)  # Liquid from reciever outlet
            RTS2R2I_C[6] = RTS2R2I[6].__copy__()

            RTS2R2I[7] = TDN(RTS2R2I_C[6].dp(dp_evap), PSI('H', 'Q', 1, 'P', RTS2R2I[6].p - dp_evap, r), 0, 0, 0, r)
            RTS2R2I_C[7] = RTS2R2I[7].__copy__()

            RTS2R2I[8] = TDN(RTS2R2I[7].p, 0, RTS2R2I_C[7].dt(dt_superheat_ig), 0, 0, r)
            RTS2R2I_C[8] = RTS2R2I[8].__copy__()

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
            RTS2R2I[9] = TDN(RTS2R2I[8].p, (RTS2R2I[5].h * q_2r_2c_lt_i + RTS2R2I[8].h * (1 - q_2r_2c_lt_i)), 0, 0, 0, r)
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

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
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

            cop_2r_2c_i = (RTS2R2I[8].h - RTS2R2I[6].h) * (1 - q_2r_2c_ht_i) * (1 - q_2r_2c_lt_i) / (
                    (RTS2R2I[11].h - RTS2R2I[10].h) * (1 - q_2r_2c_lt_i) + (RTS2R2I[14].h - RTS2R2I[13].h))

            """---------------------------------------------------------------------------------------------------------"""

            num_data_r_2r_2c_i_m = 17

            """__________________________________VCC   2 STAGE  2 Rec Int Modified  ____________________________"""

            """--------------RECIEVER ASSUMING STAEDY STATE---------- """
            """ VCC states for refrigerant """

            p_mid = (p_low + p_high) / 2

            RTS2R2IM = [TDN for i in range(num_data_r_2r_2c_i_m)]
            RTS2R2IM_C = [TDN for i in range(num_data_r_2r_2c_i_m)]

            RTS2R2IM[0] = TDN(p_high, PSI('H', 'Q', 0, 'P', p_high, r), 0, 0, 0, r)

            RTS2R2IM[0] = TDN(p_high, 0, RTS2R2IM[0].dt(dt_subcool_ig), 0, 0, r)

            RTS2R2IM[1] = TDN(p_mid, RTS2R2IM[0].h, 0, 0, 0, r)

            q_2r_2c_ht_i_m = RTS2R2IM[1].q()

            RTS2R2IM[2] = TDN(RTS2R2IM[1].p, PSI('H', 'Q', 1, 'P', RTS2R2IM[1].p, r), 0, 0, 0, r)  # gas from reciever outlet

            RTS2R2IM[3] = TDN(RTS2R2IM[1].p, PSI('H', 'Q', 0, 'P', RTS2R2IM[1].p, r), 0, 0, 0,
                              r)  # Liquid from reciever outlet

            RTS2R2IM[4] = TDN(p_low + dp_evap, RTS2R2IM[3].h, 0, 0, 0, r)
            RTS2R2IM_C[4] = RTS2R2IM[4].__copy__()

            q_2r_2c_lt_i_m = RTS2R2IM_C[4].q()

            RTS2R2IM[5] = TDN(RTS2R2IM[4].p, PSI('H', 'Q', 1, 'P', RTS2R2IM[4].p - dp_evap, r), 0, 0, 0,
                              r)  # gas from reciever outlet

            RTS2R2IM[6] = TDN(RTS2R2IM[4].p, PSI('H', 'Q', 0, 'P', RTS2R2IM[4].p, r), 0, 0, 0,
                              r)  # Liquid from reciever outlet
            RTS2R2IM_C[6] = RTS2R2IM[6].__copy__()

            RTS2R2IM[7] = TDN(RTS2R2IM_C[6].dp(dp_evap), PSI('H', 'Q', 1, 'P', RTS2R2IM[6].p - dp_evap, r), 0, 0, 0, r)
            RTS2R2IM_C[7] = RTS2R2IM[7].__copy__()

            RTS2R2IM[8] = TDN(RTS2R2IM[7].p, 0, RTS2R2IM_C[7].dt(dt_superheat_ig), 0, 0, r)
            RTS2R2IM_C[8] = RTS2R2IM[8].__copy__()

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
            RTS2R2IM[9] = TDN(RTS2R2IM[8].p*(1-q_2r_2c_ht_i_m) + RTS2R2IM[2].p * q_2r_2c_ht_i_m, (RTS2R2IM[5].h * q_2r_2c_lt_i_m + RTS2R2IM[8].h * (1 - q_2r_2c_lt_i_m)*(1-q_2r_2c_ht_i_m)+ RTS2R2IM[2].h * q_2r_2c_ht_i_m), 0, 0, 0, r)
            RTS2R2IM_C[9] = RTS2R2IM[9].__copy__()
            """-------------------------------"""

            RTS2R2IM[10] = TDN(RTS2R2IM_C[9].dp(dp_suc), RTS2R2IM[9].h, 0, 0, 0, r)
            RTS2R2IM_C[10] = RTS2R2IM[10].__copy__()

            s3_isen_2r_2c_lp_i_m = TDN(p_mid, 0, 0, RTS2R2IM[10].s, 0, r)  # Isentropic compressor
            s3_isenc_2r_2c_lp_i_m = s3_isen_2r_2c_lp_i_m.__copy__()

            RTS2R2IM[11] = TDN(p_mid, RTS2R2IM_C[10].comp(eta_isen, s3_isenc_2r_2c_lp_i_m.h), 0, 0, 0, r)
            RTS2R2IM_C[11] = RTS2R2IM[11].__copy__()

            RTS2R2IM[12] = TDN(p_mid, 0, t_amb + dt_approach_air, 0, 0, r)
            RTS2R2IM_C[12] = RTS2R2IM[12].__copy__()

            """Energy balance on superheated gas from LT compressor and gas directly coming from reciever"""
            RTS2R2IM[13] = TDN(RTS2R2IM[2].p, RTS2R2IM[12].h , 0, 0, 0, r)
            RTS2R2IM_C[13] = RTS2R2IM[13].__copy__()
            """-------------------------------"""

            s3_isen_2r_2c_hp_i_m = TDN(p_high, 0, 0, RTS2R2IM[12].s, 0, r)  # Isentropic compressor
            s3_isenc_2r_2c_hp_i_m = s3_isen_2r_2c_hp_i_m.__copy__()

            RTS2R2IM[14] = TDN(p_high, RTS2R2IM_C[13].comp(eta_isen, s3_isenc_2r_2c_hp_i_m.h), 0, 0, 0, r)
            RTS2R2IM_C[14] = RTS2R2IM[14].__copy__()

            RTS2R2IM[15] = TDN(RTS2R2IM_C[14].dp(dp_dis), RTS2R2IM_C[14].h, 0, 0, 0,
                               r)  # Pressure drop in discharge side entrance of compressor
            RTS2R2IM_C[15] = RTS2R2IM[15].__copy__()

            RTS2R2IM[16] = TDN(RTS2R2IM_C[15].dp(dp_cond), PSI('H', 'Q', 0, 'P', RTS2R2IM[14].p - dp_cond, r), 0, 0, 0,
                               r)  # H calculated by quality
            RTS2R2IM_C[16] = RTS2R2IM[16].__copy__()

            cop_2r_2c_i_m = (RTS2R2IM[8].h - RTS2R2IM[6].h) * (1 - q_2r_2c_ht_i_m) * (1 - q_2r_2c_lt_i_m) / (
                    (RTS2R2IM[11].h - RTS2R2IM[10].h)  + (RTS2R2IM[14].h - RTS2R2IM[13].h))


            """___________________________________________HX FLUID STATES_______________________________________________"""

            HXTS = [TDN for i in range(num_data_hx)]
            HXTS_C = [TDN for i in range(num_data_hx)]

            HXTS[0] = TDN(p_amb, 0, CTS[3].t , 0, 0, hx)  # starting from piston expander output
            HXTS_C[0] = HXTS[0].__copy__()
            hx_in_r = HXTS[0]

            HXTS[1] = TDN(p_amb, HXTS[0].h + (m_dot_r/m_dot_hx)*(cond_in.h-cond_out.h), 0, 0, 0, hx)  # HX fluid exchanges heat with in the condensor
            HXTS_C[1] = HXTS[1].__copy__()
            hx_out_r = HXTS[1]

            HXTS[2] = TDN(p_amb, 0, t_amb+dt_approach, 0, 0,hx)  # HX fluid exchanges heat with ambient(gas cooler)
            HXTS_C[2] = HXTS[2].__copy__()
            hx_out_amb = HXTS[2]

            HXTS[3] = TDN(p_exp_max, 0, HXTS_C[2].t, 0, 0, hx)  # HXF passes through pump to get high P
            HXTS_C[3] = HXTS[3].__copy__()  # Enough to enter piston
            hx_in_exp = HXTS[3]

            HXTS[4] = TDN(exp_c_out.p, 0, CTS[3].t, 0, 0, hx) # HXF exchanges heat with piston expander
            HXTS_C[4] = HXTS[4].__copy__()
            hx_mix_exp = HXTS[4]

            HXTS[5] = TDN(CTS_C[3].p, 0, CTS[3].t, 0, 0, hx) # HXF exchanges heat with piston expander
            HXTS_C[5] = HXTS[5].__copy__()
            hx_out_exp = HXTS[5]

            """----------------------------------------------"""

            HX0 = TDN(p_amb, 0, t_amb, 0, 0, hx)

            """----------------------------------------------------------------------------------------------------------"""

            """_____________________________________________________________________________________________________
            _____________________________________________M_Dot calculation__________________________________________
            _____________________________________________________________________________________________________"""

            """_______Per kg of Cryogen_________"""
            w_exp = eta_exp_ig * (ETSI[1].p / ETSI[1].d * np.log((ETSI[1].d / ETSI[2].d)))  # ISOTHERMAL work extraction in expander
            q_hx_c = exp_c_in.h - tank_out.h  # Heat exchanged in cryogenic HX
            q_exp_c = exp_c_out.h - exp_c_in.h  # Heat exchanged in isothermal expnsion for cryogen

            """_______Per kg of Refrigerant_________"""
            w_comp = comp_out.h - comp_in.h  # Compressor work required (J/Kg refrierant)
            w_comp_rec = comp_out_rec.h - comp_in_rec.h
            q_cond_r = cond_in.h - cond_out.h  # Heat discipated at the condensor in (J/Kg refrierant)
            q_evap_r = evap_out.h - evap_in.h  # Heat absorbed from Cooling effect  (J/Kg refrierant)
            q_evap_r_rec = evap_out_rec.h - evap_in_rec.h

            """_______Per kg of HX fluid_________"""
            q_hx_r = hx_out_r.h - hx_in_r.h  # Heat exchanged for HX fluid at condensor
            q_hx_exp = hx_in_exp.h - hx_in_r.h  # Heat exchanged for HX fluid at expander
            q_hx_amb = hx_out_r.h - hx_out_amb.h

            """-------------------------------------------------------------------------------------------------"""

            """_________________________________________Matrix solution____________________________________________"""

            a = np.array([[w_exp/eta_exp_ig + q_exp_c, 0, -q_hx_exp],

                          [0, q_cond_r, -q_hx_exp],

                          [w_exp * eta_shaft * eta_aux, -w_comp, 0]])

            b = np.array([0, q_out, 0])

            try:
                x = np.linalg.solve(a, b)

            except np.linalg.LinAlgError:
                x = np.linalg.lstsq(a, b)[0]

            m_dot_c = x[0]             # Mass flow rate for cryogen in [kg/s]
            m_dot_r = x[1]             # Mass flow rate for Refrigrant in [kg/s]
            m_dot_hx = x[2]            # Mass flow rate for HX fliud in [kg/s]
            m_dot_r_g = m_dot_r * q_req
            m_dot_r_l = m_dot_r * (1-q_req)

            for xx in x:
                m_dot.append(xx)

            """----------------------------------------------------------------------------------------------------"""

            """_____________________________________Cooling effect__________________________________________________"""

            cool_prod = (q_hx_c * m_dot_c + q_evap_r * m_dot_r)   # Total cooling effect
            q_out_cal = q_out * (cool_req / cool_prod)            # Heat discipated to ambient by gas cooler
            dq = abs(q_out_cal - q_out)                           # Check for convergance
            q_out = q_out_cal                                     # Updated value of q_out

        """----------------------------------------------------------------------------------------------------------"""

        """_____________________________________________________________________________________________________________
        _____________________________________________Performance calculation____________________________________________
        _____________________________________________________________________________________________________________"""

        for xx in x:
            m_dot.append(xx)

        cop = (q_evap_r * m_dot_r) / (w_comp * m_dot_r)            # Coeficient of performance

        cop_rec = (q_evap_r_rec * m_dot_r_l) / (w_comp_rec * m_dot_r)

        qq.append(q_evap_r * m_dot_r)
        ww.append(w_comp * m_dot_r)

        cop_l.append(cop)
        cop_r_l.append(cop_rec)
        cop_r_2_l.append(cop_r_2)
        cop_2r_2c_l.append(cop_2r_2c)
        cop_2r_2c_i_l.append(cop_2r_2c_i)
        cop_2r_2c_i_m_l.append(cop_2r_2c_i_m)

        v_dot_suc = m_dot_r / comp_in.d                            # Volume flow rate at the suction side of compressor  (M3/s)
        p_ratio = comp_out.p / comp_in.p                           # Pressure ratio of outlet over inlet of compressor
        cool_kg_c = q_hx_c + w_exp * eta_aux * eta_shaft * cop     # Cooling per kg of cryogen
        cool_kg_c_l.append(cool_kg_c)
        w_vcc_sup = w_exp * eta_aux * eta_shaft * m_dot_c          # Work supplied to VCC
        w_vcc_req = w_comp * m_dot_r                               # Work required by VCC

        """-----------------------------------------------------------------------------------------------------------"""

        """_______________________________________________________________________________________________________
        _____________________________________________Expander States_____________________________________________
        __________________________________________________________________________________________________________"""
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

        """Cryogen valve opening degree.if started <0 exert negative power on piston so 0 is the best to start"""
        theta_c_open_ig = 20

        """Degree at which the cryogen valve opens"""
        theta_exh_open_ig = 180

        """HX fluid valve opening degree"""
        theta_hx_open_ig = 330

        """Assumption : Total volume of HX at top dead state is taken as th volume of Top dead state """

        """ For the whole expantion process the volume of HX taken constant equal to that of top dead space """

        v_dead_top = (m_dot_hx / hx_in_exp.d) / f_ig
        v_dead_bot_ideal = (m_dot_c / exp_exh_c.d) / f_ig
        cos_open = mt.cos((theta_c_open_ig * mt.pi) / 180)
        v_ratio_open = mt.cos((theta_c_open_ig * mt.pi) / 180)
        v_c_o_ideal = (m_dot_c / exp_c_out.d * v_ratio_open) / f_ig
        p_power = exp_c_out.p / v_ratio_open
        p_exh = p_power * (1 - cos_open) / (2)

        v_pis_c_o = r_shaft_ig * (d_bore_ig / 2) ** 2 * mt.pi * (1 - cos_open)
        v_pis_bot = r_shaft_ig * (d_bore_ig / 2) ** 2 * mt.pi * 2

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
        RPM_f = rpm

        """ Intake stroke   """
        d_bore_intake = d_bore_ig/intake_bore_ratio_ig
        gam = exp_c_in.cp()/exp_c_in.cv()
        p_star = exp_c_in.p * (2/(gam+1))**((gam)/(gam-1))
        m_dot_lim = c_disch * (d_bore_intake**2)/4*mt.pi * mt.sqrt(gam*exp_c_in.d*exp_c_in.p*((2/(gam+1)**((gam+1)/(gam-1)))))

        """--------------------------------------------------------------------------------------------------------------"""

        """_____________________________________________________________________________________________________________
          _____________________________________________EXPANDER OPTIMIZATION____________________________________________
          _____________________________________________________________________________________________________________"""

        """-------OBJECTIVE FUNCTION TO BE MINIMIZED------"""

        """------OPTIMIZED VALUES REPLACED WITH INITIAL GUESS------"""
        p1 = 3e6
        f = 20.0
        """------OPTIMIZED VALUES REPLACED WITH INITIAL GUESS------"""
        iter=0
        while round(abs(p1-p_lim_exp*0.9), -5) != 0:
            if p1-p_lim_exp*0.9>0:
                f += .2
            else:
                f -= .2

            iter+=1

            theta_c_open = 20.0
            d_bore = 0.07245951200
            l_dead = 0.001
            l_con = 0.1
            r_shaft = 0.03993971385
            intake_bore_ratio = 4.0
            theta_c_open = int(round(theta_c_open))
            """____________________________________________________________________"""

            l_str = r_shaft * 2  # length of expander chamber

            v_dead = (l_dead * mt.pi * d_bore ** 2) / 4  # Volume of dead space at the top
            v_tot = (l_str_ig * mt.pi * d_bore ** 2) / 4  # Total volume

            theta_exp_l = []  # Angle of shaft
            z_exp_l = []  # Position of piston in expander
            v_exp_l = []  # Volume of expander
            p_exp_l = []  # Pressure withing piston expander

            """ For the whole expantion process the volume of HX taken constant equal to that of top dead space """

            v_dead_top = (m_dot_hx / hx_in_exp.d) / f
            v_dead_bot_ideal = (m_dot_c / exp_exh_c.d) / f
            cos_open = mt.cos((theta_c_open * mt.pi) / 180)
            v_ratio_open = mt.cos((theta_c_open * mt.pi) / 180)
            v_c_o_ideal = (m_dot_c / exp_c_out.d * v_ratio_open) / f
            p_power = exp_c_out.p / v_ratio_open
            p_exh = p_power * (1 - cos_open) / (2)

            v_pis_c_o = r_shaft * (d_bore / 2) ** 2 * mt.pi * (1 - cos_open)
            v_pis_bot = r_shaft * (d_bore / 2) ** 2 * mt.pi * 2

            dead_top = 360  # Dead top angle
            dead_bot = 180  # Dead bottom angle

            for th in range(0, 360):

                z_c = r_shaft * (1 - mt.cos(th / 180 * mt.pi)) + l_con * (1 - mt.sqrt(1 - (r_shaft / l_con) ** 2 *
                                                                                      mt.sin(
                                                                                                th / 180 * mt.pi) ** 2))  # Piston distance from the dead top position
                v_c = (z_c + l_dead) * mt.pi * d_bore ** 2 / 4  # Total volume of piston expander based on rot degree of shaft

                """ _1_ Cryogen in__ """
                """ _2_ Power stroke_____Isothermal process___ """
                """ _3_ Return stroke__ """
                """ _4_ HX fluid in__ """
                """ _5_ Top dead center__"""

                if 0 <= th < theta_c_open:
                    p_exp_l.append(p_exp_max)

                elif theta_c_open <= th < dead_bot:
                    p_exp_l.append((ETSI[1].p / ETSI[1].d) * m_dot_c / (v_c * f))

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
            d1 = TDN(p_exp_l[theta_c_open], 0, t_exp_ig, 0, 0, c).d  # Highest density in the piston expander
            d2 = p2 / p1 * d1  # Lowest density in the piston expander
            intake = TDN(p1, 0, t_amb, 0, 0, c)
            exhaust = TDN(p2, 0, t_amb, 0, 0, c)
            w_exp_real = m_dot_c * (p1 / d1) * np.log(d1 / d2)  # Work produced based on pressures that achieved in expander
            dead_min = 1 / HXTS[3].d * f_ig * m_dot_hx  # Min dead volume required based on HX Fluid volume

        print(iter)
        f_l.append(f)
        w_out_exp.append(w_exp_real)

        # ISOTHERMAL work extraction in expander
        RPM = (1 / d2) * m_dot_c / max(v_exp_l) * 60  # Revolution per second for piston expander
        RPM_f = rpm

        """ Intake stroke   """
        d_bore_intake = d_bore / intake_bore_ratio
        gam = exp_c_in.cp() / exp_c_in.cv()
        p_star = exp_c_in.p * (2 / (gam + 1)) ** ((gam) / (gam - 1))
        m_dot_lim = c_disch * (d_bore_intake ** 2) / 4 * mt.pi * mt.sqrt(gam * exp_c_in.d * exp_c_in.p * ((2 / (gam + 1) ** ((gam + 1) / (gam - 1)))))

        w_exp_real_TDN = (TDN(p1, 0, t_amb, 0, 0, c).u() - TDN(p2, 0, t_amb, 0, 0, c).u() - TDN(p1, 0, t_amb, 0, 0, c).t * (
                    TDN(p1, 0, t_amb, 0, 0, c).s - TDN(p2, 0, t_amb, 0, 0, c).s)) * m_dot_c

        """ Exhaust stroke   """
        d_bore_exh = d_bore_ig / exhaust_bore_ratio_ig
        gam_exh = exhaust.cp()/exhaust.cv()
        p_star_exh = exhaust.p * (2/(gam_exh+1))**((gam_exh)/(gam_exh-1))
        m_dot_lim_exh = c_disch * (d_bore_exh**2)/4*mt.pi * mt.sqrt(gam_exh*exhaust.d*exhaust.p*((2/(gam_exh+1)**((gam_exh+1)/(gam_exh-1)))))
        p_drop_exh = 0.5 * (HXTS[4].d+exhaust.d)*(1-(1/exhaust_bore_ratio_ig)**4)*((m_dot_hx/HXTS[4].d)/(((180 - theta_c_open) / 360)*0.75* d_bore_exh**2 *mt.pi/4))**2 +\
                     0.5 * (exhaust.d)*(1-(1/exhaust_bore_ratio_ig)**4)*((m_dot_c/exhaust.d)/(((180 - theta_c_open) / 360)*0.75* d_bore_exh**2 *mt.pi/4))**2

        """--------------------------------------------------------------------------------------------------------------"""
        """--------------------------------------------------------------------------------------------------------------"""
        """--------------------------------------------------------------------------------------------------------------"""

        """_____________________________________________________________________________________________________________________
        _____________________________________________2nd Law Analysis____________________________________________________
        _____________________________________________________________________________________________________________________"""

        """_____________________________________________CRYOGEN_____________________________________________________________"""
        """Separation of Cryogen from HX fluid is not considered here, since two fluids are inimicible and at different phases """

        """Cryogenic heat exchanger-----Assuming constant CP"""
        s_q_c_hx = q_hx_c * (1/mlog(tank_out.t, exp_c_in.t) - 1 / t_evap) * m_dot_c
        s_tot_c_hx = - (tank_out.s - exp_c_in.s) * m_dot_c

        """Ass : for simpilicity the expansion process is devided into two different process:
        1- Heat exchange of Cryogen with HX fluid up until reaching the equilibrium TEMP
        2- Isothermal expansion of cryogen in the expander(HX fluid keeps the temp constant)
        """

        """Cryogenic HX fluid mixture in expander-----Assuming cons,tant CP"""
        s_q_exp_c = q_exp_c * (1/mlog(exp_c_in.t, exp_c_out.t) - 1 / mlog(hx_in_exp.t, hx_mix_exp.t)) * m_dot_c
        s_tot_exp_c = (exp_c_out.s - exp_c_in.s) * m_dot_c

        """Cryogenic HX fluid mixture expansion in expander-----Assuming constant temperature"""
        """Ass :
        Cryogen behaves idealy in the range p_max to p_exh
        No heat loss occurs in the expander
        """

        s_exp_c = (exp_c_out.z() + exp_exh_c.z()) / 2  * (1-eta_exp_ig) * rgas / exp_c_in.m() * mt.log(p_ratio) *m_dot_c # Isothermal expander entropy
        s_tot_exh_c = (exp_exh_c.s - exp_c_out.s) * m_dot_c

        """Expansion from lowest p reached within expander to ambient"""
        s_exp_exh = (exp_amb_c.z() + exp_exh_c.z()) / 2 * rgas / exp_c_in.m() * mt.log(p_exp_l[dead_bot - 1] / p_amb) * m_dot_c # Isothermal expander entropy
        s_tot_exp_exh = (exp_amb_c.s - exp_exh_c.s) * m_dot_c

        """_____________________________________________HX FLUID_____________________________________________________________"""

        """HX Fluid exchanging heat in condensor"""
        """ The heat exchange entropy is already calculated in refrigeration condensor"""
        s_tot_hx_cond = (hx_out_r.s - hx_in_r.s) * m_dot_hx

        """HX Fluid exchanging heat within the air-cooled heat exchanger"""
        s_hx_cond = q_out * m_dot_hx *(1/t_amb -1/mlog(hx_out_r.t, hx_out_amb.t))
        s_tot_hx_cond = (hx_out_r.s - hx_out_amb.s) * m_dot_hx

        """HX Fluid pumped to reach the pressure high enough to enter expader"""
        s_w_pump_hx = m_dot_hx / hx_in_exp.d * (hx_in_exp.p - hx_out_amb.p)/hx_in_exp.t * ((1-eta_pump)/(eta_pump))
        s_tot_pump_hx = -m_dot_hx * (hx_in_exp.s-hx_out_amb.s)

        """HX Fluid exchanging heat cryogen in expander"""
        """This consist of 2 phases of heat exchange :
        1-Heat transfer for getting cryogen from t_evap to expander temperature
        2-Heat transfer to keep expander at isothermal condition during expansion
        
        ASS: All heat transfer occurs with entropy gen as if all happened within first process
        """
        s_hx_c_exp = m_dot_hx * (hx_in_exp.h - hx_out_exp.h) * (1/mlog(exp_exh_c.t, exp_c_in.t)-1/mlog(hx_in_exp.t, hx_out_exp.t))
        s_tot_hx_c_exp = m_dot_hx * (hx_in_exp.s - hx_out_exp.s)

        """------------------------------------------------------------------------------------------------------------"""

        """_____________________________________________ Refrigeration cycle __________________________________________"""

        """Pressure drops"""
        s_dp_suc = m_dot_r * rgas/R0.m() * comp_in.z() * mt.log(1 + dp_suc/comp_in.p)
        s_dp_dis = m_dot_r * rgas/R0.m() * comp_out.z() * mt.log(1 + dp_dis/cond_in.p)
        s_dp_evap = m_dot_r * rgas/R0.m() * comp_in.z()* mt.log(1 + dp_evap/evap_out.p)
        s_dp_cond = m_dot_r * rgas/R0.m() * (cond_in.z() + cond_out.z()) / 2 * mt.log(1 + dp_cond/cond_out.p)

        """Compressor """
        s_comp_r = m_dot_r * rgas / comp_in.m() * (comp_in.z() + comp_out.z()) / 2 * ((1-eta_isen)/eta_isen) * mt.log(p_ratio)

        """Condensor heat exchange entropy gen"""
        s_cond_r = q_cond_r * m_dot_r *(1/mlog(hx_out_r.t, hx_in_r.t)-1/mlog(cond_in.t, cond_out.t))

        """TXV isentalpic throttling"""
        s_trot_r = m_dot_r * rgas / cond_out.m() * mt.log(evap_in.p / evap_out.p)

        """Evaporator heat exchange entropy"""
        s_evap_r = q_evap_r * m_dot_r *(1/mlog(evap_in.t, evap_out.t) - 1/t_reefer)


        """-----------------------------------------------------------------------------------------------------------"""
        """__Sorted List of entropy generations : __"""
        s = [s_q_c_hx, s_q_exp_c, s_hx_cond,s_exp_exh, s_w_pump_hx, s_dp_suc,  s_dp_dis, s_comp_r, s_dp_evap,
                 s_dp_cond,  s_cond_r, s_evap_r, s_trot_r]
        for value in s:
            ss.append(value)

        s_list = [['s_q_c_hx', s_q_c_hx], ['s_q_exp_c', s_q_exp_c], ['s_exp_exh', s_exp_exh], ['s_hx_cond', s_hx_cond], ['s_w_pump_hx', s_w_pump_hx], ['s_dp_suc', s_dp_suc], ['s_dp_dis', s_dp_dis], ['s_dp_evap', s_dp_evap], ['s_dp_cond', s_dp_cond], ['s_comp_r', s_comp_r], ['s_cond_r', s_cond_r], ['s_trot_r', s_trot_r], ['s_evap_r', s_evap_r]]
        sort_s = sort(s_list)
        six = ['s_q_c_hx', 's_exp_exh', 's_comp_r', 's_cond_r', 's_dp_evap', 's_q_exp_c']

        for i in range(len(six)):
            for j in range(len(s_list)):
                a = s_list[j]
                if a[0] == six[i]:
                    l.append(a[1])

        s_sum.append(sum([s_q_c_hx,s_q_exp_c,s_exp_exh,s_hx_cond,s_w_pump_hx,s_dp_suc,s_dp_dis,s_dp_evap,s_dp_cond,s_comp_r,s_cond_r,s_trot_r,s_evap_r]))
        s_tot.append((CTS[4].s - CTS[0].s) * m_dot_c)
        """------------------------------------------------------------------------------------------------------------"""

        m_l.append(m_dot_c)
        eta_2_law.append(1-(sum([s_q_c_hx,s_q_exp_c,s_exp_exh,s_hx_cond,s_w_pump_hx,s_dp_suc,s_dp_dis,s_dp_evap,s_dp_cond,s_comp_r,s_cond_r,s_trot_r,s_evap_r]))/((CTS[4].s - CTS[0].s) * m_dot_c))
        cop_o_l.append(cop)

    m_l_r.append(m_l)
    eta_l_r.append(eta_2_law)
    cop_o_l_r.append(cop_o_l)

    """______________________________________________________________________________________________________________
    ____________________________________________________PRINT________________________________________________________
    ______________________________________________________________________________________________________________"""

    # print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Refrigerant'],
    #       "\n", TDN.prop(test), "\n", TDN.prop(test))
    print("Expansion ratio of cryogenic expander = " + str(p1/p2))

    if p1 > p_lim_exp or dead_min < v_dead:
        print("!!!! WARNING!!!! :Wrong selection top dead volume(not enough space for HX or Cryogen), try larger values")

    print("Lowest pressure in the piston expander is  :  " + str(p2 / 1e5))
    print("Highest pressure in the piston expander is  : " + str(p1 / 1e5))

    print("Mass flow rate for Cryogen in [kg/s] is :     " + str(m_dot_c) + " of " + c)
    print("Mass flow rate for Refrigerant in [kg/s] is : " + str(m_dot_r) + " of " + r)
    print("Mass flow rate for HX fluid in [kg/s] is :    " + str(m_dot_hx) + " of " + hx)

    print("Work supplied to VCC is :  " + str(w_vcc_sup))
    print("Work required for VCC is :  " + str(w_vcc_req))

    print("Cooling effect (of VCC) per Kg of refrigerant is :  " + str(q_evap_r))
    print("TOTAL COOLING effect per Kg of Cryogen is :         " + str(cool_kg_c))
    print("TOTAL COOLING output :                              " + str(cool_prod))

    print("COP is :                             " + str(cop))
    print("COP with use of reciever is :        " + str(cop_rec))
    print("COP 2 stage + rec :                  " + str(cop_r_2))
    print("COP 2 stage + 2 rec :                " + str(cop_2r_2c))
    print("COP 2 stage + 2 rec + 1 int :        " + str(cop_2r_2c_i))
    print("COP 2 stage + 2 rec + 1 int + mod :  " + str(cop_2r_2c_i_m))

    print("Volume flow rate at the suction side of compressor is :  " + str(v_dot_suc))
    print("Pressure ratio in compressor is :  " + str(p_ratio))

    print("Heat rejected to ambient at gas cooler  :         " + str(q_hx_amb*m_dot_hx))
    print("Heat transfered to piston expander by HX fluid  : " + str((HXTS[3].h - HXTS[4].h) * m_dot_hx))

    print("RPM is :             " + str(RPM))
    print("W_exp_real is :      " + str(w_exp_real))
    print("W_exp_real_TDN is :  " + str(w_exp_real_TDN))

    print("Heat exchanged at condenser is  : " + str(q_cond_r * m_dot_r))
    print("Heat expander : " + str((exp_c_out.h - exp_c_in.h) * m_dot_c))

    print("Sum of entropies : " + str(sum([s_q_c_hx,s_q_exp_c,s_exp_exh,s_hx_cond,s_w_pump_hx,s_hx_c_exp,s_dp_suc,s_dp_dis,s_dp_evap,s_dp_cond,s_comp_r,s_cond_r,s_trot_r,s_evap_r])))
    print('Total work that can be exploited from cryogen : ' + str((CTS[4].ex() - CTS[0].ex()) * m_dot_c))

    print()
    print('"""--------------------------------------------"""')
    print()
    print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Refrigerant'])
    print('-------------------Single stage VCC------------')
    i = 0
    for til in RTS:
        print(i, TDN.prop(til))
        i += 1
    i = 0
    print('-------------------Single stage + receiver VCC------------')
    for tilr in RTSR:
        print(i, TDN.prop(tilr))
        i += 1
    i = 0

    print('-------------------2 stage + receiver VCC------------')
    for tilr2 in RTSR2:
        print(i, TDN.prop(tilr2))
        i += 1
    i = 0

    print('-------------------2 stage + 2 receiver VCC------------')
    for til2r2 in RTS2R2:
        print(i, TDN.prop(til2r2))
        i += 1
    i = 0

    print('-------------------2 stage + 2 receiver + 1 intercooler VCC------------')
    for til2r2i in RTS2R2I:
        print(i, TDN.prop(til2r2i))
        i += 1
    i = 0

    print('-------------------2 stage + 2 receiver + 1 intercooler modified VCC ------------')
    for til2r2im in RTS2R2IM:
        print(i, TDN.prop(til2r2im))
        i += 1
    i = 0

    print('-------------------Cryogen TDN states------------')
    print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Cryogenic Fluid'])
    for cts in CTS:
        print(i, TDN.prop(cts))
        i += 1
    i = 0
    print('-------------------HXF TDN states------------')
    print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'HX Fluid'])
    for hxt in HXTS:
        print(i, TDN.propl(hxt))
        i += 1
    i = 0

    """--------------------------------------------------------------------------------------------------------------"""
    print()
    print('"""-----------------------------------------------------------------------------------------------"""')
    print()

    """___________________________________________________________________________________________________________
    ____________________________________________________PLOTS_____________________________________________________
    ___________________________________________________________________________________________________________"""

    """__________________________________Thermodynamic states plot____________________________________________ : """

    """ 
    First plots the lines and then 'ro' spesifies red dots as the TDN states __ r ::  for red and o :: for circle
        Color               Shape                  shape
        b : blue            "8"	: octagon          "," : pixel 
        g : green           "s"	: square           "o" : circle 
        r : red             "p"	: pentagon         "v" : triangle_down 
        c : cyan            "P"	: plus (filled)    "x" : x
        m : magenta         "*"	: star             "X" : x (filled)
        y : yellow          "h"	: hexagon1         "D" : diamond
        k : black           "H"	: hexagon2         "d" : thin_diamond  
        w : white           "+"	: plus                       
    
    """

    def make_patch_spines_invisible(ax):
        ax.set_frame_on(True)
        ax.patch.set_visible(False)
        for sp in ax.spines.values():
            sp.set_visible(False)


    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()
    par3 = host.twinx()
    par4 = host.twinx()
    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.24))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(par2)

    par3.spines["right"].set_position(("axes", 1.08))
    make_patch_spines_invisible(par3)
    par4.spines["right"].set_position(("axes", 1.16))
    make_patch_spines_invisible(par4)

    # Second, show the right spine.

    par2.spines["right"].set_visible(True)

    par3.spines["right"].set_visible(True)

    par4.spines["right"].set_visible(True)

    cool_kg_c_l[:] = [x/1000 for x in cool_kg_c_l]

    p1, = host.plot(t_l, m_l_r[refrig], "bD", label= R[refrig] +" Cryogen mass flow rate [Kg/s]")
    p2, = par1.plot(t_l, cool_kg_c_l, "rd", label=R[refrig] +" Cooling per kg of cryogen [Kj/kg]")
    p3, = par2.plot(t_l, cop_o_l_r[refrig], "g*", label=R[refrig] +" VCC coeficient of performance")
    p4, = par3.plot(t_l, f_l, "c+", label=R[refrig] + " Frequency of expander")
    p5, = par4.plot(t_l, w_out_exp, "mx", label=R[refrig] + " Work produced by expander")

    host.set_xlim(0, max(t_l))
    host.set_ylim(min(min(m_l_r))*0.9, max(max(m_l_r))  *1.1)
    par1.set_ylim(min(cool_kg_c_l)*0.8, max(cool_kg_c_l)*1.08)
    par2.set_ylim(min(min(cop_o_l_r))*0.9, max(max(cop_o_l_r))*1.1)
    par3.set_ylim(min(f_l) * 0.85, max(f_l) * 1.15)
    par4.set_ylim(min(w_out_exp) * 0.8, max(w_out_exp) * 1.2)

    host.set_xlabel("Temperature variations [$^\circ$C]")
    host.set_ylabel("Cryogen mass flow rate [Kg/s]")
    par1.set_ylabel(" Cooling per kg of cryogen [Kj/kg]")
    par2.set_ylabel(" VCC coeficient of performance")
    par3.set_ylabel(" Frequency of expander [1/S]")
    par4.set_ylabel(" Work produced by expander")

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())
    par3.yaxis.label.set_color(p4.get_color())
    par4.yaxis.label.set_color(p5.get_color())

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    par3.tick_params(axis='y', colors=p4.get_color(), **tkw)
    par4.tick_params(axis='y', colors=p5.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)

    lines = [p1, p2, p3, p4, p5]

    host.legend(lines, [l.get_label() for l in lines])
    plt.title('Off-design condition of ambient temperature for ' + str(R[refrig]))

    """------------------------------------------------------"""
    h_plot = []
    p_plot = []

    h_plot_r = []
    p_plot_r = []

    h_plot_r_2 = []
    p_plot_r_2 = []

    h_plot_2r_2c = []
    p_plot_2r_2c = []

    h_plot_2r_2c_i = []
    p_plot_2r_2c_i = []

    h_plot_2r_2c_i_m = []
    p_plot_2r_2c_i_m = []

    for i in range(len(RTS)):
        h_plot.append(RTS[i].h)
        p_plot.append(RTS[i].p)

    h_plot.append(RTS[0].h)
    p_plot.append(RTS[0].p)

    for i in range(len(RTSR)):
        h_plot_r.append(RTSR[i].h)
        p_plot_r.append(RTSR[i].p)

    h_plot_r.append(RTSR[0].h)
    p_plot_r.append(RTSR[0].p)

    for i in range(len(RTSR2)):
        if i == 3:
            h_plot_r_2.append(RTSR2[9].h)
            p_plot_r_2.append(RTSR2[9].p)
        h_plot_r_2.append(RTSR2[i].h)
        p_plot_r_2.append(RTSR2[i].p)

    h_plot_r_2.append(RTSR2[0].h)
    p_plot_r_2.append(RTSR2[0].p)

    for i in range(len(RTS2R2)):
        # if i == 3:
        #     h_plot_2r_2c.append(RTS2R2[9].h)
        #     p_plot_2r_2c.append(RTS2R2[9].p)
        h_plot_2r_2c.append(RTS2R2[i].h)
        p_plot_2r_2c.append(RTS2R2[i].p)

    h_plot_2r_2c.append(RTS2R2[0].h)
    p_plot_2r_2c.append(RTS2R2[0].p)

    for i in range(len(RTS2R2I)):
        h_plot_2r_2c_i.append(RTS2R2I[i].h)
        p_plot_2r_2c_i.append(RTS2R2I[i].p)

    h_plot_2r_2c_i.append(RTS2R2I[0].h)
    p_plot_2r_2c_i.append(RTS2R2I[0].p)


    for i in range(len(RTS2R2IM)):
        h_plot_2r_2c_i_m.append(RTS2R2IM[i].h)
        p_plot_2r_2c_i_m.append(RTS2R2IM[i].p)

    h_plot_2r_2c_i_m.append(RTS2R2IM[0].h)
    p_plot_2r_2c_i_m.append(RTS2R2IM[0].p)

    plt.figure('Vapor Compressionn cycle TDN states of refrigerant ' + str(R[refrig]))
    h_p = plt.semilogy(h_plot, p_plot,'b')
    h_p_r = plt.semilogy(h_plot_r, p_plot_r,'lime')
    h_p_r_2 = plt.semilogy(h_plot_r_2, p_plot_r_2,'m')
    h_p_2r_2c = plt.semilogy(h_plot_2r_2c, p_plot_2r_2c,'c')
    h_p_2r_2c_i = plt.semilogy(h_plot_2r_2c_i, p_plot_2r_2c_i,'gold')
    h_p_2r_2c_i_m = plt.semilogy(h_plot_2r_2c_i_m, p_plot_2r_2c_i_m,'salmon')

    h_p_ro = plt.semilogy(h_plot, p_plot, 'bd')
    h_p_ro_r = plt.semilogy(h_plot_r, p_plot_r, 'lime','+')
    h_p_ro_r_2 = plt.semilogy(h_plot_r_2, p_plot_r_2, 'mx')
    h_p_ro_2r_2c = plt.semilogy(h_plot_2r_2c, p_plot_2r_2c, 'cd')
    h_p_ro_2r_2c_i = plt.semilogy(h_plot_2r_2c_i, p_plot_2r_2c_i, 'gold','+')
    h_p_ro_2r_2c_i_m = plt.semilogy(h_plot_2r_2c_i_m, p_plot_2r_2c_i_m, 'salmon','+')


    plt.legend(['1 Stage','1 Stage + 1 Reciever','2 stage + 1 Reciever ', '2 stage + 2 Reciever ', '2 stage + 2 Reciever + 1 intercooler', '2 stage + 2 Reciever + 1 intercooler modified'])

    plt.title('Vapor Compressionn cycle TDN states of refrigerant  ' + str(R[refrig]))
    plt.xlabel('Entalpy [Kj/kg]')
    plt.ylabel('Log Pressure, Log(P) [Pa]')


    """----------------------------------------------------------------------------------------------------------"""
    try:
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
    except:

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

    plt.figure('Thermodynamic states plot of cryogenic flow')

    h_p_c = plt.semilogy(h_plot_c, p_plot_c)
    h_p_ro_c = plt.semilogy(h_plot_ro_c, p_plot_ro_c, 'ro')
    plt.title('TDN states of cryogen ' + c)
    plt.xlabel('Entalpy [Kj/kg]')
    plt.ylabel('Log Pressure, Log(P) [Pa]')

    """--------------------------------------------------------------------------------------------------------------"""

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
        t_plot_hx.append(HXTS[i].t)
        h_plot_ro_hx.append(HXTS[i].h)
        t_plot_ro_hx.append(HXTS[i].t)

    plt.figure(' Thermodynamic states plot of HX fluid')
    h_p_hx = plt.plot(h_plot_hx, t_plot_hx)
    h_p_ro_hx = plt.plot(h_plot_ro_hx, t_plot_ro_hx, 'ro')
    plt.title('TDN states of HX fluid  Ethylene Glycol 60% ')
    plt.xlabel('Entalpy [Kj/kg]')
    plt.ylabel('Temperature,  [K]')

    plt.grid()  # Grid on


    """--------------------------------------------------------------------------------------------------------------"""

    """ Cryogen PV plot """

    v_plot_cc = []
    v_plot_ro_cc = []  # First plots the lines and then 'ro' spesifies red dots as the TDN states
    p_plot_cc = []
    p_plot_ro_cc = []

    """    Expander PV diagram"""

    plt.figure('Expander PV diagram')
    plt.semilogy(v_exp_l, p_exp_l)
    plt.title('Expander  P-V plot for various shaft angle ')
    plt.xlabel('Volume [m3/kg]')
    plt.ylabel('Log Pressure, Log(P)  [pa]')


    """--------------------------------------------------------------------------------------------------------------"""

# """Entropy generations for Nitrogen"""
#
# plt.figure('Entropy generations for Nitrogen')
# labels = 's_q_c_hx','s_q_exp_c','s_hx_cond','s_exp_exh','s_w_pump_hx','s_dp_suc','s_dp_dis','s_comp_r','s_dp_evap','s_dp_cond','s_cond_r','s_trot_r','s_evap_r'
# sizes = ss[0:int(len(ss)/4)]
# colors = ['orange', 'brown', 'lightcoral', 'coral', 'indianred', 'silver', 'rosybrown', 'lightcoral', 'darksalmon', 'salmon', 'grey', 'rosybrown','orangered', 'gold']
# explode = (0.15, 0.3, 0, 0, 0, 0.3, 0.15, 0, 0, 0.14, 0, 0, 0)  # explode 1st slice
#
# print(len(labels))
# print(len(sizes))
#
# plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)
#
# plt.axis('equal')
# plt.legend(bbox_to_anchor=(0.1, 0.3))

# """--------------------------------------------------------------------------------------------------------------"""
# """   Top 6 Entropy generations for various cryogens  """
# Cryogen0 =(s_sum[0], l[0], l[1], l[2], l[3], l[4], l[5])
# Cryogen1 =(s_sum[1], l[6], l[7], l[8], l[9], l[10], l[11])
# Cryogen2 =(s_sum[2], l[12], l[13], l[14], l[15], l[16], l[17])
# Cryogen3 =(s_sum[3], l[18], l[19], l[20], l[21], l[22], l[23])
#
# ind = np.arange(len(Cryogen0))  # the x locations for the groups
# width = 0.4  # the width of the bars
# fig, ax = plt.subplots()
# rects0 = ax.bar(ind - width/2, Cryogen0, width/2, label=R[0])
# rects1 = ax.bar(ind - width/4, Cryogen1, width/2, label=R[1])
# rects2 = ax.bar(ind + width/4, Cryogen2, width/2, label=R[2])
# rects3 = ax.bar(ind + width/2, Cryogen3, width/2, label=R[3])
#
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Entropy generations [J/kg/K/s]')
# ax.set_xlabel('Components causing highest losses')
# ax.set_title('Entropy generation rates for various cryogen used(Losses) in [J/kg/K/s]')
# ax.set_xticks(ind)
# ax.set_xticklabels(('Total losses', 'Cryogenic HX', 'Expander exhaust', 'Compressor', 'Condensor HX', 'Evaporator HX', 'Expander mixing HX'))
# ax.legend()
#
# """--------------------------------------------------------------------------------------------------------------"""
# """   Mass flow rate in [kg / s] for various refrigerant  """
# refrig0 = (m_dot[0], m_dot[1], m_dot[2])
# refrig1 = (m_dot[3], m_dot[4], m_dot[5])
# refrig2 = (m_dot[6], m_dot[7], m_dot[8])
# refrig3 = (m_dot[9], m_dot[10], m_dot[11])
#
# ind = np.arange(len(refrig0))  # the x locations for the groups
# width = 0.5  # the width of the bars
#
# fig, ax = plt.subplots()
#
# rects0 = ax.bar(ind - width / 2, refrig0, width / 2, label=R[0])
# rects1 = ax.bar(ind - width / 4, refrig1, width / 2, label=R[1])
# rects2 = ax.bar(ind + width / 4, refrig2, width / 2, label=R[2])
# rects3 = ax.bar(ind + width / 2, refrig3, width / 2, label=R[3])
#
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Mass flow rate in [kg / s]')
# ax.set_xlabel('Mass flow rate for various flows')
# ax.set_title('Mass flow rate in [kg / s] for ')
# ax.set_xticks(ind)
# ax.set_xticklabels(('Cryogenic', 'Refrigerant', 'HX fluid'))
# ax.legend()
#
# """   Cooling effect and work required   """
# rect0 = (qq[0], ww[0])
# rect1 = (qq[1], ww[1])
# rect2 = (qq[2], ww[2])
# rect3 = (qq[3], ww[3])
#
# ind = np.arange(len(rect0))  # the x locations for the groups
# width = 0.5  # the width of the bars
#
# fig, ax = plt.subplots()
# rects0 = ax.bar(ind - width / 2, rect0, width / 2, label=R[0])
# rects1 = ax.bar(ind - width / 4, rect1, width / 2, label=R[1])
# rects2 = ax.bar(ind + width / 4, rect2, width / 2, label=R[2])
# rects3 = ax.bar(ind + width / 2, rect3, width / 2, label=R[3])
#
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Cooling and work in [W]')
# ax.set_xlabel('Cooling and work in [W]')
# ax.set_title('Cooling effect in Evapotator and Work consumed in compressor')
# ax.set_xticks(ind)
# ax.set_xticklabels(('Cooling effect in Evapotator', 'Work consumed in compressor '))
# ax.legend()
#
# """   COP  """
# rects0 = (cop_c[0], cop_c[1], cop_c[2], cop_c[3])
# rects1 = (cop_c_r[0], cop_c_r[1], cop_c_r[2], cop_c_r[3])
# rects2 = (cop_r_2_l[0], cop_r_2_l[1], cop_r_2_l[2], cop_r_2_l[3])
# rects3 = (cop_2r_2c_l[0], cop_2r_2c_l[1], cop_2r_2c_l[2], cop_2r_2c_l[3])
# rects4 = (cop_2r_2c_i_l[0], cop_2r_2c_i_l[1], cop_2r_2c_i_l[2], cop_2r_2c_i_l[3])
# rects5 = (cop_2r_2c_i_m_l[0], cop_2r_2c_i_m_l[1], cop_2r_2c_i_m_l[2], cop_2r_2c_i_m_l[3])
#
# ind = np.arange(len(rects0))  # the x locations for the groups
# width = 0.5  # the width of the bars
#
# cop_l = [cop_c[0], cop_c_r[0], cop_r_2_l[0], cop_2r_2c_l[0], cop_2r_2c_i_l[0], cop_2r_2c_i_m_l[0]]
#
# fig, ax = plt.subplots()
# rect0 = ax.bar(ind - width / 3, rects0, width / 4, color=(1, (cop_c[0]-min(cop_l))/(max(cop_l)-min(cop_l)), 0.3), label='1 stage - simple')
# rect1 = ax.bar(ind - width / 6, rects1, width / 4, color=(1, (cop_c_r[0]-min(cop_l))/(max(cop_l)-min(cop_l)), .3 ), label='1 stage - with receiver')
# rect2 = ax.bar(ind , rects2, width / 4,            color=(1, (cop_r_2_l[0]-min(cop_l))/(max(cop_l)-min(cop_l)), 0.3), label='2 stage + rec ')
# rect3 = ax.bar(ind + width / 6, rects3, width / 4, color=(1, (cop_2r_2c_l[0]-min(cop_l))/(max(cop_l)-min(cop_l)), 0.3), label='2 stage + 2 rec ')
# rect4 = ax.bar(ind + width / 3, rects4, width / 4, color=(1, (cop_2r_2c_i_l[0]-min(cop_l))/(max(cop_l)-min(cop_l)), 0.3), label='2 stage + 2 rec + 1 int ')
# rect5 = ax.bar(ind + width / 2, rects5, width / 5, color=(1, (cop_2r_2c_i_m_l[0]-min(cop_l))/(max(cop_l)-min(cop_l)), 0.3), label='2 stage + 2 rec + 1 int + mod ')
#
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('COP(Coeficient of performance)')
# ax.set_xlabel('Various refrigerant in with and without receiver configuration')
# ax.set_title('COP changes of various refrigerant with and without receiver configuration')
# ax.set_xticks(ind)
# ax.set_xticklabels((R[0], R[1], R[2], R[3]))
# ax.legend()
#
# """--------------------------------------------------------------------------------------------------------------"""
#
# fig, axs = plt.subplots(2, 2)
#
# labels = 's_hx_cond',  's_dp_suc', 's_dp_dis', 's_comp_r', 's_dp_evap', 's_dp_cond', 's_cond_r', 's_trot_r', 's_evap_r'
# sizes0 = ss[0:int(len(ss)/4)]
# sizes1 = ss[int(len(ss)/4):int(len(ss)/4) * 2]
# sizes2 = ss[int(len(ss)/4) * 2:int(len(ss)/4)*3]
# sizes3 = ss[int(len(ss)/4) * 3:int(len(ss)/4)*4]
#
# sizes0.pop(0), sizes0.pop(0), sizes0.pop(1), sizes0.pop(1)
# sizes1.pop(0), sizes1.pop(0), sizes1.pop(1), sizes1.pop(1)
# sizes2.pop(0), sizes2.pop(0), sizes2.pop(1), sizes2.pop(1)
# sizes3.pop(0), sizes3.pop(0), sizes3.pop(1), sizes3.pop(1)
# explode = (0.3, 0, 0.3, 0, 0, 0, 0, 0.3, 0)  # explode  slice
#
# col0 = []
# for a in range(len(sizes0)):
#     col0.append((1, 1-(sizes0[a]/max(sizes0)), a/20))
#
# axs[0, 0].pie(sizes0, explode=explode, labels=labels, colors=col0, autopct='%1.1f%%', shadow=False, startangle=0)
# axs[0, 0].set_title(R[0])
# axs[0, 1].pie(sizes1, explode=explode, labels=labels, colors=col0, autopct='%1.1f%%', shadow=False, startangle=0)
# axs[0, 1].set_title(R[1])
# axs[1, 0].pie(sizes2, explode=explode, labels=labels, colors=col0, autopct='%1.1f%%', shadow=False, startangle=0)
# axs[1, 0].set_title(R[2])
# axs[1, 1].pie(sizes3, explode=explode, labels=labels, colors=col0, autopct='%1.1f%%', shadow=False, startangle=0)
# axs[1, 1].set_title(R[3])
#
# plt.legend(bbox_to_anchor=(0.5, 0.7), bbox_transform=ax.transAxes)

"""_________________________________________________PLOT SHOW________________________________________________________"""

""" PLT.SHOW() :: shows multiple figures at the same time, if the figures are intended to be shown as one closes there
should be plt.show after each individual plot, Any comment after plt.show will be ran as soon as the figure is closed"""

plt.show()

"""--------------------------------------------------------------------------------------------------------------"""







