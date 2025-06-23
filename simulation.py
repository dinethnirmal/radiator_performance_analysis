import math as m
import matplotlib.pyplot as plt
import numpy as np



DRZ400 = {
    'name': 'DRZ400',
    'K_f':237, # fin heat transfer coefficient
    'N_f': 633, #633,
    'N_p': 12, #12
    'N_ct': 11, #11
    'N_r': 2,
    'B_w': 118.1e-3, #118
    'B_h': 225.3e-3, #225
    'T_s':2.5e-3,
    'B_t': 40.1e-3,
    'F_p': 1.58e-3,
    'F_t': 0.2e-3, #0.1e-3
    'alpha_f': 0,
    'R_f': 0.79e-3,
    'F_h': 7e-3,
    'Y_cl': 14.5e-3,
    'Y_cw': 1.5e-3,
    'Y_t': 0.3e-3,
    'Y_l': 230.3e-3,
    'L_p': 3e-3, #1e-3,
    'L_l': 6e-3,
    'L_h': 0.4e-3,
    'L_a': 35,
}

specimen_04 = {
    'name': 'test_04',
    'K_f':237, # fin heat transfer coefficient
    'N_f': 633, #633,
    'N_p': 12, #12
    'N_ct': 11, #11
    'N_r': 2,
    'B_w': 118.1e-3, #118
    'B_h': 225.3e-3, #225
    'T_s':2.5e-3,
    'B_t': 36.6e-3,
    'F_p': 2.25e-3,
    'F_t': 0.2e-3, #0.1e-3
    'alpha_f': 0,
    'R_f': 0.79e-3,
    'F_h': 8e-3,
    'Y_cl': 14.5e-3,
    'Y_cw': 1.5e-3,
    'Y_t': 0.3e-3,
    'Y_l': 230.3e-3,
    'L_p': 1.2e-3, #1e-3,
    'L_l': 6e-3,
    'L_h': 0.4e-3,
    'L_a': 28,
}

YFM700R = {
    'name': 'YFM700R',
    'K_f':237, # fin heat transfer coefficient
    'N_f': 633,
    'N_p': 19,
    'N_ct': 18,
    'N_r': 2,
    'B_w': 186e-3,
    'B_h':304e-3,
    'T_s':2.5e-3,
    'B_t': 40.1e-3,
    'F_p': 1.58e-3,
    'F_t': 0.2e-3, #0.1e-3
    'alpha_f': 0,
    'R_f': 0.79e-3,
    'F_h': 7e-3,
    'Y_cl': 14.5e-3,
    'Y_cw': 1.5e-3,
    'Y_t': 0.3e-3,
    'Y_l': 310e-3,
    'L_p': 3e-3, #1e-3,
    'L_l': 6e-3,
    'L_h': 0.4e-3,
    'L_a': 35,
}

# at 30
air = {
    'rho':1.16,#1.16, # dry air at 1 atm
    'mu': 1.9e-5, #pas
    'c_pa': 1007,  #J/m2K
    'v': 20, # m/s
    'Pr': 0.7,
    'k': 0.0270,
    'T_i_a': 33,
}

# at 50C
water = {
    'rho': 1000, #kg/m-3
    'mu': 0.0005474, #Pa.s
    'Pr': 4,
    'm_dot': 0.2, #kg/s
    'c_pc': 4180, #j/kgK
    'k': 0.64060,  #W/mK
    'T_o_c': 45,
    'T_i_c': 55,
}

# CONSTANTS



class Radiator:

    def __init__(self, data):
        self.name = data['name']
        self.K_f = data['K_f']
        self.N_f = data['N_f']
        self.N_p = data['N_p']
        self.N_ct = data['N_ct']
        self.N_r = data['N_r']
        self.B_h = data['B_h']
        self.B_w = data['B_w']
        self.B_t = data['B_t']
        self.F_p = data['F_p']
        self.F_t = data['F_t']
        self.F_h = data['F_h']
        # self.T_s = data['T_s']
        self.alpha_f = data['alpha_f']
        self.R_f = data['R_f']
        self.Y_cl = data['Y_cl']
        self.Y_cw = data['Y_cw']
        self.Y_t = data['Y_t']
        self.Y_l = data['Y_l']
        self.R_t = self.Y_cw/2
        self.L_p = data['L_p']
        self.L_l = data['L_l']
        self.L_a = data['L_a']
        self.L_h = self.L_p * m.sin(self.L_a*m.pi/180) # louver height
        self.T_d = self.Y_cl * self.N_r # tube depth
        self.T_p = self.F_h + self.Y_cw # tube pitch
        # self.L_d = 36.6e-3

    def calculate_radiator_parameters(self):
        # F_l = total length of a fin
        # A_fr_r = radiator frontal area
        # A_fr_t = frontal area of tubes
        # A_fr_f = frontal area of fins
        # A_f = total fin area
        # A_a = air side heat transfer area
        # A_c = coolant side heat transfer area
        # A_pa = air flow minimum area

        self.F_l = m.pi*self.R_f + (self.F_h - 2*self.R_f)/m.cos(self.alpha_f*m.pi/180)
        self.A_fr_r = self.B_w * self.B_h
        self.A_fr_t = self.Y_l * self.Y_cw * self.N_ct
        self.A_fr_f = self.F_t * self.F_l * self.N_f * self.Y_l * self.N_p
        self.A_f = 2 * self.B_t * self.F_l * self.N_f * self.Y_l * self.N_p
        self.A_a = self.A_f + 2 * self.N_ct * self.Y_l * self.N_r * ((self.Y_cl - 2*self.R_t) + 2*m.pi*self.R_t)
        self.A_c = (2*m.pi*(self.R_t - self.Y_t) + 2*(self.Y_cl - 2*self.R_t))*self.Y_l*self.N_r*self.N_ct
        self.A_pa = self.A_fr_r - self.A_fr_f - self.A_fr_t


    def set_U(self, air, coolant):

    # AIR SIDE CALCULATIONS

        rho_a = air['rho']
        c_p_a = air['c_pa']
        mu_a = air['mu']
        v_a = air['v']
        Pr_a = air['Pr']
        k_a = air['k']



        G = (self.A_fr_r*rho_a*v_a)/ self.A_pa
        Re_lp = G*self.L_p/mu_a

        # fan

        # Using Davenport's correlation 1983
        # J1 = 0.249 * Re_lp**(-0.42) * self.L_h**0.33 * (self.L_h/self.F_h)**1.1 * (self.F_h**0.26)

        # Using Dong, Chen correlation 2007
        J = 0.26712*Re_lp**(-0.1944) * (self.L_a/90)**0.257 * (self.F_p/self.L_p)**(-0.5177) * (self.F_h/self.L_p)**(-1.9045) * (self.L_l/self.L_p)**(1.7159) * (self.B_t/self.L_p)**(-0.2147) * (self.F_t/self.L_p)**(-0.05)
        f = 0.54486*Re_lp**(-0.3068) * (self.L_a/90)**0.444 * (self.F_p/self.L_p)**(-0.9925) * (self.F_h/self.L_p)**0.5458 * (self.L_l/self.L_p)**(-0.2003) * (self.B_t/self.L_p)**(0.0688)
        self.J = J
        self.Re_lp = Re_lp

        # using Chang's 1997 correlation
        # J3 = Re_lp**(-0.49) * (self.L_a/90)**0.27 *  (self.F_p/self.L_p)**(-0.14) * (self.F_l/self.L_p)**(-0.29) * (self.T_d/self.L_p)**(-0.23) * (self.L_l/self.L_p)**(0.68) * (self.T_p/self.L_p)**(-0.28) * (self.F_t/self.L_p)**(-0.05)

        h_a = J * G * c_p_a / Pr_a**(2/3)
        self.h_air = h_a

        l = self.F_h/2
        m_eff = m.sqrt(2*h_a/(self.K_f*self.F_t))
        # print('meff', m_eff)
        eff_f = m.tanh(m_eff*l)/(m_eff*l)
        eff_o_a = 1 - (1-eff_f)*self.A_f/self.A_a

        print('air mass flowrate: ', v_a*self.A_fr_r*rho_a)
        print('n0 =', eff_o_a)
        print('=' * 20)
        print('G:', G)
        print("Re_Lp:", Re_lp)
        print('J:', J)
        print('f:', f)
        print('h_a:', h_a)
        print('=' * 20)

    # coolant side calculations

        rho_c = coolant['rho']
        mu_c = coolant['mu']
        Pr_c = coolant['Pr']
        m_dot_c = coolant['m_dot']
        k_c = coolant['k']
        m_dot_c_t = m_dot_c/(self.N_ct*self.N_r)


        A_ff_t = (self.Y_cl - self.Y_cw)*(self.Y_cw - 2*self.Y_t) + m.pi*((self.Y_cl - 2*self.Y_t)/2)**2 # freeflow area per tube
        P_ff_t = 2*(self.Y_cl - self.Y_cw) + m.pi*(self.Y_cw - self.Y_t*2)
        D_h_t = 4*A_ff_t/P_ff_t

        G_c = m_dot_c_t/A_ff_t
        Re_c = G_c * D_h_t/mu_c

        # Hansen equations --> mean nusselt number
        Nu_c = 3.66 + (0.0668 * (D_h_t/self.Y_l) * Re_c * Pr_c)/(1 + 0.04*((D_h_t/self.Y_l)*Re_c*Pr_c)**(2/3))
        # Nu_c = 0.023*Re_c**0.8*Pr_c**0.4


        h_c = k_c * Nu_c/D_h_t

        print('G_c:', G_c)
        print('D_h:', D_h_t)
        print('Re_c:', Re_c)
        print('Nu_c:', Nu_c)
        # print('Nu_c1:' , Nu_c1)
        print('h_c:', h_c)
        print('='*20)
    # overall heat tranfer coefficent
        R = 1/(eff_o_a*h_a*self.A_a) + (1/(h_c*self.A_c))
        self.UA = 1/R
        print('UA --->', self.UA)
        # print('U -->', self.UA/self.A_fr_r)


    # coolant pressure drop across red.

        # for laminar flow in rectangular pipe l/w = 8

        f1 = (0.79*m.log(Re_c, m.e) - 1.64)**(-2)
        f2 = 83/Re_c;
        print(f2)
        V_c = m_dot_c_t/(rho_c*A_ff_t)
        print(V_c)
        h_drop_c = rho_c*f1*self.Y_l*V_c**2 / (D_h_t*2)
        print(h_drop_c)

    def NTU(self, air, coolant):
        rho_a = air['rho']
        c_p_a = air['c_pa']
        mu_a = air['mu']
        v_a = air['v']
        Pr_a = air['Pr']
        k_a = air['k']
        T_i_a = air['T_i_a']
        # T_o_a = air['T_o_a']

        rho_c = coolant['rho']
        mu_c = coolant['mu']
        Pr_c = coolant['Pr']
        m_dot_c = coolant['m_dot']
        k_c = coolant['k']
        c_p_c = coolant['c_pc']
        T_i_c = coolant['T_i_c']


        m_dot_a = self.A_fr_r * v_a * rho_a
        C_a = m_dot_a * c_p_a
        C_c = m_dot_c * c_p_c

        C_min = min(C_a, C_c)
        C_max = max(C_a, C_c)
        q_max = C_min * (T_i_c - T_i_a)

        NTU = self.UA/C_min
        Cr = C_min/C_max

        eff = 1 - m.exp((NTU**0.22)*(m.exp(-Cr * NTU**0.78) - 1)/Cr)
        self.q = q_max*eff
        print('eff',eff)
        print('q_max', q_max)
        print('q', self.q)



#config 1
def config1(water, rad1):

    motor_heat = 1500
    inv_heat = 1000
    T0 = water['T_o_c']
    T1 = T0 + motor_heat/(water['m_dot'] * water['c_pc'])
    T3 = T1 + inv_heat/(water['m_dot'] * water['c_pc'])

    water['T_i_c'] = T3

    r1 = Radiator(rad1)
    r1.calculate_radiator_parameters()
    r1.set_U(air, water)
    r1.NTU(air, water)
    h = r1.h_air
    J = r1.J
    Re_lp = r1.Re_lp
    T4 = T3 - r1.q/(water['m_dot'] * water['c_pc'])

    print('T0', T0)
    print('T1', T1)
    print('T3', T3)
    print('T4', T4)

    return [r1.q, T4, h, J, Re_lp]



def config2(water, rad1, rad2):

    motor_heat = 1500
    inv_heat = 1000
    T0 = water['T_o_c']
    T1 = T0 + motor_heat/(water['m_dot'] * water['c_pc'])
    T3 = T1 + inv_heat/(water['m_dot'] * water['c_pc'])

    water['T_i_c'] = T3

    r1 = Radiator(rad1)
    r1.calculate_radiator_parameters()
    r1.set_U(air, water)
    r1.NTU(air, water)
    h = r1.h_air
    J = r1.J
    Re_lp = r1.Re_lp
    T4 = T3 - r1.q/(water['m_dot'] * water['c_pc'])

    water['T_i_c'] = T4
    r2 = Radiator(rad2)
    r2.calculate_radiator_parameters()
    r2.set_U(air, water)
    r2.NTU(air, water)


    T5 = T4 - r2.q/(water['m_dot'] * water['c_pc'])

    # print("++++++++++++++++++++++++++++++++++")
    # print('T0', T0)
    # print('T1', T1)
    # print('T3', T3)
    # print('T4', T4)
    # print('T5', T5)
    #
    # print('q1:', r1.q)
    # print('q2:', r2.q)
    return [r1.q, r2.q, T5, h, J, Re_lp]

def steady_state_run_config1(rad):
    water = {
        'rho': 1000,  # kg/m-3
        'mu': 0.0005474,  # Pa.s
        'Pr': 5,
        'm_dot': 0.2,  # kg/s
        'c_pc': 4180,  # j/kgK
        'k': 0.64060,  # W/mK
        'T_o_c': 45,
        'T_i_c': 55,
        # 'm_dot_fan':
    }
    temps = []

    for i in range(100):
        q, T4, h, J, Re = config1(water, rad)
        water['T_o_c'] = T4
        temps.append(T4)
    else:
        print("steady state heat tranfer:", q)

    plt.plot(temps, label=f"L: {rad['name']}")
    plt.title("coolent radiator outlet temperature variation with each cycle")
    plt.xlabel("Cycle number")
    plt.ylabel("temperature (C)")
    plt.legend()
    plt.grid()
    plt.show()

def steady_state_run(rad1, rad2):
    water = {
        'rho': 1000,  # kg/m-3
        'mu': 0.0005474,  # Pa.s
        'Pr': 5,
        'm_dot': 0.2,  # kg/s
        'c_pc': 4180,  # j/kgK
        'k': 0.64060,  # W/mK
        'T_o_c': 45,
        'T_i_c': 55,
        # 'm_dot_fan':
    }
    temps = []

    for i in range(50):
        q1, q2, T5, h, J, Re = config2(water, rad1, rad2)
        water['T_o_c'] = T5
        temps.append(T5)
    else:
        print("steady state heat tranfer:", q1 + q2)

    plt.plot(temps, label=f"L: {rad1['name']}, R: {rad2['name']}")
    plt.title("System temperature variation with each cycle")
    plt.xlabel("Cycle number")
    plt.ylabel("temperature (C)")
    plt.legend()
    plt.grid()
    plt.show()




def dissipation_rate_with_velocity(water, rad1, rad2, max_vel):

    vs = list(range(1,max_vel))

    q1s = []
    q2s = []
    qT = []
    for v in range(1,max_vel):
        air['v'] = v
        q1, q2, T5, h, J, Re = config2(water, rad1, rad2)
        q1s.append(q1)
        q2s.append(q2)
        qT.append(q1+q2)

    plt.plot(vs, q1s, label=f"left: {rad1['name']}")
    plt.plot(vs, q2s, label=f"right: {rad2['name']}")
    plt.plot(vs, qT, 'red', label='both')
    plt.title('heat dissapation (W) vs velocity(m/s)')
    plt.xlabel('velocity (m/s)')

    plt.ylabel('q (W)')
    plt.legend()
    plt.grid()
    plt.show()


def test1(water, rad, max_vel):

    vs = list(range(1, max_vel))
    q1s = []

    for v in range(1,max_vel):
        air['v'] = v
        q1, q2, T5, h, J, Re = config2(water, rad, rad)
        q1s.append(h)


    plt.scatter(vs, q1s, label=f"left: {rad['name']}")

    plt.title('h vs velocity(m/s)')
    plt.xlabel('velocity (m/s)')

    plt.ylabel('h_air')
    plt.legend()
    plt.grid()
    plt.show()

def J_vs_Re(water, rad, max_vel):

    vs = list(range(1, max_vel))
    # q1s = []
    Js = []
    Res = []

    for v in range(1,max_vel):
        air['v'] = v
        q1, q2, T5, h, J, Re = config2(water, rad, rad)
        # q1s.append(h)
        Js.append(J)
        Res.append(Re)


    plt.scatter(Res, Js, label=f"left: {rad['name']}")

    plt.title('J vs Re_lp')
    plt.xlabel('Re_lp')

    plt.ylabel('J')
    plt.legend()
    plt.grid()
    plt.show()


# dissipation_rate_with_velocity(water, DRZ400, DRZ400, 30)
steady_state_run(DRZ400, DRZ400)
# steady_state_run_config1(DRZ400)

# r1 = Radiator(DRZ400)
# r1.calculate_radiator_parameters()
# r1.set_U(air, water)
# r1.NTU(air, water)
# print(r1.q)


