import math
import numpy as np


class LoadCondition:

    def __init__(self, torque, misalignment):

        # Torque - N*m
        # Misalignment - um

        self.torque = torque
        self.misalignment = misalignment

        print("----------")
        print(self.torque)
        print(self.misalignment)


class Relief:

    def __init__(self, lead_linear, lead_crowning, profile_linear, profile_barrelling, profile_tip):

        # Relief - um
        self.lead_1 = lead_linear
        self.lead_2 = lead_crowning
        self.profile_1 = profile_linear
        self.profile_2 = profile_barrelling
        self.profile_3 = profile_tip

        # range of relief : 0.0 to 1.0
        self.lead_start = 0.0
        self.lead_end = 1.0
        self.profile_start = 0.0
        self.profile_end = 1.0
        self.profile_start_tip = 0.95

    def lead_start_end(self, lead_start, lead_end):

        self.lead_start = lead_start
        self.lead_end = lead_end

    def profile_start_end(self, profile_start, profile_end):

        self.profile_start = profile_start
        self.profile_end = profile_end

    def profile_tip_start(self, profile_start_tip):

        self.profile_start_tip = profile_start_tip

    def lead_relief(self):

        pass

    def profile_relief(self):

        pass

    def total_relief(self):

        pass


class Gear:

    def __init__(self, z, mn, alpha_n, beta, b, tip_form, root_form, relief):

        # Number of Teeth
        self.z = z
        # Normal Module - mm
        self.mn = mn
        # Normal Pressure Angle - deg
        self.alpha_n = alpha_n
        self.alpha_n_rad = math.radians(alpha_n)
        # Helix Angle - deg
        self.beta = beta
        self.beta_rad = math.radians(beta)
        # Face Width - mm
        self.b = b
        # Tip Form and Root Form using Rolling Distance - mm
        self.tip_form = tip_form
        self.root_form = root_form

        # Calculated in geometry()
        self.pt = 0.0
        self.alpha_t_rad = 0.0
        # Transverse Base Pitch - mm
        self.pbt = 0.0
        # Base Diameter - mm
        self.diameter_base = 0.0
        # Reference Diameter - mm
        self.diameter_ref = 0.0
        # Base Helix Angle
        self.beta_base_rad = 0.0

    def geometry(self):

        self.pt = math.pi * self.mn / math.cos(self.beta_rad)
        self.alpha_t_rad = math.atan(math.tan(self.alpha_n_rad)/math.cos(self.beta_rad))
        self.pbt = self.pt * math.cos(self.alpha_t_rad)
        self.diameter_ref = self.z * self.mn / math.cos(self.beta_rad)
        self.diameter_base = self.diameter_ref * math.cos(self.alpha_t_rad)
        self.beta_base_rad = math.asin(math.sin(self.beta_rad)*math.cos(self.alpha_n_rad))

        # print("----------")
        # print(self.pt)
        # print(math.degrees(self.alpha_t_rad))
        # print(self.diameter_ref)
        # print(self.pbt)
        # print(self.diameter_base)

    def calculate_relief(self, coord_in_active_plane):

        pass


class GearPair:

    def __init__(self, pinion, wheel, central_distance, profile_shift_coefficient_pinion):

        self.a = central_distance
        self.x1 = profile_shift_coefficient_pinion

        self.z1 = pinion.z
        self.z2 = wheel.z
        self.mn = pinion.mn
        self.alpha_n_rad = pinion.alpha_n_rad
        self.beta_rad = pinion.beta_rad
        self.b_eff = min(pinion.b, wheel.b)
        self.alpha_t_rad = pinion.alpha_t_rad
        self.base_diameter_pinion = pinion.diameter_base
        self.base_diameter_wheel = wheel.diameter_base
        self.tip_form_pinion = pinion.tip_form
        self.tip_form_wheel = wheel.tip_form
        self.beta_base_rad = pinion.beta_base_rad
        self.pbt = pinion.pbt

        # Calculated in mesh_geometry()
        self.alpha_wt_rad = 0.0
        self.x_sum = 0.0
        self.x2 = 0.0
        self.c1 = 0.0
        self.c2 = 0.0
        self.sap_pinion = 0.0
        self.sap_wheel = 0.0

        # Used and Calculated in single_stiffness()
        self.C_M = 0.8
        self.C_R = 1.0
        self.C_B = 0.975
        self.zn1 = 0.0
        self.zn2 = 0.0
        self.q_prime = 0.0
        self.c_prime = 0.0

        # Used and Calculated in calculate_te()
        self.te = []
        # self.pos_x = np.zeros(1)
        # self.pos_y = np.zeros(1)

    def mesh_geometry(self):

        self.alpha_wt_rad = math.acos((self.z1+self.z2)*self.mn*math.cos(self.alpha_t_rad)/2./self.a/math.cos(self.beta_rad))
        self.x_sum = (self.z1+self.z2)*(math.tan(self.alpha_wt_rad)-self.alpha_wt_rad-math.tan(self.alpha_t_rad)+self.alpha_t_rad)/(2.*math.tan(self.alpha_n_rad))
        self.x2 = self.x_sum - self.x1

        self.c1 = 0.5 * self.base_diameter_pinion * math.tan(self.alpha_wt_rad)
        self.c2 = 0.5 * self.base_diameter_wheel * math.tan(self.alpha_wt_rad)
        self.sap_pinion = self.c1 + self.c2 - self.tip_form_pinion
        self.sap_wheel = self.c1 + self.c2 - self.tip_form_wheel

        print("----------")
        # print(self.x_pinion)
        # print(self.x_wheel)
        print(self.x_sum)
        print(self.sap_pinion)
        print(self.sap_wheel)

    def single_stiffness(self):

        self.zn1 = self.z1 / math.pow(math.cos(self.beta_base_rad), 2) / math.cos(self.beta_rad)
        self.zn2 = self.z2 / math.pow(math.cos(self.beta_base_rad), 2) / math.cos(self.beta_rad)

        self.q_prime = 0.04723 + 0.15551/self.zn1 + 0.25791/self.zn2 \
                               + (-0.00635)*self.x1 + (-0.11654)*self.x1/self.zn1 \
                               + (-0.00193)*self.x2 + (-0.24188)*self.x2/self.zn2 \
                               + 0.00529*self.x1**2 + 0.00182*self.x2**2
        self.c_prime = self.C_M * self.C_R * self.C_B * math.cos(self.beta_rad) / self.q_prime

        print("----------")
        print(self.c_prime)
        print(self.zn1)
        print(self.zn2)

    def calculate_te(self, load_condition, n_steps, n_width, n_contact_lines=5, n_iteration=20, tol=0.01):

        np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

        self.misalignment = load_condition.misalignment
        self.pos_x = np.linspace(0.0, self.b_eff, n_width+1, endpoint=True, dtype=float)
        print(self.pos_x)
        self.pos_y = np.zeros((n_contact_lines, n_width+1))
        print(self.pos_y)

        for i in range(n_steps):

            for j in range(n_contact_lines):

                self.pos_y[j] = self.pos_x * math.tan(self.beta_base_rad) + i * self.pbt / n_steps + (j - 2) * self.pbt

            print(self.pos_y)


            for ite in range(n_iteration):

                pass


load_condition_1 = LoadCondition(100, 10)
relief_1 = Relief(0.0, 0.0, 0.0, 0.0,0.0)
pinion_1 = Gear(31, 2.271, 20.0, 25.0, 20.0, 19.775, 8.06, relief_1)
pinion_1.geometry()
wheel_1 = Gear(41, 2.271, 20.0, 25.0, 20.0, 24.142, 12.175, relief_1)
wheel_1.geometry()
gear_pair_1 = GearPair(pinion_1, wheel_1, 90.0, 0.0)
gear_pair_1.mesh_geometry()
gear_pair_1.single_stiffness()
gear_pair_1.calculate_te(load_condition_1, 1, 10)