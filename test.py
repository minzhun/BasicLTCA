import math
import numpy as np


class LoadCondition:

    def __init__(self, torque_pinion, misalignment):

        # Pinion Torque - N*m
        # Misalignment - um

        self.torque = torque_pinion
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

        # Range of Relief : 0.0 to 1.0
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

    def lead_relief(self, x):

        if self.lead_1 >= 0.0:
            temp1 = (1.0-x) * self.lead_1 / (self.lead_end-self.lead_start)
        else:
            temp1 = x * abs(self.lead_1) / (self.lead_end-self.lead_start)
        temp2 = abs(self.lead_2) / (0.5*(self.lead_end-self.lead_start))**2 * (x-0.5*(self.lead_end-self.lead_start))**2
        temp = temp1 + temp2

        # print(temp)

        return temp

    def profile_relief(self, y):

        temp1 = y * abs(self.profile_1) / (self.profile_end-self.profile_start)
        temp2 = abs(self.profile_2)/ (0.5*(self.profile_end-self.profile_start))**2 * (y-0.5*(self.profile_end-self.profile_start))**2
        if y >= self.profile_start_tip:
            temp3 = (y-self.profile_start_tip) / (self.profile_end-self.profile_start_tip) * self.profile_3
        else:
            temp3 = 0.0
        temp = temp1 + temp2 + temp3

        # print(temp)

        return temp

    def total_relief(self, x, y):

        temp1 = self.lead_relief(x)
        temp2 = self.profile_relief(y)
        temp = temp1 + temp2

        # print(temp)

        return temp


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
        # Relief
        self.relief = relief

        # Calculated in geometry()
        # Transverse Pitch - mm
        self.pt = 0.0
        # Transverse Pressure Angle - rad
        self.alpha_t_rad = 0.0
        # Reference Diameter - mm
        self.diameter_ref = 0.0
        # Transverse Base Pitch - mm
        self.pbt = 0.0
        # Base Diameter - mm
        self.diameter_base = 0.0
        # Base Helix Angle - rad
        self.beta_base_rad = 0.0

    def geometry(self):

        self.pt = math.pi * self.mn / math.cos(self.beta_rad)
        self.alpha_t_rad = math.atan(math.tan(self.alpha_n_rad)/math.cos(self.beta_rad))
        self.diameter_ref = self.z * self.mn / math.cos(self.beta_rad)
        self.pbt = self.pt * math.cos(self.alpha_t_rad)
        self.diameter_base = self.diameter_ref * math.cos(self.alpha_t_rad)
        self.beta_base_rad = math.asin(math.sin(self.beta_rad)*math.cos(self.alpha_n_rad))

        print("----------")
        print(self.pt)
        print(math.degrees(self.alpha_t_rad))
        print(self.diameter_ref)
        print(self.pbt)
        print(self.diameter_base)

    def calculate_relief(self, x, y):

        if y == -1.0:
            temp_relief = -1.0
        else:
            temp_x = x / self.b
            temp_y = (y-self.root_form) / (self.tip_form-self.root_form)
            temp_relief = self.relief.total_relief(temp_x,temp_y)

        # print(temp_relief)

        return temp_relief


class GearPair:

    def __init__(self, pinion, wheel, central_distance, profile_shift_coefficient_pinion, pinion_driving=1):

        self.a = central_distance
        self.x1 = profile_shift_coefficient_pinion
        self.pinion_driving = pinion_driving

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
        self.eap_pinion = 0.0
        self.eap_wheel = 0.0
        self.c = 0.0

        # Used and Calculated in single_stiffness()
        self.C_M = 0.8
        self.C_R = 1.0
        self.C_B = 0.975
        self.zn1 = 0.0
        self.zn2 = 0.0
        self.q_prime = 0.0
        self.c_prime = 0.0

        # Used and Calculated in calculate_te()
        self.cal_relief_pinion = pinion.calculate_relief
        self.cal_relief_wheel = wheel.calculate_relief

    def mesh_geometry(self):

        self.alpha_wt_rad = math.acos((self.z1+self.z2)*self.mn*math.cos(self.alpha_t_rad)/2./self.a/math.cos(self.beta_rad))
        self.x_sum = (self.z1+self.z2)*(math.tan(self.alpha_wt_rad)-self.alpha_wt_rad-math.tan(self.alpha_t_rad)+self.alpha_t_rad)/(2.*math.tan(self.alpha_n_rad))
        self.x2 = self.x_sum - self.x1

        self.c1 = 0.5 * self.base_diameter_pinion * math.tan(self.alpha_wt_rad)
        self.c2 = 0.5 * self.base_diameter_wheel * math.tan(self.alpha_wt_rad)
        self.c = self.c1 + self.c2
        self.sap_pinion = self.c1 + self.c2 - self.tip_form_wheel
        self.sap_wheel = self.c1 + self.c2 - self.tip_form_pinion
        self.eap_pinion = self.tip_form_pinion
        self.eap_wheel = self.tip_form_wheel

        print("----------")
        print(self.pbt)
        print(self.x1)
        print(self.x2)
        print(self.x_sum)
        print(self.sap_pinion)
        print(self.eap_pinion)
        print(self.sap_wheel)
        print(self.eap_wheel)

    def single_stiffness(self):

        self.zn1 = self.z1 / math.pow(math.cos(self.beta_base_rad), 2) / math.cos(self.beta_rad)
        self.zn2 = self.z2 / math.pow(math.cos(self.beta_base_rad), 2) / math.cos(self.beta_rad)

        self.q_prime = 0.04723 + 0.15551/self.zn1 + 0.25791/self.zn2 \
                               + (-0.00635)*self.x1 + (-0.11654)*self.x1/self.zn1 \
                               + (-0.00193)*self.x2 + (-0.24188)*self.x2/self.zn2 \
                               + 0.00529*self.x1**2 + 0.00182*self.x2**2
        self.c_prime = self.C_M * self.C_R * self.C_B * math.cos(self.beta_rad) / self.q_prime

        # TODO
        # self.c_prime = 17.0927

        print("----------")
        print(self.c_prime)
        print(self.zn1)
        print(self.zn2)

    def calculate_te(self, load_condition, n_steps, n_width, n_contact_lines=5, n_iteration=20, tol=0.001):

        np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

        te = []

        torque = load_condition.torque
        base_force = torque * 1000.0 / self.base_diameter_pinion * 2.
        print("----------")
        print(base_force)

        misalignment = load_condition.misalignment

        print("----------")

        for i in range(n_steps):

            # Define Contact Line Related

            x_pos = np.linspace(0.0, self.b_eff, n_width, endpoint=False, dtype=float)
            x_pos = x_pos + 0.5 * self.b_eff / n_width
            y_pos = np.zeros((n_contact_lines, n_width))

            # print("----------")
            # print(x_pos)

            for j in range(n_contact_lines):
                y_pos[j] = x_pos * math.tan(self.beta_base_rad) + i * self.pbt / n_steps + (j - 2) * self.pbt

            if self.pinion_driving == 1:
                y_pinion = y_pos + self.sap_pinion
                y_wheel = self.pbt - y_pos + self.c - self.sap_pinion - self.pbt
            else:
                y_pinion = y_pos + self.c - self.sap_wheel - self.pbt
                y_wheel = self.pbt - y_pos + self.sap_wheel

            # print(y_pos)

            y_pinion = np.where((y_pinion - self.sap_pinion) >= -1e-6, y_pinion, -1.)
            y_pinion = np.where((y_pinion - self.eap_pinion) <= 1e-6, y_pinion, -1.)
            y_wheel = np.where((y_wheel - self.sap_wheel) >= -1e-6, y_wheel, -1.)
            y_wheel = np.where((y_wheel - self.eap_wheel) <= 1e-6, y_wheel, -1.)

            # print("----------")
            # print(y_pinion)
            # print(y_wheel)

            relief_pinion = np.zeros((n_contact_lines, n_width))
            relief_wheel = np.zeros((n_contact_lines, n_width))
            relief_contact_line = np.zeros((n_contact_lines, n_width))

            for j in range(n_contact_lines):
                for k in range(n_width):
                    relief_pinion[j][k] = self.cal_relief_pinion(x_pos[k], y_pinion[j][k])
                    relief_wheel[j][k] = self.cal_relief_wheel(x_pos[k], y_wheel[j][k])
                    relief_contact_line[j][k] = relief_pinion[j][k] + relief_wheel[j][k]

            # print("----------")
            # print(relief_contact_line)

            misalignment_contact_line = np.zeros((n_contact_lines, n_width))

            for j in range(n_contact_lines):
                for k in range(n_width):
                    if y_pinion[j][k] == -1.:
                        misalignment_contact_line[j][k] = -1.
                    else:
                        if misalignment >= 0.0:
                            misalignment_contact_line[j][k] = x_pos[k] / self.b_eff * misalignment
                        else:
                            misalignment_contact_line[j][k] = (self.b_eff-x_pos[k]) / self.b_eff * abs(misalignment)

            # print("----------")
            # print(misalignment_contact_line)

            relief_and_misalignment_contact_line = relief_contact_line + misalignment_contact_line
            relief_and_misalignment_contact_line_1 = relief_and_misalignment_contact_line[relief_and_misalignment_contact_line > -1e-6]

            if len(relief_and_misalignment_contact_line_1) == 0:
                relief_and_misalignment_contact_line_1 = y_pinion[y_pinion > 0.0]
                relief_and_misalignment_contact_line_1[:] = 0.0

            # print(relief_and_misalignment_contact_line_1)

            temp_te = base_force / self.b_eff / self.c_prime
            converge_or_not = 0
            ite_number = 0

            for ite in range(n_iteration):

                deflection_contact_line_1 = temp_te - relief_and_misalignment_contact_line_1
                deflection_contact_line_1 = deflection_contact_line_1[deflection_contact_line_1 > 0.0]
                temp_deflection = np.sum(deflection_contact_line_1) * math.cos(self.beta_base_rad)
                diff = base_force - temp_deflection * self.b_eff / n_width * self.c_prime

                # print("te  ", temp_te)
                # print("deflection  ", temp_deflection)
                # print("diff  ", diff)
                # print(deflection_contact_line_1)
                # print("\n")

                if abs(diff) < base_force * tol:
                    te.append(temp_te)
                    converge_or_not = 1
                    ite_number = ite
                    break
                else:
                    if diff < 0.0:
                        delta = min(abs(temp_te*0.5), abs(diff/self.b_eff/self.c_prime))
                        delta = -1.*delta
                    else:
                        delta = min(abs(temp_te*2.), abs(diff/self.b_eff/self.c_prime))
                    temp_te = temp_te + delta

            if converge_or_not == 1:
                print("step ", i, " converged at ", ite_number)
            else:
                print("step ", i, "not converged")

        te = np.array(te)
        te_peak = np.max(te) - np.min(te)

        print(te)
        print(len(te))
        print(te_peak)


load_condition_1 = LoadCondition(500.0, 0.0)
#
relief_1 = Relief(0.0, 0.0, 0.0, 0.0,0.0)
relief_2 = Relief(10.0, 0.0, 10.0, 0.0, 0.0)
#
pinion_1 = Gear(31, 2.271, 20.0, 25.0, 20.0, 19.775, 8.06, relief_1)
wheel_1 = Gear(41, 2.271, 20.0, 25.0, 20.0, 24.142, 12.175, relief_1)
pinion_1.geometry()
wheel_1.geometry()
#
gear_pair_1 = GearPair(pinion_1, wheel_1, 90.0, 0.0)
gear_pair_1.mesh_geometry()
gear_pair_1.single_stiffness()
gear_pair_1.calculate_te(load_condition_1, 32, 40)
