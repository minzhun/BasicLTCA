import math
import numpy as np
from matplotlib import pyplot as plt


class Load:

    def __init__(self, torque_pinion, misalignment):

        # Pinion Torque - N*m
        # Misalignment - um

        self.torque = torque_pinion
        self.misalignment = misalignment

        # print("----- Load -----")
        # print(self.torque)
        # print(self.misalignment)


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

        # print("----- Relief -----")
        # print(self.lead_1)
        # print(self.lead_2)
        # print(self.profile_1)
        # print(self.profile_2)
        # print(self.profile_3)

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
        self.total_relief = relief.total_relief

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

        # print("----- Gear Geometry -----")
        # print(self.pt)
        # print(math.degrees(self.alpha_t_rad))
        # print(self.diameter_ref)
        # print(self.pbt)
        # print(self.diameter_base)

    def calculate_relief(self, x, y):

        if y < 0.0:
            temp_relief = -1.0
        else:
            temp_x = x / self.b
            temp_y = (y-self.root_form) / (self.tip_form-self.root_form)
            temp_relief = self.total_relief(temp_x, temp_y)

        # print(temp_relief)

        return temp_relief


class GearPair:

    def __init__(self, pinion, wheel, central_distance, profile_shift_coefficient_pinion, active_flank="right", pinion_driving=1):

        # Input
        # Centre Distance - mm
        self.a = central_distance
        # Profile Shift Coefficient of Pinion
        self.x1 = profile_shift_coefficient_pinion
        # If Pinion is driving or not
        self.pinion_driving = pinion_driving
        # Active Flank
        self.active = active_flank

        # Property of Gear Pair
        self.z1 = pinion.z
        self.z2 = wheel.z
        self.mn = pinion.mn
        self.alpha_n_rad = pinion.alpha_n_rad
        self.alpha_t_rad = pinion.alpha_t_rad
        self.beta_rad = pinion.beta_rad
        self.beta_base_rad = pinion.beta_base_rad
        self.b_eff = min(pinion.b, wheel.b)
        self.b_pinion = pinion.b
        self.b_wheel = wheel.b
        self.base_diameter_pinion = pinion.diameter_base
        self.base_diameter_wheel = wheel.diameter_base
        self.tip_form_pinion = pinion.tip_form
        self.tip_form_wheel = wheel.tip_form
        self.pbt = pinion.pbt

        # Calculated in mesh_geometry()
        self.alpha_wt_rad = 0.0
        self.x_sum = 0.0
        self.x2 = 0.0
        self.c1 = 0.0
        self.c2 = 0.0
        self.c = 0.0
        self.sap_pinion = 0.0
        self.sap_wheel = 0.0
        self.eap_pinion = 0.0
        self.eap_wheel = 0.0
        self.transverse_contact_ratio = 0.0

        # Used and Calculated in mesh_stiffness()
        self.C_M = 0.8
        self.C_R = 1.0
        self.C_B = 0.975
        self.zn1 = 0.0
        self.zn2 = 0.0
        self.q_prime = 0.0
        self.c_prime = 0.0
        self.c_gamma = 0.0

        # Used in calculate_te()
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
        self.transverse_contact_ratio = (self.eap_pinion-self.sap_pinion) / self.pbt

        print("----- Mesh Geometry -----")
        print(self.x1)
        print(self.x2)
        print(self.x_sum)
        print(self.sap_pinion)
        print(self.eap_pinion)
        print(self.sap_wheel)
        print(self.eap_wheel)
        print(self.transverse_contact_ratio)

    def mesh_stiffness(self):

        self.zn1 = self.z1 / math.pow(math.cos(self.beta_base_rad), 2) / math.cos(self.beta_rad)
        self.zn2 = self.z2 / math.pow(math.cos(self.beta_base_rad), 2) / math.cos(self.beta_rad)
        self.q_prime = 0.04723 + 0.15551/self.zn1 + 0.25791/self.zn2 \
                               + (-0.00635)*self.x1 + (-0.11654)*self.x1/self.zn1 \
                               + (-0.00193)*self.x2 + (-0.24188)*self.x2/self.zn2 \
                               + 0.00529*self.x1**2 + 0.00182*self.x2**2
        self.c_prime = self.C_M * self.C_R * self.C_B * math.cos(self.beta_rad) / self.q_prime
        self.c_gamma = self.c_prime * (0.75 * self.transverse_contact_ratio + 0.25)

        print("----- Mesh Stiffness -----")
        print(self.c_prime)
        print(self.c_gamma)

    def calculate_te(self, load_condition, n_steps, n_width, n_contact_lines=5, n_iteration=20, tol=0.001):

        # Print Setting
        np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

        # Input
        torque = load_condition.torque
        misalignment = load_condition.misalignment
        base_force = torque * 1000.0 / self.base_diameter_pinion * 2.

        print("----- Base Force -----")
        print(base_force)

        # Output
        te = []
        contact_length = []

        print("----- Steps -----")

        for i in range(n_steps):

            # Define Contact Line Related

            # Position in Active Plane

            x_pos = np.linspace(0.0, self.b_eff, n_width, endpoint=False, dtype=float)
            x_pos = x_pos + 0.5 * self.b_eff / n_width
            y_pos = np.zeros((n_contact_lines, n_width))

            for j in range(n_contact_lines):
                if self.active == "right":
                    y_pos[j] = x_pos * math.tan(self.beta_base_rad) - i * self.pbt / n_steps + (j - 2) * self.pbt
                else:
                    y_pos[j] = x_pos * math.tan(self.beta_base_rad) + i * self.pbt / n_steps + (j - 2) * self.pbt

            # print("----- Position in Active Plane Coord -----")
            # print(x_pos)
            # print(y_pos)

            # Position in Gear Coord
            # Outside the range of SAP and EAP will set to -1.

            if self.pinion_driving == 1:
                y_pinion = y_pos + self.sap_pinion
                y_wheel = self.pbt - y_pos + self.c - self.sap_pinion - self.pbt
            else:
                y_pinion = self.pbt - y_pos + self.c - self.sap_wheel - self.pbt
                y_wheel = y_pos + self.sap_wheel

            y_pinion = np.where((y_pinion - self.sap_pinion) >= -1e-3, y_pinion, -1.)
            y_pinion = np.where((y_pinion - self.eap_pinion) <= 1e-3, y_pinion, -1.)
            y_wheel = np.where((y_wheel - self.sap_wheel) >= -1e-3, y_wheel, -1.)
            y_wheel = np.where((y_wheel - self.eap_wheel) <= 1e-3, y_wheel, -1.)

            # print("----- Position in Gear Coord -----")
            # print(y_pinion)
            # print(y_wheel)

            # Relief in Active Plane
            # Relief is positive
            # Outside of the Range will set to negative

            relief_pinion = np.zeros((n_contact_lines, n_width))
            relief_wheel = np.zeros((n_contact_lines, n_width))
            relief_contact_line = np.zeros((n_contact_lines, n_width))

            for j in range(n_contact_lines):
                for k in range(n_width):
                    if self.b_eff < self.b_pinion:
                        relief_pinion[j][k] = self.cal_relief_pinion(x_pos[k]+0.5*(self.b_pinion-self.b_eff), y_pinion[j][k])
                    else:
                        relief_pinion[j][k] = self.cal_relief_pinion(x_pos[k], y_pinion[j][k])
                    if self.b_eff < self.b_wheel:
                        relief_wheel[j][k] = self.cal_relief_wheel(x_pos[k]+0.5*(self.b_wheel-self.b_eff), y_wheel[j][k])
                    else:
                        relief_wheel[j][k] = self.cal_relief_wheel(x_pos[k], y_wheel[j][k])
                    relief_contact_line[j][k] = relief_pinion[j][k] + relief_wheel[j][k]

            # print("----- Relief in Active Plane -----")
            # print(relief_contact_line)

            # Misalignment in Active Plane
            # Misalignment is positive
            # Outside the range will set to -1.

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

            # print("----- Misalignment in Active Plane -----")
            # print(misalignment_contact_line)

            # Sum of the Relief and Misalignment and drop the outside part

            relief_and_misalignment_contact_line = relief_contact_line + misalignment_contact_line
            relief_and_misalignment_contact_line_1 = relief_and_misalignment_contact_line[relief_and_misalignment_contact_line > -1e-6]

            # print("----- Sum of Relief and Misalignment -----")
            # print(relief_and_misalignment_contact_line_1)

            # Initial Value of Iteration
            temp_te = base_force / self.b_eff / self.c_gamma
            converge_or_not = 0
            ite_number = 0

            # Calculate TE in One Step

            for ite in range(n_iteration):

                deflection_contact_line_1 = temp_te - relief_and_misalignment_contact_line_1
                deflection_contact_line_1 = deflection_contact_line_1[deflection_contact_line_1 > 0.0]
                temp_deflection = np.sum(deflection_contact_line_1) * math.cos(self.beta_base_rad)
                diff = base_force - temp_deflection * self.b_eff / n_width * self.c_gamma

                # print("Step  ", ite)
                # print("te  ", temp_te)
                # print("deflection  ", temp_deflection)
                # print("diff  ", diff)
                # print(deflection_contact_line_1)
                # print("\n")

                if abs(diff) < base_force * tol:
                    te.append(temp_te)
                    temp_contact_length = self.b_eff / n_width * len(deflection_contact_line_1) / math.cos(self.beta_base_rad)
                    contact_length.append(temp_contact_length)
                    converge_or_not = 1
                    ite_number = ite
                    break
                else:
                    if diff < 0.0:
                        delta = min(abs(temp_te*0.5), abs(diff/self.b_eff/self.c_gamma))
                        delta = -1.*delta
                    else:
                        delta = min(abs(temp_te*2.), abs(diff/self.b_eff/self.c_gamma))
                    temp_te = temp_te + delta

            if converge_or_not == 1:
                print("step ", i, " converged at ", ite_number)
            else:
                print("step ", i, "not converged")

        te = np.array(te)
        te_peak = np.max(te) - np.min(te)
        contact_length = np.array(contact_length)

        print("----- TE -----")
        print(te)
        print(contact_length)
        print("Number of Steps ", len(te))
        print("Peak-Peak TE ", te_peak)
        print("Contact Length Max ", np.max(contact_length))
        print("Contact Length Min ", np.min(contact_length))

        return te, contact_length


# Define Load Condition - Pinion Torque and Misalignment
load_condition_1 = Load(500.0, 0.0)
load_condition_2 = Load(500.0, 10.0)
# Define Relief
relief_1 = Relief(0.0, 0.0, 0.0, 0.0, 0.0)
relief_2 = Relief(10.0, 10.0, 10.0, 10.0, 0.0)
# Define Gears
pinion_1 = Gear(31, 2.271, 20.0, 25.0, 20.0, 19.775, 8.06, relief_2)
wheel_1 = Gear(41, 2.271, 20.0, 25.0, 20.0, 24.142, 12.175, relief_2)
pinion_1.geometry()
wheel_1.geometry()
# Define Gear Pair
gear_pair_1 = GearPair(pinion_1, wheel_1, 90.0, 0.0)
gear_pair_1.mesh_geometry()
gear_pair_1.mesh_stiffness()
y1, y2 = gear_pair_1.calculate_te(load_condition_2, 32, 1000)
y1 = np.concatenate((y1,y1))
y2 = np.concatenate((y2,y2))
x = np.arange(64)
# Plot
plt.subplot(2, 1, 1)
plt.plot(x, y1)
plt.title("Basic LTCA")
plt.ylabel("TE")
plt.subplot(2, 1, 2)
plt.plot(x,y2)
plt.xlabel("Step")
plt.ylabel("Length of Contact")
plt.show()
