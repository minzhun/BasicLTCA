import numpy as np


class Gear:

    def __init__(self, base_pitch, b, beta, Sap, Eap, stiffness, contact_ratio):

        # base pitch , mm
        self.base_pitch = base_pitch
        # face width , mm
        self.b = b
        # helix angle , degree
        self.beta_deg = beta
        # start of active profile , mm , rolling distance
        self.Sap = Sap
        # end of active profile , mm , rolling distance
        self.Eap = Eap
        # single stiffness , N/um/mm
        self.stiffness = stiffness
        # contact ratio , define the number of contact lines
        self.contact_ratio = contact_ratio

        # defined in lead_relief
        self.lead_relief_crowning = 0.0 # um
        self.lead_relief_crowning_start = 0.0 # mm
        self.lead_relief_crowning_end = 0.0 # mm

    def define_lead_relief_crowning(self, relief, start, end):
        self.lead_relief_crowning = relief
        self.lead_relief_crowning_start = start
        self.lead_relief_crowning_end = end

    def cal_lead_relief_crowing(self, x):
        temp1 = x - 0.5 * (self.lead_relief_crowning_end - self.lead_relief_crowning_start)
        temp2 = (self.lead_relief_crowning_end - self.lead_relief_crowning_start) * 0.5
        return temp1 * temp1 / temp2 / temp2 * self.lead_relief_crowning

    def profile_relief(self):
        pass


class BasicLTCA:

    # strip number
    n_strip = 30
    # step number in one mesh cycle
    n_step = 16
    n_cycle = 2

    def __init__(self, Gear, misalignment,):

        # TE result
        self.te = np.zeros(BasicLTCA.n_step * BasicLTCA.n_cycle)
        self.n_contact = Gear.contact_ratio





