
class Load:

    def __init__(self, torque, misalignment):

        self.torque = torque
        self.misalignment = misalignment


class Relief:

    def __init__(self):

        pass


class Gear:

    def __init__(self, z, mn, alphan, beta, b):

        self.z = z
        self.mn = mn
        self.alphan = alphan
        self.beta = beta
        self.b = b


class GearPair:

    def __init__(self, pinion, wheel, a, x_pinion):

        self.a = a
        self.x_pinion = x_pinion

    def single_stiffness(self):

        pass

    def calculate_te(self, load, relief, n_contact_lines, n_steps, tol):

        pass

