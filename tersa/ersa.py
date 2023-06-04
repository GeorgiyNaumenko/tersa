from eisenstein import EisensteinIntegers


class ERSA:

    def __init__(self, public_key: EisensteinIntegers, public_mod: EisensteinIntegers):
        self.gamma = public_mod
        self.eta = public_key

    def encrypt(self, message: EisensteinIntegers):
        return EisensteinIntegers.pow(message, self.eta.norm()) % self.gamma.norm()

    def decrypt(self, message: EisensteinIntegers, private_key: EisensteinIntegers):
        return EisensteinIntegers.pow(message, private_key.norm()) % self.gamma.norm()
