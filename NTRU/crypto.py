from eisenstein import EisensteinPoly, EisensteinIntegers


class ETRU:
    """
  ETRU initialization
  :param: ETRU params dictionary
  """

    def __init__(self, params):
        self.N = params['N']
        self.p = params['p']
        self.q = params['q']
        self.f = params['f']
        self.g = params['g']
        self.phi = params['phi']
        polyModCoefs = [EisensteinIntegers(0, 0)] * (self.N + 1)
        polyModCoefs[0] = EisensteinIntegers(-1, 0)
        polyModCoefs[-1] = EisensteinIntegers(1, 0)
        self.polyMod = EisensteinPoly(polyModCoefs)

    """
  Key Generarion
  :return: {PK, SK} - public and secret keys
  """

    def keygen(self):
        self.Fp = (self.f).modinv(self.polyMod, self.p)
        self.Fq = (self.f).modinv(self.polyMod, self.q)
        self.h = (self.Fq * self.g).mod(self.polyMod, self.q)
        print('Fp * f', (self.Fp * self.f).mod(self.polyMod, self.p))
        print('Fq * f', (self.Fq * self.f).mod(self.polyMod, self.q))
        self.h.normalizeMod(self.q)
        public = {'h': self.h, 'p': self.p, 'q': self.q, 'N': self.N}
        secret = {'f': self.f}
        return {'public': public, 'secret': secret}

    """
  Encryption
  :param polyMessage: polynom for encrypt
  :return: polynom - encrypted message
  """

    def encrypt(self, polyMessage):
        polyP = EisensteinPoly([self.p])
        e = (polyP * self.phi * self.h + polyMessage).mod(self.polyMod, self.q)
        return e

    """
  Decryption
  :param encryptedPoly: - polynom for decrypt
  :return: decrypted polynom
  """

    def decrypt(self, encryptedPoly):
        a = (self.f * encryptedPoly).mod(self.polyMod, self.q)
        decrypted = (self.Fp * a).mod(self.polyMod, self.p)
        decrypted.normalizeMod(self.p)
        return decrypted
