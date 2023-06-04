import copy

import numpy as np

from .integers import EisensteinIntegers


class EisensteinPoly:
    """
  Polynom with Eisenstein coefs initialization
  :param coefficients: coefficients
  """

    def __init__(self, coefficients=[]):
        self.coefficients = coefficients
        self.deg = max(len(self.coefficients) - 1, 0)
        self.normalize()

    """
  Check if zero
  :return: True if zero else False
  """

    def isZero(self):
        self.normalize()
        if len(self.coefficients) == 0:
            return True
        if self.coefficients[0] == EisensteinIntegers(0, 0) and self.deg <= 1:
            return True
        return False

    """
  String representation of polynom
  :return: polynom as string
  """

    def __str__(self):
        deg = self.deg
        coefs = self.coefficients
        s = ""
        for i in range(deg + 1):
            str_repr = str(coefs[i])
            if not coefs[i].isZero():
                if i == 1:
                    if coefs[i] != EisensteinIntegers(1, 0):
                        if s != "":
                            if ' ' in str_repr:
                                s += " + (" + str_repr + ")x"
                            else:
                                s += " + " + str_repr + "x"
                        else:
                            if ' ' in str_repr:
                                s += "(" + str_repr + ")x"
                            else:
                                s += str_repr + "x"
                    else:
                        if s != "":
                            s += " + x"
                        else:
                            s += "x"
                elif i == 0:
                    if s != "":
                        if ' ' in str_repr:
                            s += " + (" + str_repr + ")"
                        else:
                            s += " + " + str_repr
                    else:
                        if ' ' in str_repr:
                            s += "(" + str_repr + ")"
                        else:
                            s += str_repr
                else:
                    if coefs[i] != EisensteinIntegers(1, 0):
                        if s != "":
                            if ' ' in str_repr:
                                s += " + (" + str_repr + ")x^" + str(i)
                            else:
                                s += " + " + str_repr + "x^" + str(i)
                        else:
                            if ' ' in str_repr:
                                s += "(" + str_repr + ")x^" + str(i)
                            else:
                                s += str_repr + "x^" + str(i)
                    else:
                        if s != "":
                            s += " + x^" + str(i)
                        else:
                            s += "x^" + str(i)
        return s

    """
  Addition of polynoms
  :param other: other polynom
  :return: addition result
  """

    def __add__(self, other):
        coefs = []
        if self.deg > other.deg:
            poly1 = self
            poly2 = other
        else:
            poly1 = other
            poly2 = self
        coefs1 = poly1.coefficients
        coefs2 = poly2.coefficients
        n = poly1.deg
        m = poly2.deg
        for i in range(n - m):
            coefs.insert(0, coefs1[n - i])
        for i in range(m + 1):
            coefs.insert(0, coefs2[m - i] + coefs1[m - i])
        return EisensteinPoly(coefs)

    """
  Substraction of polynoms
  :param other: other polynom
  :return: substraction result
  """

    def __sub__(self, other):
        coefs = []
        if self.deg > other.deg:
            mul = EisensteinIntegers(1, 0)
            poly1 = self
            poly2 = other
        else:
            mul = EisensteinIntegers(-1, 0)
            poly1 = other
            poly2 = self
        coefs1 = poly1.coefficients
        coefs2 = poly2.coefficients
        n = poly1.deg
        m = poly2.deg
        for i in range(n - m):
            coefs.insert(0, mul * coefs1[n - i])
        for i in range(m + 1):
            coefs.insert(0, mul * (coefs1[m - i] - coefs2[m - i]))
        poly = EisensteinPoly(coefs)
        poly.normalize()
        return poly

    """
  Multiplication of polynoms
  :param other: other polynom
  :return: multiplication result
  """

    def __mul__(self, other):
        self.normalize()
        other.normalize()
        poly1 = self
        poly2 = other
        coefs1 = np.array(poly1.coefficients)
        coefs2 = np.array(poly2.coefficients)
        coefs1 = coefs1.reshape((1, coefs1.shape[0]))
        coefs2 = coefs2.reshape((1, coefs2.shape[0]))
        m1 = coefs1.T @ coefs2
        m2 = np.rot90(m1)
        start = -coefs2.shape[1] + 1
        coefs = np.array([np.sum(m2.diagonal(offset=start + i)) for i in range(coefs1.shape[1] + coefs2.shape[1] - 1)])
        return EisensteinPoly(list(coefs))

    """
  Modular division of polynoms
  :param other: other polynom
  :param p: prime Eisenstein (coefficients' ideal generator)
  :return: quotient, remainder
  """

    def divmod(self, other, p):
        f = copy.deepcopy(self)
        g = copy.deepcopy(other)
        if f.deg < g.deg:
            return EisensteinPoly([EisensteinIntegers(0, 0)]), f
        coefsF = np.flip(f.coefficients)
        coefsG = np.flip(g.coefficients)
        q = [EisensteinIntegers(0, 0)] * (f.deg - g.deg + 1)
        q[-1] = coefsF[0] * (coefsG[0].modinv(p))
        d = [EisensteinIntegers(0, 0)] * (f.deg - g.deg + 1)
        d[-1] = q[-1]
        temp = g * EisensteinPoly(d)
        temp.normalizeMod(p)
        r = f - temp
        r.normalizeMod(p)
        while r.deg >= g.deg:
            f = r
            coefsF = np.flip(f.coefficients)
            coefsG = np.flip(g.coefficients)
            d = [EisensteinIntegers(0, 0)] * (f.deg - g.deg + 1)
            q[len(d) - 1] = coefsF[0] * (coefsG[0].modinv(p))
            d[-1] = q[len(d) - 1]
            r = f - g * EisensteinPoly(d)
            r.normalizeMod(p)
            if f.deg * f.deg + g.deg * g.deg == 0:
                break
        return EisensteinPoly(q), r

    """
  Modular division of polynoms
  :param other: other polynom
  :param p: prime Eisenstein (coefficients' ideal generator)
  :return: quotient
  """

    def div(self, other, p):
        f = copy.deepcopy(self)
        g = copy.deepcopy(other)
        q, r = f.divmod(g, p)
        return q

    """
  Modular division of polynoms
  :param other: other polynom
  :param p: prime Eisenstein (coefficients' ideal generator)
  :return: remainder
  """

    def mod(self, other, p):
        f = copy.deepcopy(self)
        g = copy.deepcopy(other)
        q, r = f.divmod(g, p)
        return r

    """
  Extended Euclidean Algorithm
  :param other: other polynom
  :param p: prime Eisenstein (coefficients' ideal generator)
  :return: gcd, x, y such x*self + y*other = gcd
  """

    def gcdExtended(self, other, p):
        a = copy.deepcopy(self)
        b = copy.deepcopy(other)
        lastremainder = a
        remainder = b
        x = EisensteinPoly([EisensteinIntegers(0, 0)])
        lastx = EisensteinPoly([EisensteinIntegers(1, 0)])
        y = EisensteinPoly([EisensteinIntegers(1, 0)])
        lasty = EisensteinPoly([EisensteinIntegers(0, 0)])
        while not remainder.isZero():
            lastremainder, (quotient, remainder) = remainder, lastremainder.divmod(remainder, p)
            x, lastx = lastx - quotient * x, x
            y, lasty = lasty - quotient * y, y
        return lastremainder, lastx, lasty

    """
  Searching inverse modulo p
  :param other: other polynom, generator of polynomial ring
  :param p: prime Eisenstein (coefficients' ideal generator)
  :return: self^(-1) if exists else raise Exception
  """

    def modinv(self, other, p):
        g, x, y = self.gcdExtended(other, p)
        if g.deg == 0:
            inv = g.coefficients[0].modinv(p)
            mul = EisensteinPoly([inv])
            g = g * mul
            x = x * mul
            y = y * mul
        if g.deg != 0 or g.coefficients[-1] % p != EisensteinIntegers(1, 0):
            print('gcd:', g)
            raise Exception('modular inverse does not exist')
        else:
            res = x.mod(other, p)
            res.normalizeMod(p)
            return res

    """
  Normalizing of polynom all coefficients modulo
  :param mod: modulo
  """

    def normalizeMod(self, mod):
        for i in range(self.deg + 1):
            self.coefficients[i] = self.coefficients[i] % mod
        self.normalize()

    """
  Normalizing of polynom (deleting all first zeros)
  """

    def normalize(self):
        while self.coefficients and self.coefficients[-1].isZero():
            self.coefficients.pop()
            self.deg -= 1
        if len(self.coefficients) == 0:
            self.coefficients = [EisensteinIntegers(0, 0)]
            self.deg = 0

    def __call__(self, x):
        res = EisensteinIntegers(0, 0)
        for i, coef in enumerate(self.coefficients):
            res += coef * EisensteinIntegers.pow(x, i)
        return res
