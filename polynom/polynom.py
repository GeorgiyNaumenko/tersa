from math import ceil, log

import numpy as np

from coefficient import *


class Polynom:

    def __init__(self, field, *coefficients, allowedFirstZero=False):
        """
        Polynom structure with coefficients over field

        :param field: coefficients' field
        :param coefficient: array of coefficients
        :param allowedFirstZero: boolean, if True - first zero coefficients will be sliced
        """
        if len(coefficients) == 0:
            self.__deg = 0
            self.__coefficients = np.array([Coefficient(0, field)])
            self.__field = field
        else:
            self.__deg = len(coefficients) - 1
            self.__coefficients = np.flip(np.array([Coefficient(coefficient, field) for coefficient in coefficients]))
            if self.__coefficients[-1].isZero() and not allowedFirstZero:
                raise ZeroCoefficientError()
            if allowedFirstZero:
                while self.__coefficients[-1].isZero():
                    if self.__deg > 0:
                        self.__coefficients = np.delete(self.__coefficients, -1)
                        self.__deg -= 1
                    else:
                        self.__deg = 0
                        self.__coefficients = np.array([Coefficient(0, field)])
                        break
            self.__field = field

    def getDeg(self):
        """
        Method for taking the degree of polynom

        :return: degree
        """
        return self.__deg

    def getField(self):
        """
        Method for taking the coefficients' field

        :return: field
        """
        return self.__field

    def getCoefficients(self):
        """
        Method for taking the coefficients

        :return: coefficients array
        """
        return self.__coefficients

    def __getitem__(self, deg):
        """
        Overriding indexing

        :return: deg degree coefficient
        """
        return self.getCoefficients()[deg]

    def __str__(self):
        """
        Overriding str()

        :return: string representation of polynom
        """
        deg = self.getDeg()
        s = ""
        for i in range(deg + 1):
            if self.getCoefficients()[i].getNumber() != 0:
                if i == 1:
                    if self.getCoefficients()[i].getNumber() > 0:
                        if self.getCoefficients()[i].getNumber() != 1:
                            if s != "":
                                s += " + " + self.getCoefficients()[i].__str__() + "x"
                            else:
                                s += self.getCoefficients()[i].__str__() + "x"
                        else:
                            if s != "":
                                s += " + x"
                            else:
                                s += "x"
                    else:
                        if self.getCoefficients()[i].getNumber() != -1:
                            if s != "":
                                s += " - " + self.getCoefficients()[i].__str__()[1:] + "x"
                            else:
                                s += self.getCoefficients()[i].__str__() + "x"
                        else:
                            if s != "":
                                s += " - x"
                            else:
                                s += "x"

                elif i == 0:
                    if self.getCoefficients()[i].getNumber() > 0:
                        if s != "":
                            s += " + " + self.getCoefficients()[i].__str__()
                        else:
                            s += self.getCoefficients()[i].__str__()
                    else:
                        if s != "":
                            s += " - " + self.getCoefficients()[i].__str__()[1:]
                        else:
                            s += self.getCoefficients()[i].__str__()

                else:
                    if self.getCoefficients()[i].getNumber() > 0:
                        if self.getCoefficients()[i].getNumber() != 1:
                            if s != "":
                                s += " + " + self.getCoefficients()[i].__str__() + "x^" + str(i)
                            else:
                                s += self.getCoefficients()[i].__str__() + "x^" + str(i)
                        else:
                            if s != "":
                                s += " + x^" + str(i)
                            else:
                                s += "x^" + str(i)
                    else:
                        if self.getCoefficients()[i].getNumber() != -1:
                            if s != "":
                                s += " - " + self.getCoefficients()[i].__str__()[1:] + "x^" + str(i)
                            else:
                                s += self.getCoefficients()[i].__str__() + "x^" + str(i)
                        else:
                            if s != "":
                                s += " - x^" + str(i)
                            else:
                                s += "-x^" + str(i)
        return s

    def isZero(self):
        """
        Zero polynom or not

        :return: True if polynom is zero else False
        """
        return self.getDeg() == 0 and self.getCoefficients()[0] == Coefficient(0, self.getField())

    def valueInField(self, x):
        """
        Calculating value in given point over field

        :param x: point
        :return: polynom value in given point over field
        """
        x_array = x * np.ones(self.getDeg() + 1)
        coefs = np.array(list(map(lambda x: x.getNumber(), self.getCoefficients())))
        val = np.sum(np.power(x_array, np.arange(self.getDeg() + 1)) * coefs)
        return self.getField().rdivmod(val)

    def value(self, x):
        """
        Calculating value in given point

        :param x: point
        :return: polynom value over field
        """
        x_array = x * np.ones(self.getDeg() + 1)
        coefs = np.array(list(map(lambda x: x.getNumber(), self.getCoefficients())))
        val = np.sum(np.power(x_array, np.arange(self.getDeg() + 1)) * coefs)
        return val

    def derivative(self):
        """
        Calculating the direvative of polynom

        :retrun: polynom - derivative
        """
        field = self.getField()
        coefs_current = np.array(list(map(lambda x: x.getNumber(), self.getCoefficients())))
        coefs = np.flip(np.arange(self.getDeg() + 1) * coefs_current)
        poly = Polynom(field, *coefs, allowedFirstZero=True)
        return poly

    def __add__(self, other):
        """
        Overriding "+" for polynoms over field

        :param other: other polynom
        :return: new polynom - result of addition
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self, other)
        field = self.getField()

        poly1 = self
        poly2 = other
        if self.getDeg() < other.getDeg():
            poly2 = self
            poly1 = other

        m = poly2.getDeg()
        n = poly1.getDeg()
        r = n - m
        coefs1 = np.array(list(map(lambda x: x.getNumber(), poly1.getCoefficients())))
        coefs2 = np.array(list(map(lambda x: x.getNumber(), poly2.getCoefficients())))
        coefs2 = np.append(coefs2, np.zeros(r))
        coefs = np.flip(coefs1 + coefs2)
        return Polynom(field, *coefs, allowedFirstZero=True)

    def __sub__(self, other):
        """
        Overriding "-" for polynoms over field

        :param other: other polynom
        :return: new polynom - result of substraction
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self, other)
        field = self.getField()

        poly1 = self
        poly2 = other
        mul = 1
        if self.getDeg() < other.getDeg():
            poly2 = self
            poly1 = other
            mul = -1

        m = poly2.getDeg()
        n = poly1.getDeg()
        r = n - m
        coefs1 = np.array(list(map(lambda x: x.getNumber(), poly1.getCoefficients())))
        coefs2 = np.array(list(map(lambda x: x.getNumber(), poly2.getCoefficients())))
        coefs2 = np.append(coefs2, np.zeros(r))
        coefs = np.flip(mul * (coefs1 - coefs2))
        return Polynom(field, *coefs, allowedFirstZero=True)

    def __mul__(self, other):
        """
        Overriding "*" for polynoms over field

        :param other: other polynom
        :return: new polynom - result of multiplication
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self, other)
        poly1 = self
        poly2 = other
        field = self.getField()
        coefs1 = np.flip(np.array(list(map(lambda x: x.getNumber(), poly1.getCoefficients()))))
        coefs2 = np.flip(np.array(list(map(lambda x: x.getNumber(), poly2.getCoefficients()))))
        coefs1 = coefs1.reshape((1, coefs1.shape[0]))
        coefs2 = coefs2.reshape((1, coefs2.shape[0]))
        m1 = coefs1.T @ coefs2
        m2 = np.rot90(m1)
        start = -coefs2.shape[1] + 1
        coefs = np.array([np.sum(m2.diagonal(offset=start + i)) for i in range(coefs1.shape[1] + coefs2.shape[1] - 1)])
        return Polynom(field, *coefs, allowedFirstZero=True)

    def scale(self, n):
        """
        Computation q(x) = [p(x) / x^n] from http://web.cs.iastate.edu/~cs577/handouts/polydivide.pdf

        :return: new polynom - result of computation
        """
        coefs = np.flip(list(map(lambda x: x.getNumber(), self.getCoefficients())))
        if n < 0:
            return Polynom(self.getField(), *(coefs[:n]), allowedFirstZero=True)
        coefs = coefs.reshape((1, coefs.shape[0]))
        assert int(n) == n
        n = int(n)
        z = np.zeros(n).reshape((1, n))
        c = np.concatenate((coefs, z), axis=1)[0]
        return Polynom(self.getField(), *c, allowedFirstZero=True)

    def recip(self):
        """
        Realization of Reciprocal() from http://web.cs.iastate.edu/~cs577/handouts/polydivide.pdf
        It computes [x^(2k) / p(x)], deg(p(x)) = k

        :return: new polynom - result of Reciprocal()
        """
        field = self.getField()
        coefs = np.flip(np.array(list(map(lambda x: x.getNumber(), self.getCoefficients()))))
        k = self.getDeg() + 1
        if k == 1:
            return Polynom(field, Coefficient(1, field) / Coefficient(coefs[0], field), allowedFirstZero=True)
        poly = Polynom(field, *(coefs[:int(k / 2)]), allowedFirstZero=True)
        q = poly.recip()
        p1 = Polynom(field, 2, allowedFirstZero=True) * q
        r = p1.scale((3 / 2) * k - 2) - ((q * q) * self)
        return r.scale(-k + 2)

    def inverse(self, f):
        """
        Fast algorithm for finding the inverse polynom over field: field[x] / (f), (f) - ideal, generated by f
        https://www.lirmm.fr/arith18/papers/kobayashi-AlgorithmInversionUsingPolynomialMultiplyInstruction.pdf

        :return: inverse polynom
        """
        if self.getField() != f.getField():
            raise DifferentFieldsError(self, f)
        field = self.getField()
        s = f
        r = self
        v = Polynom(field, 0, allowedFirstZero=True)
        u = Polynom(field, 1, allowedFirstZero=True)
        while r.getDeg() != 0:
            d = s.getDeg() - r.getDeg()
            if d < 0:
                temp = s
                s = r
                r = temp
                temp = v
                v = u
                u = temp
                d = -d
            tempCoefs = list(np.zeros(d + 1))
            tempCoefs[0] = 1
            tempPoly = Polynom(field, *tempCoefs, allowedFirstZero=True)
            s1 = s - tempPoly * r
            v1 = v - tempPoly * u
            s = s1
            v = v1
        return u

    def __divmod__(self, other):
        """
        Overriding divmod() for polynoms over field

        Realization of fast polynomial division
        http://web.cs.iastate.edu/~cs577/handouts/polydivide.pdf

        :param other: other polynom
        :return: new polynom - result of addition
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self, other)
        poly1 = self
        poly2 = other
        field = self.getField()
        coefs1 = np.flip(list(map(lambda x: x.getNumber(), poly1.getCoefficients())))
        coefs2 = np.flip(list(map(lambda x: x.getNumber(), poly2.getCoefficients())))
        m = poly1.getDeg()
        n = poly2.getDeg()

        d = int(2 ** ceil(log(n + 1, 2))) - n - 1
        poly1_ = poly1.scale(d)
        poly2_ = poly2.scale(d)
        m_ = m + d
        n_ = n + d
        s = poly2_.recip()
        temp = poly1_ * s
        q = temp.scale(-2 * n_)
        if m_ > 2 * n_:
            unit = Polynom(field, 1)
            t1 = s * poly2_
            t = unit.scale(2 * n_) - t1
            t2 = poly1_ * t
            t3 = t2.scale(-2 * n_)
            q2, r2 = divmod(t3, poly2)
        r = self - (other * q)
        return q, r

    def __floordiv__(self, other):
        """
        Overriding "//" for polynoms over field

        :param other: other polynom
        :return: new polynom - result of floor division
        """
        quotient, remainder = divmod(self, other)
        return quotient

    def __mod__(self, other):
        """
        Overriding "%" for polynoms over field

        :param other: other polynom
        :return: new polynom - remainder of division
        """
        quotient, remainder = divmod(self, other)
        return remainder