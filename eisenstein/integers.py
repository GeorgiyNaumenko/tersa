import math


class EisensteinIntegers:
    """
    Eisenstein Integer initialization
    :param a: real part
    :param b: imaginary part
    """

    def __init__(self, a=0, b=0):
        self.a = a
        self.b = b
        self.units = None
        self.zero = None

    """
    Norm function
    :return: norm of current integer
    """

    def norm(self):
        return self.a * self.a - self.a * self.b + self.b * self.b

    """
    String representation of Eisenstein Integer
    :return: integer as string
    """
    def __str__(self):
        if self.b > 0:
            if self.b == 1:
                if self.a != 0:
                    return str(self.a) + ' + ' + 'w'
                else:
                    return 'w'
            else:
                if self.a != 0:
                    return str(self.a) + ' + ' + str(self.b) + 'w'
                else:
                    return str(self.b) + 'w'
        elif self.b < 0:
            if self.b == 1:
                if self.a != 0:
                    return str(self.a) + ' - w'
                else:
                    return '- w'
            else:
                if self.a != 0:
                    return str(self.a) + ' - ' + str(-self.b) + 'w'
                else:
                    return str(self.b) + 'w'
        else:
            return str(self.a)

    """
    Addition over Eisenstein Integers
    :param other: other integer
    :return: addition result
    """

    def __add__(self, other):
        return EisensteinIntegers(self.a + other.a, self.b + other.b)

    """
    Substraction over Eisenstein Integers
    :param other: other integer
    :return: substraction result
    """

    def __sub__(self, other):
        return EisensteinIntegers(self.a - other.a, self.b - other.b)

    def __neg__(self):
        return EisensteinIntegers(-self.a, -self.b)

    """
    Multipication over Eisenstein Integers
    :param other: other integer
    :return: multiplication result
    """

    def __mul__(self, other):
        return EisensteinIntegers(self.a * other.a - self.b * other.b,
                                  self.b * other.a + self.a * other.b - self.b * other.b)

    """
    Divisoin (for modular division over Eisenstein Integers)
    :param other: other integer
    :return: true division result
    """

    def __truediv__(self, other):
        return (self.a * other.a + self.b * other.b - self.a * other.b) / other.norm(), (
                    self.b * other.a - self.a * other.b) / other.norm()

    """
    Convertation Eisenstein to complex 
    :param a: real Eisenstein part
    :param b: imaginary Eisenstein part
    :return: complex number
    """

    def convertToComplex(self, a, b):
        return a - b / 2, b * math.sqrt(3) / 2

    """
    Convertation complex to Eisenstein
    :param a: real part
    :param b: imaginary part
    :return: Eisenstein
    """

    def convertFromComplex(self, a, b):
        return a + (b / math.sqrt(3)), 2 * b / math.sqrt(3)

    """
    Modular division over Eisenstein Integers
    :param other: other integer
    :return: quotient and remainder
    """

    def __divmod__(self, other):
        alpha = self
        beta = other

        temp1, temp2 = alpha / beta
        a1, b1 = self.convertToComplex(temp1, temp2)
        a2, b2 = a1 + 1 / 2, b1 - math.sqrt(3) / 2

        ro1_ = complex(a1 - self.nearestInteger(a1), b1 - math.sqrt(3) * self.nearestInteger(b1 / math.sqrt(3)))
        ro2_ = complex(a2 - self.nearestInteger(a2), b2 - math.sqrt(3) * self.nearestInteger(b2 / math.sqrt(3)))
        temp1, temp2 = self.convertToComplex(beta.a, beta.b)
        ro1 = complex(temp1, temp2) * ro1_
        ro2 = complex(temp1, temp2) * ro2_
        ro1a, ro1b = self.convertFromComplex(ro1.real, ro1.imag)
        ro2a, ro2b = self.convertFromComplex(ro2.real, ro2.imag)
        ro1 = EisensteinIntegers(self.nearestInteger(ro1a), self.nearestInteger(ro1b))
        ro2 = EisensteinIntegers(self.nearestInteger(ro2a), self.nearestInteger(ro2b))

        gamma1a, gamma1b = self.convertFromComplex(self.nearestInteger(a1),
                                                   math.sqrt(3) * self.nearestInteger(b1 / math.sqrt(3)))
        gamma2a, gamma2b = self.convertFromComplex(self.nearestInteger(a2),
                                                   math.sqrt(3) * self.nearestInteger(b2 / math.sqrt(3)))
        gamma1 = EisensteinIntegers(self.nearestInteger(gamma1a), self.nearestInteger(gamma1b))
        gamma2 = EisensteinIntegers(self.nearestInteger(gamma2a), self.nearestInteger(gamma2b)) + EisensteinIntegers(0, 1)

        if ro1.norm() < ro2.norm():
            return gamma1, ro1
        elif ro1.norm() > ro2.norm():
            return gamma2, ro2
        else:
            re1, im1 = self.convertToComplex(gamma1.a, gamma1.b)
            re2, im2 = self.convertToComplex(gamma2.a, gamma2.b)
            if re1 > re2:
                return gamma1, ro1
            else:
                return gamma2, ro2

    """
    Modular division over Eisenstein Integers
    :param other: other integer
    :return: quotient
    """

    def __floordiv__(self, other):
        q, r = divmod(self, other)
        return q

    """
    Modular division over Eisenstein Integers
    :param other: other integer
    :return: remainder
    """

    def __mod__(self, other):
        q, r = divmod(self, other)
        return r

    """
    Check primarity
    :return: True if primary else False
    """

    def isPrimary(self):
        return self.a % 3 == 2 and self.b % 3 == 0

    """
    Chek if unit
    :return: True if unit else False
    """

    def isUnit(self):
        return self.norm() == 1

    """
    Chek if zero
    :return: True if zero else False
    """

    def isZero(self):
        return self.a == 0 and self.b == 0

    """
    Conjugate calculation
    :return: conjugate to current
    """

    def conjugate(self):
        return EisensteinIntegers(self.a - self.b, -self.b)

    """
    Check if current is integer
    :return: True if integer else False
    """

    def isInteger(self):
        return self.b == 0

    """
    Finding reciprocal if it exists
    :return: current^(-1) if exists else None
    """

    def recip(self):
        if self == EisensteinIntegers(0, 1):
            return EisensteinIntegers(-1, -1)
        elif self == EisensteinIntegers(0, -1):
            return EisensteinIntegers(1, 1)
        elif self == EisensteinIntegers(-1, -1):
            return EisensteinIntegers(0, 1)
        elif self == EisensteinIntegers(1, 1):
            return EisensteinIntegers(0, -1)
        elif self == EisensteinIntegers(1, 0):
            return EisensteinIntegers(1, 0)
        elif self == EisensteinIntegers(-1, 0):
            return EisensteinIntegers(-1, 0)
        else:
            return None

    """
    Convertation to integer
    :return: integer if current is integer else False
    """

    def toInteger(self):
        if self.isInteger():
            return self.a
        return None

    """
    Searching for the nearest integer
    :param a: some complex
    :return: nearest integer
    """

    def nearestInteger(self, a):
        if abs(a - math.floor(a)) < 0.5:
            return math.floor(a)
        return math.ceil(a)

    """
    Extended Euclidean Algorithm
    :param other: other Eisenstein
    :return: gcd and pair (x, y) such self * x + other * y = gcd
    """

    def gcdExtended(self, other):
        a = self
        b = other
        if a.isZero():
            return (b, EisensteinIntegers(0, 0), EisensteinIntegers(1, 0))
        else:
            g, y, x = (b % a).gcdExtended(a)
            return g, x - (b // a) * y, y

    """
    Searching for inverse modulo m
    :param m: modulo
    :return: inverse modulo m if it exists else raise Exception
    """

    def modinv(self, m):
        units = [
            EisensteinIntegers(1, 0),
            EisensteinIntegers(-1, 0),
            EisensteinIntegers(0, 1),
            EisensteinIntegers(0, -1),
            EisensteinIntegers(-1, -1),
            EisensteinIntegers(1, 1),
        ]
        g, x, y = self.gcdExtended(m)
        if g not in units:
            raise Exception('modular inverse does not exist')
        else:
            return (x * g.recip()) % m


    """
    Equality check
    :param other: other integer
    :return: True if self == other else False
    """

    def __eq__(self, other):
        if self.a == other.a and self.b == other.b:
            return True
        return False

    @staticmethod
    def pow(x, y):
        number = EisensteinIntegers(1, 0)
        while y:
            if y & 1:
                number = number * x
            y >>= 1
            x = x * x
        return number

    def divide(self, other):
        return EisensteinIntegers(
            (self.a * other.a - self.a * other.b + self.b * other.b) / other.norm(),
            (self.b * other.a - self.a * other.b) / other.norm()
        )

    def __hash__(self):
        return hash(str(self))
