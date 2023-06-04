import math
import random

from .integers import EisensteinIntegers


class PrimalityTestEisenstein:
    """
  Primality test initialization
  """

    def __init__(self):
        self.units = [
            EisensteinIntegers(1, 0),
            EisensteinIntegers(-1, 0),
            EisensteinIntegers(0, 1),
            EisensteinIntegers(0, -1),
            EisensteinIntegers(-1, -1),
            EisensteinIntegers(1, 1),
        ]
        self.primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                       101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157]
        self.pseudocubes = [
            [5, 643, EisensteinIntegers(29, 18)],
            [7, 5113, EisensteinIntegers(71, 72)],
            [11, 13507, EisensteinIntegers(23, 126)],
            [13, 39199, EisensteinIntegers(227, 90)],
            [17, 107803, EisensteinIntegers(-181, 198)],
            [19, 360007, EisensteinIntegers(653, 126)],
            [23, 3904969, EisensteinIntegers(443, 2160)],
            [29, 6107191, EisensteinIntegers(-1669, 1170)],
            [31, 10318249, EisensteinIntegers(3617, 2520)],
            [37, 27333067, EisensteinIntegers(6023, 3366)],
            [41, 99179467, EisensteinIntegers(4973, 11466)],
            [43, 532997833, EisensteinIntegers(-15451, 11088)],
            [47, 2278522747, EisensteinIntegers(54017, 17514)],
            [53, 2741702809, EisensteinIntegers(47477, 56160)],
            [59, 18500766499, EisensteinIntegers(66887, 156510)],
            [61, 41547553813, EisensteinIntegers(235061, 107172)],
            [67, 41547553813, EisensteinIntegers(235061, 107172)],
            [71, 119423348797, EisensteinIntegers(-139813, 253764)],
            [73, 824621013649, EisensteinIntegers(-267733, 744120)],
            [79, 1151810360731, EisensteinIntegers(1227419, 761670)],
            [83, 1151810360731, EisensteinIntegers(1227419, 761670)],
            [89, 25079082769801, EisensteinIntegers(5052689, 4961880)],
            [97, 33932637528481, EisensteinIntegers(-2127709, 4462200)],
            [101, 91756768829893, EisensteinIntegers(10322861, 8601732)],
            [103, 214089061932079, EisensteinIntegers(3056387, 15918570)],
            [107, 214089061932079, EisensteinIntegers(3056387, 15918570)],
            [109, 812216615153761, EisensteinIntegers(-27791551, 1366560)],
            [113, 10706700434813749, EisensteinIntegers(109364777, 13014540)],
            [127, 15846955654747279, EisensteinIntegers(-114717193, 19952010)],
            [131, 21448509758341459, EisensteinIntegers(160585853, 126202050)],
            [137, 596036690644131739, EisensteinIntegers(845355437, 667764090)],
            [139, 2127627080411019739, EisensteinIntegers(-724036477, 954969030)],
            [149, 5736341949347177659, EisensteinIntegers(696254903, 2666049750)],
            [151, 9708823441723568077, EisensteinIntegers(2979509543, 3236384556)],
            [157, 14102281783170625921, EisensteinIntegers(3671532959, 3833807040)]
        ]

    """
  Fast power calculation (z ^ pow)
  :param z: some Eisenstein
  :param pow: power 
  :return: result of powering
  """

    def power(self, z, pow):
        if pow == 0:
            return EisensteinIntegers(1, 0)
        elif pow == 1:
            return z
        elif pow % 2 == 0:
            return self.power(z * z, pow // 2)
        elif pow % 2 != 0:
            return (z * self.power(z * z, (pow - 1) // 2))

    """
  Fast power calculation modulo mod (z ^ pow) % mod
  :param z: some Eisenstein
  :param pow: power
  :param mod: modulo
  :return: result of powering
  """

    def powerMod(self, z, pow, mod):
        res = EisensteinIntegers(1, 0)
        z = z % mod
        if z.isZero():
            return EisensteinIntegers(0, 0)
        while pow > 0:
            if (pow & 1) == 1:
                res = (res * z) % mod
            pow >>= 2
            z = (z * z) % mod
        return res

    """
  Check if alpha is cubic residue moulo pi
  :param alpha: some Eisenstein
  :param pi: prime Eisenstein
  :return: True if cubic residue else False
  """

    def isCubicResidue(self, alpha, pi):
        mod = pi
        return self.powerMod(alpha, (pi.norm() - 1) // 3, mod) == self.units[0]

    """
  Residue calculation
  :param alpha: some Eisenstein
  :param pi: prime Eisenstein
  :return: residue
  """

    def calculateResidue(self, alpha, pi):
        return self.powerMod(alpha, (pi.norm() - 1) // 3, pi)

    """
  Cubic Residue Character calculation
  :param alpha: some Eisenstein
  :param mod: some prime Eisenstein
  :return: CRC
  """

    def cubicResidueCharacter(self, alpha, mod):
        if alpha % mod == EisensteinIntegers(0, 0):
            return EisensteinIntegers(0, 0)
        if self.isCubicResidue(alpha, mod):
            return self.units[0]
        else:
            return self.units[2], self.units[4]

    """
  Test for power
  :param N: integer
  :return: True if N is perfect power else False
  """

    def perfectPowerTest(self, N):
        b_max = int(math.log2(N))
        a_max = int(math.sqrt(N))
        for b in range(2, b_max + 1):
            left = 0
            right = N - 1
            while left < right:
                a = left + (right - left) // 2
                if N > pow(a, b):
                    left = a + 1
                else:
                    right = a
            if pow(a, b) == N:
                return True
        return False

    """
  Extended Euclidean Algorithm
  :param a: integer
  :param b: integer
  :return: gcd and pair (x, y) such a * x + b * y = gcd
  """

    def egcd(self, a, b):
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = self.egcd(b % a, a)
            return (g, x - (b // a) * y, y)

    """
  Searching for inverse modulo
  :param a: integer
  :param m: modulo
  :return: a^(-1) if it exists else raise Exception
  """

    def modinv(self, a, m):
        g, x, y = self.egcd(a, m)
        if g != 1:
            raise Exception('modular inverse does not exist')
        else:
            return x % m

    """
  Searching for quadratic non-residue
  :param p: modulo
  :param k: number of trials
  :return: quadratic non-residue modulo p
  """

    def QNR(self, p, k):
        while k > 0:
            k -= 1
            z = random.randint(1, p - 1)
            if pow(z, (p - 1) // 2, p) == p - 1:
                return z

    """
  TonelliShanks algorithm
  :param n: integer
  :param p: modulo
  :return: r such r^2 = n mod p if exists else None
  """

    def TonelliShanks(self, n, p):
        s = 0
        q = p - 1
        while q % 2 == 0:
            s += 1
            q //= 2
        z = self.QNR(p, 10000000)
        M = s
        c = pow(z, q, p)
        t = pow(n, q, p)
        R = pow(n, (q + 1) // 2, p)
        if t == 0:
            return 0
        elif t == 1:
            return [R, p - R]
        while t != 1:
            i = 0
            t1 = t
            while t1 != 1:
                i += 1
                t1 = (t1 * t1) % p
                if i == M:
                    return None
            b = pow(c, pow(2, M - i - 1, p - 1), p)
            M = i
            c = pow(b, 2, p)
            t = (t * b * b) % p
            R = (R * b) % p
        return [R, p - R]

    """
  Cornacchia algoritnm (f, g, m) - pairly coprime
  :param f: integer
  :param g: integer
  :param m: integer
  :return: pair (x, y) such fx^2 + gy^2 = m if exists else None
  """

    def cornacchia(self, f, g, m):
        f1 = self.modinv(f, m)
        if f1 is None:
            return None
        K = self.TonelliShanks((-g * f1) % m, m)
        k = K[0]
        M = m
        x = k
        y = 1
        B = 0
        limit = math.sqrt(m / f)
        while x > limit:
            q, r = divmod(M, x)
            M = x
            x = r
            t = y
            y = q * y + B
            B = t
        if y < math.sqrt(m / g):
            return x, y
        return None

    """
  Searching for Eisenstein primary by norm
  :param N: norm
  :return: primary x such N(z) = N if exists else False
  """

    def findEisensteinByNorm(self, N):
        sol = self.cornacchia(1, 3, N)
        if sol is None:
            return False
        s, t = sol[0], sol[1]
        alpha1 = EisensteinIntegers(2 * t, s + t)
        alpha2 = EisensteinIntegers(s + t, 2 * t)
        alpha3 = EisensteinIntegers(-s - t, -2 * t)
        alpha4 = EisensteinIntegers(-2 * t, -s - t)
        if alpha1.isPrimary():
            return alpha1
        if alpha2.isPrimary():
            return alpha2
        if alpha3.isPrimary():
            return alpha3
        if alpha4.isPrimary():
            return alpha4
        return False

    """
  :param:
  :return:
  """

    def findPseudocube(self, N):
        left = 0
        right = len(self.pseudocubes) - 1
        while left < right:
            m = left + (right - left) // 2
            if N > self.pseudocubes[m][1]:
                left = m + 1
            else:
                right = m
        if self.pseudocubes[right][1] == N:
            return self.pseudocubes[right + 1]
        return self.pseudocubes[right]

    """
  Primality test
  :param N: integer congruent 1 modulo 3
  :return: True if prime else False
  """

    def primalityTest(self, N):
        if self.perfectPowerTest(N):
            return False
        nu = self.findEisensteinByNorm(N)
        if not nu:
            return False
        pseudocube = self.findPseudocube(N)
        p = pseudocube[0]
        muNorm = pseudocube[1]
        mu = pseudocube[2]
        q = self.primes[0]
        ind = 0
        while q <= p:
            val1 = self.calculateResidue(EisensteinIntegers(q, 0), nu)
            val2 = self.powerMod(EisensteinIntegers(q, 0), (N - 1) // 3, nu)
            if val1 != val2:
                return False
            ind += 1
            q = self.primes[ind]
        return True
