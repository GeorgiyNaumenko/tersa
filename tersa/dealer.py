import math
import random

from eisenstein import EisensteinIntegers, PrimalityTestEisenstein, EisensteinPoly

MAX_VAL = 10000


class Dealer:

    def __init__(self, t: int = None, n: int = None):
        self.n = n
        self.t = t
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.eta = None
        self.mu = None
        self.p = None
        self.q = None
        self.N = None
        self.phi = None
        self.poly = EisensteinPoly([])
        self.partial_secrets = None

    def generate_params(self, **kwargs):
        if not self.n or not self.t:
            raise AttributeError("Incorrect (t, n) specification.")
        if len(kwargs) > 0:
            pass            # maintenance of pre-calculated parameters
        else:
            while not (self.alpha and self.beta):
                self.p = Dealer.generate_int_prime_1mod3(MAX_VAL)
                self.q = Dealer.generate_int_prime_1mod3(MAX_VAL)
                while self.p == self.q:
                    self.q = Dealer.generate_int_prime_1mod3(MAX_VAL)
                self.alpha = Dealer.find_eis_by_norm(self.p)
                self.beta = Dealer.find_eis_by_norm(self.q)
                print(self.alpha, self.beta)
            self.gamma = self.alpha * self.beta
            self.N = self.gamma.norm()
            self.phi = (self.p - 1) * (self.q - 1)
            while not (self.eta and self.mu):
                eta, mu = self.generate_exponent(self.phi)
                self.eta = self.find_eis_by_norm(eta)
                self.mu = self.find_eis_by_norm(mu)
            self.poly = self.generate_polynom()
            self.partial_secrets = self.generate_partial_keys()

    @staticmethod
    def isprime(a):
        if a < 2: return False
        for x in range(2, int(math.sqrt(a)) + 1):
            if a % x == 0:
                return False
        return True

    @staticmethod
    def generate_int_prime_1mod3(maxval):
        p = random.randint(1, maxval//3 - 1)
        while not Dealer.isprime(3 * p + 1):
            p = random.randint(1, maxval//3 - 1)
        return 3 * p + 1

    def generate_polynom(self):
        if not self.mu:
            raise AttributeError("Common secret must be generated.")
        coefs = [self.mu]
        for i in range(1, self.t):
            coefs.append(Dealer.generate_random_ei_bounded(self.N))
        return EisensteinPoly(coefs)

    @staticmethod
    def find_eis_by_norm(n: int):
        test = PrimalityTestEisenstein()
        return test.findEisensteinByNorm(n)

    @staticmethod
    def generate_random_ei_bounded(n):
        def f():
            bound = int(math.sqrt(n))
            b = random.randint(-bound, bound)
            a_1 = int(1 - math.sqrt(4 * n) / 2)
            a_2 = int(math.sqrt(6 * n) / 2)
            a = random.randint(a_1, a_2)
            return a, b
        a, b = f()
        ei = EisensteinIntegers(a, b)
        while ei.norm() >= n:
            a, b = f()
            ei = EisensteinIntegers(a, b)
        return ei

    def generate_common_keys(self):
        pass

    def generate_partial_keys(self):
        points = []
        used = []
        for i in range(self.n):
            ei = Dealer.generate_random_ei_bounded(self.N)
            while ei in used or ei == EisensteinIntegers(0, 0):
                ei = Dealer.generate_random_ei_bounded(self.N)
            points.append([ei, self.poly(ei)])
        return points

    def reproduce_common_private_key(self, points):
        if len(points) < self.t:
            return None
        for u in points:
            print(f"used: {u[0]}, {self.N - u[0].norm()}")
        return self.lagrange_free_term(points)

    @staticmethod
    def ei_to_nearest_int(ei):
        return EisensteinIntegers(ei.nearestInteger(ei.a), ei.nearestInteger(ei.b))

    def lagrange_free_term(self, points):
        mu = EisensteinIntegers(0, 0)
        for j in range(self.t):
            partial_mu = EisensteinIntegers(1, 0)
            for k in range(self.t):
                if k != j:
                    partial_mu *= points[k][0].divide((points[k][0] - points[j][0]))
            mu += points[j][1] * partial_mu
        return Dealer.ei_to_nearest_int(mu)

    @staticmethod
    def extendedGCD(a, b):
        x, x1, y, y1 = 1, 0, 0, 1
        while b:
            q = a // b
            a, b = b, a % b
            x, x1 = x1, x - x1 * q
            y, y1 = y1, y - y1 * q
        return x, y, a

    def generate_exponent(self, phiN):
        while True:
            # e = random.randint(1, phiN)
            e = Dealer.generate_int_prime_1mod3(self.p + self.q + 1)
            d, _, a = Dealer.extendedGCD(e, phiN)
            if a == 1:
                return e % phiN, d % phiN
