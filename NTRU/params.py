import random
import math

from eisenstein import PrimalityTestEisenstein, EisensteinIntegers, EisensteinPoly


class ParamsETRU:

  """
  ETRU Params initialization
  :param N: degree
  :param p: prime or prime power
  :param q: prime or prime power
  :param d_f: for f generation
  :param d_g: for g generation
  :param d_phi: for phi generation
  """
  def __init__(self, N, p, q, d_f=-1, d_g=-1, d_phi=-1):
    self.units = PrimalityTestEisenstein().units
    self.primalityTest = PrimalityTestEisenstein().primalityTest
    self.paramsDict = {'N': None, 'p': p, 'q': q, 'f': None, 'g': None, 'phi': None}
    if not N:
      self.paramsDict['N'] = self.GenerateN()
    else:
      self.paramsDict['N'] = N
    self.d_f = d_f
    if d_f == -1:
      self.d_f = math.floor(self.paramsDict['N']/7)
    self.d_g = d_g
    if d_g == -1:
      self.d_g = math.floor(self.paramsDict['N']/7)
    self.d_phi = d_phi
    if d_phi == -1:
      self.d_phi = math.floor(self.paramsDict['N']/7)
    self.paramsDict['f'] = self.GenerateF()
    self.paramsDict['g'] = self.GenerateG()
    self.paramsDict['phi'] = self.GeneratePhi()

  """
  Getting params
  :return: params dictionary
  """
  def getParams(self):
    return self.paramsDict

  """
  Generating integer N
  :return: N
  """
  def GenerateN(self):
    while True:
      k = random.randint(4, 2**8)
      N = 3 * k + 1
      if self.primalityTest(N):
        return N

  """
  Generating polynom f
  :return: f
  """
  def GenerateF(self):
    coefs = []
    for i in range(self.d_f):
      for el in self.units:
        coefs.append(el)
    coefs.append(EisensteinIntegers(1, 0))
    for i in range(self.paramsDict['N'] - self.d_f*6 - 1):
      coefs.append(EisensteinIntegers(0, 0))
    random.shuffle(coefs)
    return EisensteinPoly(coefs)

  """
  Generating polynom g
  :return: g
  """
  def GenerateG(self):
    coefs = []
    for i in range(self.d_g):
      for el in self.units:
        coefs.append(el)
    for i in range(self.paramsDict['N'] - self.d_g*6):
      coefs.append(EisensteinIntegers(0, 0))
    random.shuffle(coefs)
    return EisensteinPoly(coefs)

  """
  Generating Phi
  :return: phi
  """
  def GeneratePhi(self):
    coefs = []
    for i in range(self.d_phi):
      for el in self.units:
        coefs.append(el)
    for i in range(self.paramsDict['N'] - self.d_phi*6):
      coefs.append(EisensteinIntegers(0, 0))
    random.shuffle(coefs)
    return EisensteinPoly(coefs)