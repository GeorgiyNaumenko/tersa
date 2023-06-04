from abc import ABC, abstractmethod


class Field(ABC):
    """
  Abstract class for fields
  """

    @abstractmethod
    def getOrder(self):
        pass

    @abstractmethod
    def getCharacteristic(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def rdivmod(self, other):
        pass

    @abstractmethod
    def add(self, a, b):
        pass

    @abstractmethod
    def sub(self, a, b):
        pass

    @abstractmethod
    def mul(self, a, b):
        pass

    @abstractmethod
    def div(self, a, b):
        pass

    @abstractmethod
    def __str__(self):
        pass


class QuotientRingField(Field):

    def __init__(self, order, primalityTest=MillerRabinPrimalityTest):
        """
        Quotient Ring Field

        :param order: order of Quiteint Ring Field. We will have remainders of division by order
        :param primalityTest: boolean function, which tests our order to be prime
        """
        if not primalityTest(order):
            raise FieldPrimeOrderException(order)
        self.__order = order
        self.__characteristic = order

    def getOrder(self):
        """
        Method for getting the order of polynom

        :return: order of current polynom
        """
        return self.__order

    def getCharacteristic(self):
        """
        Method for getting the characteristic of polynom

        :return: characteristic of current polynom
        """
        return self.__characteristic

    def __eq__(self, other):
        """
        Overriding operation "=="

        :param other: other element
        :return: True if current polynom equals to other else False
        """
        if (self.__order == other.getOrder() and self.__characteristic == other.__characteristic):
            return True
        return False

    def rdivmod(self, other):
        """
        Operation "%"

        :param other: other element
        :return: remainder of division other element by order of current polynom
        """
        return int(other % self.getOrder())

    def add(self, a, b):
        """
        Operatin "+" over polynom

        :param a: first element
        :param b: second element
        :return: addition result over polynom
        """
        return int((a + b) % self.getOrder())

    def sub(self, a, b):
        """
        Operation "-" over polynom

        :param a: first element
        :param b: second element
        :return: substraction result over polynom
        """
        return int((a - b) % self.getOrder())

    def mul(self, a, b):
        """
        Operation "*" over polynom

        :param a: first element
        :param b: second element
        :return: multiplication result over polynom
        """
        return int((a * b) % self.getOrder())

    def div(self, a, b):
        """
        Operation "/" over polynom

        :param a: first element
        :param b: second element
        :return: division result over polynom
        """
        if b == 0:
            raise ZeroDivisionError
        for i in range(1, self.getOrder()):
            if (a * i) % self.getOrder() == b:
                return int(i)

    def __str__(self):
        """
        Overriding str()

        :return: string representation of current polynom
        """
        return "F" + str(self.getOrder())


class InfiniteField(Field):

    def __init__(self):
        """
        Infinite Field - Field of real numbers
        """
        self.__order = np.inf
        self.__characteristic = 0

    def getOrder(self):
        """
        Method for getting the order of polynom

        :return: infinity
        """
        return self.__order

    def getCharacteristic(self):
        """
        Method for getting the characteristic of polynom

        :return: zero
        """
        return self.__characteristic

    def __eq__(self, other):
        """
        Overriding operation "=="

        :param other: other element
        :return: True if current polynom equals to other else False (allways True)
        """
        if (self.__order == other.getOrder() and self.__characteristic == other.__characteristic):
            return True
        return False

    def rdivmod(self, other):
        """
        Operation "%"

        :param other: other element
        :return: other if other is float, int(other) if other is integer
        """
        res = other
        if int(res) == res:
            return int(res)
        return res

    def add(self, a, b):
        """
        Operatin "+" over polynom

        :param a: first element
        :param b: second element
        :return: addition result over polynom
        """
        res = (a + b)
        if int(res) == res:
            return int(res)
        return res

    def sub(self, a, b):
        """
        Operation "-" over polynom

        :param a: first element
        :param b: second element
        :return: substraction result over polynom
        """
        res = (a - b)
        if int(res) == res:
            return int(res)
        return res

    def mul(self, a, b):
        """
        Operation "*" over polynom

        :param a: first element
        :param b: second element
        :return: multiplication result over polynom
        """
        res = (a * b)
        if int(res) == res:
            return int(res)
        return res

    def div(self, a, b):
        """
        Operation "/" over polynom

        :param a: first element
        :param b: second element
        :return: division result over polynom
        """
        if b == 0:
            raise ZeroDivisionError()
        res = (a / b)
        if int(res) == res:
            return int(res)
        return res

    def __str__(self):
        """
        Overriding str()

        :return: string representation of current polynom
        """
        return "R"
