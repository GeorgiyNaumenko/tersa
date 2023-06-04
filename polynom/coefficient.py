from .exceptions import *


class Coefficient:

    def __init__(self, number, field):
        """
        Coefficient object

        :param field: field for coefficient
        :param number: element from field with method rdivmod()
        """
        self.__field = field
        self.__number = field.rdivmod(number)

    def getNumber(self):
        """
        Method for taking the element

        :return: element
        """
        return self.__number

    def getField(self):
        """
        Method for takong the field

        :return: field
        """
        return self.__field

    def isZero(self):
        """
        Zero element or not

        :return: True if element is zero else False
        """
        return self.getNumber() == 0

    def __add__(self, other):
        """
        Overriding "+" over field

        :param other: other coefficient
        :return: new coefficient - result of addition
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self.getNumber(), other.getNumber(), self.getField(), other.getField())
        return self.getField().add(self.getNumber(), other.getNumber())

    def __sub__(self, other):
        """
        Overriding "-" over field

        :param other: other coefficient
        :return: new coefficient - result of substraction
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self.getNumber(), other.getNumber(), self.getField(), other.getField())
        return self.getField().sub(self.getNumber(), other.getNumber())

    def __mul__(self, other):
        """
        Overriding "*" over field

        :param other: other coefficient
        :return: new coefficient - result of multiplication
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self.getNumber(), other.getNumber(), self.getField(), other.getField())
        return self.getField().mul(self.getNumber(), other.getNumber())

    def __truediv__(self, other):
        """
        Overriding "/" over field

        :param other: other coefficient
        :return: new coefficient - result of division
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self.getNumber(), other.getNumber(), self.getField(), other.getField())
        return self.getField().div(self.getNumber(), other.getNumber())

    def __str__(self):
        """
        Overriding str()

        :return: string representation of the coefficient
        """
        return str(self.getNumber())

    def __eq__(self, other):
        """
        Overriding "==" over field

        :param other: other coefficient
        :return: True if coefficients are equal else False
        """
        if self.getField() != other.getField():
            raise DifferentFieldsError(self.getNumber(), other.getNumber(), self.getField(), other.getField())
        return self.getNumber() == other.getNumber()