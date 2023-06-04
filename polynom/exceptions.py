class FieldOrderException(Exception):

    def __init__(self, order, message="Incorrect order. Order must be a power of a prime number."):
        """
        Exception for incorrect order value

        :param order: incorrect order
        :param message: message for output
        """
        self.message = message
        self.order = order
        super().__init__(self.message)

    def __str__(self):
        """
        Overriding str()

        :return: string for output
        """
        return f'{self.order} -> {self.message}'


class FieldPrimeOrderException(Exception):

    def __init__(self, order, message="Incorrect order. Order must be a prime number."):
        """
        Exception for incorrect order value

        :param order: incorrect order
        :param message: message for output
        """
        self.message = message
        self.order = order
        super().__init__(self.message)

    def __str__(self):
        """
        Overriding str()

        :return: string for output
        """
        return f'{self.order} -> {self.message}'


class ZeroCoefficientError(Exception):

    def __init__(self, message="The coefficient of leading term is zero."):
        """
        Exception for first coefficient is zero

        :param message: message for output
        """
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        """
        Overriding str()

        :return: string for output
        """
        return f'{self.message}'


class DifferentFieldsError(Exception):

    def __init__(self, a, b, f1, f2, message="Not available to operate with elements from different fields"):
        """
        Exception for operations with elements from different fields

        :param a: first element
        :param b: second element
        :param a: first element's field
        :param a: second element's field
        :param message: message for output
        """
        self.message = message
        self.a = a
        self.b = b
        super().__init__(self.message)

    def __str__(self):
        """
        Overriding str()

        :return: string for output
        """
        return f'{self.a} from {f1} and {self.b} from {f2} -> {self.message}'