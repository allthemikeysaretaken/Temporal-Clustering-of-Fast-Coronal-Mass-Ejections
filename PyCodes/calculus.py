import numpy as np
from scipy.integrate import quad

class Calculus():

    """
    purpose:
        (1) --  compute derivatives via centered-differences generally and
                edge-differences at the boundaries
        (2) --  integrate f(x) along iterable lower and
                upper bounds
        *(3) --  compute Jacobian
    """

    @staticmethod
    def bounded_integral(f, lbound, ubound, args=(), points=None):
        """
        f               :   type <function>
        lbounds         :   type <array>
        ubounds         :   type <array>
        args            :   type <tuple>
        points          :   type <sequence> or None
        """
        return quad(f, lbound, ubound, args=args, points=points)[0]

    @staticmethod
    def derivative(x, y):
        """
        x               :   type <array>
        y               :   type <array>
        """
        res = np.zeros_like(y)
        res[1:-1] = (x[2:] - x[:-2]) / (y[2:] - y[:-2])
        res[0] = (y[1] - y[0]) / (x[1] - x[0])
        res[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
        return res

    def nth_derivative(self, x, y, n):
        """
        x               :   type <array>
        y               :   type <array>
        n               :   type <int>
        """
        if not isinstance(n, int):
            raise ValueError("n = type <int>")
        dev = self.derivative(x, y)
        for i in range(n-1):
            dev = self.derivative(x, dev)
        return dev

    @staticmethod
    def compute_jacobian(f, args=()):
        """
        **************
        """
        raise ValueError("not yet implemented")
