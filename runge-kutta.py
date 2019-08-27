import numpy as np
import sympy as sym
from sympy.integrals.quadrature import gauss_lobatto
import matplotlib.pyplot as plt



class RungeKutta:

    def __init__(self, lamb, start_point, order, delta_t, t_stop, t_start=0):
        self.lamb = lamb
        self.start = start_point
        self.order = order
        self.delta_t = delta_t
        self.t_start = 0
        self.t_stop = t_stop
        self.s = int((order + 2.0) / 2.0)
        print("Stage is")
        print(self.s)
        self.amount_steps = int((t_stop - t_start) / delta_t)

    def create_time_step_vector(self, t_start, delta_t, t_stop):
        '''
        This function creates an numpy vector of all time steps used for
        integration.
        '''
        if t_stop < t_start:
            raise ValueError("The stop time is smaller then start time!")

        t_vector = np.linspace(t_start, t_stop, num=self.amount_steps)
        return t_vector


    def get_next_result(self, A, b, u_n):
        k = self.create_k_vector(self.s, self.lamb, A, u_n, self.delta_t)
        u_next = u_n + self.delta_t * np.dot(k, b)
        return u_next

    def run(self):
        c_vec = self.create_c_vector(self.s)
        A, b = self.create_A_and_b(self.s, c_vec)
        t_vector = self.create_time_step_vector(self.t_start, self.delta_t, self.t_stop)
        u_n = self.start
        solution = np.zeros(self.amount_steps)
        for count, time in enumerate(t_vector):
            if count == 0:
                solution[count] = u_n
            else:
                u_n = self.get_next_result(A, b, u_n)
                solution[count] = u_n
        return t_vector, solution


    def create_b_vector(self, s):
        '''
        Creates the b vector, size is 1xS
        '''
        pass

    def create_c_vector(self, s):
        c, w = gauss_lobatto(s, s)
        c = np.array(c)
        c = c + 1
        c = c * 0.5
        print("C is")
        print(c)
        return c

    def get_lagrange_term(self, a, b, t):
        return (t - a)/(b - a)

    def get_lagrange_polynomial(self, c_vec, j):
        t = sym.Symbol('t')
        Polynomial = sym.Rational(1)
        for a in range(len(c_vec)):
            if a != j:
                Polynomial *= self.get_lagrange_term(c_vec[a], c_vec[j], t)
        return Polynomial, t

    def create_A_and_b(self, s, c_vec):
        b = np.ones(s)
        A = np.ones((s,s))
        for count in range(s):
            lagrange, t = self.get_lagrange_polynomial(c_vec, count)
            t = sym.Symbol('t')
            b[count] = sym.Integral(lagrange, (t, 0, 1)).evalf()
            for count2 in range(s):
                A[count2][count] = sym.Integral(lagrange, (t, 0, c_vec[count2])).evalf()

        print("A is")
        print(A)
        print("b is")
        print(b)
        return A, b

    def create_k_vector(self, s, lamb, A, u_n, t_delta):
        b = -1.0 * lamb * u_n * np.ones(s)
        A_solve = t_delta * lamb * A - np.identity(s)
        k = np.linalg.solve(A_solve, b)
        return k


rung = RungeKutta(-2.0 * np.pi, 3.0, 1, 0.01, 7.183)
x, y = rung.run()
rung = RungeKutta(-2.0 * np.pi, 3.0, 4, 0.01, 7.183)
x,y_2 = rung.run()
exact = 3.0 * np.exp(-2.0 * np.pi * x)
error = np.abs(y_2[-1] - exact[-1])/np.abs(exact[-1])
print("Error for order 4 is: %f" % error)
print(error)
plt.plot(x, exact, 'bx--', label="exact")
plt.plot(x,y, 'ro', label="order 1")
plt.plot(x,y_2, 'go', label="order 4")
plt.xlabel('time [s]')
plt.ylabel('value')
plt.legend()
plt.show()
