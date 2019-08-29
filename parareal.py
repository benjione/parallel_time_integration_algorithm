import numpy as np
from RungeKutta import RungeKutta
import matplotlib.pyplot as plt


class Parareal:

    def __init__(self, coarse_function, fine_function, p, K, delta_t, initial_value, start_time=0.0, stop_time=1.0):
        self.coarse_func = coarse_function
        self.fine_func = fine_function
        self.delta_t = delta_t
        self.u_0 = initial_value
        self.K = K
        self.p = p
        self.amount_steps = int((stop_time - start_time)/delta_t)
        # allign time steps to number of p
        self.amount_steps = self.amount_steps - (self.amount_steps % p)
        self.coarse_vector = np.zeros(p)
        self.delta_t_coarse = (stop_time - start_time) / p       # coarse time steps
        self.time_intervals = self.split_time_in_intervals(start_time, stop_time, p)
        self.t_start = start_time
        self.t_stop = stop_time

    def split_time_in_intervals(self, t_start, t_stop, p):
        length_interval = self.delta_t_coarse
        interval_list = []
        for c in range(p):
            int_start = c*length_interval + t_start
            int_stop = (c+1)*length_interval + t_start
            interval_list.append((int_start, int_stop))
        return interval_list

    def run(self):
        vec = self.coarse_func(self.t_start, self.t_stop+self.delta_t, self.u_0, self.delta_t_coarse)
        for k in range(self.K):
            tmp = np.zeros(self.p+1)
            tmp[0] = vec[0] # initial value
            for p, time in zip(range(self.p), self.time_intervals):
                start, stop = time
                fine_calc = self.fine_func(start, stop, vec[p], self.delta_t)
                # print("Fine vector is")
                # print(fine_calc)
                delta_q = fine_calc[-1] - vec[p+1]
                # print("Delta q is")
                # print(delta_q)
                new_q = self.coarse_func(start, stop+self.delta_t, tmp[p], self.delta_t_coarse)
                # print("new coarse is")
                # print(new_q)
                new_q = new_q[-1]
                tmp[p+1] = new_q + delta_q
            vec = tmp
        return vec

    def run_with_new(self, delta_t):
        self.delta_t = delta_t
        self.amount_steps = int((stop_time - start_time)/delta_t)



def coarse_func(start, stop, initial_value, delta_t):
    x, y = rung.run_with_new(initial_value, delta_t, start, stop)
    return y


def fine_func(start, stop, initial_value, delta_t):
    x, y = rung4.run_with_new(initial_value, delta_t, start, stop)
    return y


if __name__ == "__main__":
    p = 30
    K = 10
    delta_t = 0.1

    rung4 = RungeKutta(-2.0 * np.pi, 1, 4, 1, 1, t_start=0)
    rung = RungeKutta(-2.0 * np.pi, 1, 2, 1, 1, t_start=0)
    para = Parareal(coarse_func, fine_func, p, 2, delta_t, 3.0, start_time=0.0, stop_time=7.183)
    y = para.run()
    para2 = Parareal(coarse_func, fine_func, p, 5, delta_t, 3.0, start_time=0.0, stop_time=7.183)
    y2 = para2.run()
    time = np.linspace(0.0, 7.183, p+1)
    exact = 3.0 * np.exp(-2.0 * np.pi * time)
    plt.plot(time, exact, 'bx--', label="exact")
    plt.plot(time, y, 'ro', label="K=2")
    plt.plot(time, y2, 'go', label="K=5")
    plt.xlabel('time [s]')
    plt.ylabel('value')
    plt.legend()
    plt.show()

    start_time = 0.0
    stop_time = 7.0

    compare_num = 10
    error_1 = np.zeros(compare_num)
    error_2 = np.zeros(compare_num)
    error_3 = np.zeros(compare_num)
    dt_axis = []
    dt_axis.append(0.25)
    for count in range(compare_num):
        delta_t = dt_axis[count]
        tmp, y = para.run_with_new(3.0, delta_t, start_time, stop_time)
        exact = 3.0 * np.exp(-2.0 * np.pi * tmp)
        # error_1[count] = np.average(np.abs(y - exact) / np.abs(exact))
        error_1[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
        tmp, y = rung_2.run_with_new(3.0, delta_t, start_time, stop_time)
        # error_2[count] = np.average(np.abs(y - exact) / np.abs(exact))
        error_2[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
        tmp, y = rung_3.run_with_new(3.0, delta_t, start_time, stop_time)
        # error_3[count] = np.average(np.abs(y - exact) / np.abs(exact))
        error_3[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
        if count != compare_num-1:
            dt_axis.append(delta_t/2.0)
    plt.loglog(dt_axis, error_1, 'rx--', label="order 2")
    plt.loglog(dt_axis, error_2, 'gx--', label="order 4")
    plt.loglog(dt_axis, error_3, 'yx--', label="order 6")
    plt.legend()
    plt.show()
