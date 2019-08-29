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
        vec = self.coarse_func(self.t_start, self.t_stop, self.u_0, self.delta_t_coarse)
        last_coarse = self.coarse_func(self.t_start, self.t_stop, self.u_0, self.delta_t_coarse)
        for k in range(self.K):
            updated_coarse = np.zeros(self.p+1)
            tmp = np.zeros(self.p+1)
            tmp[0] = vec[0] # initial value
            for p, time in zip(range(self.p), self.time_intervals):
                start, stop = time
                fine_calc = self.fine_func(start, stop, vec[p], self.delta_t)
                # print("Fine vector is")
                # print(fine_calc[-1])
                delta_q = fine_calc[-1] - last_coarse[p+1]
                # print("Delta q is")
                # print(delta_q)
                new_q = self.coarse_func(start, stop, tmp[p], self.delta_t_coarse)
                # print("new coarse is")
                # print(new_q)
                updated_coarse[p+1] = new_q[-1]
                tmp[p+1] = (updated_coarse[p+1] - last_coarse[p+1]) + fine_calc[-1]

                # print("old coarse was")
                # print(last_coarse[p+1])
                # print("fine calculation gives")
                # print(fine_calc[-1])
                # print("delta q gives")
                # print(delta_q)
                # print("new coarse gives")
                # print(updated_coarse[p+1])
                # print("Therefore new value is")
                # print(updated_coarse[p+1] + delta_q)
            vec = tmp
            last_coarse = updated_coarse
        return vec

    def run_result_per_iter(self):
        vec = self.coarse_func(self.t_start, self.t_stop, self.u_0, self.delta_t_coarse)
        last_coarse = self.coarse_func(self.t_start, self.t_stop, self.u_0, self.delta_t_coarse)
        per_iter_value = np.zeros((self.K, self.p+1))
        for k in range(self.K):
            updated_coarse = np.zeros(self.p+1)
            tmp = np.zeros(self.p+1)
            tmp[0] = vec[0] # initial value
            for p, time in zip(range(self.p), self.time_intervals):
                start, stop = time
                fine_calc = self.fine_func(start, stop, vec[p], self.delta_t)
                # print("Fine vector is")
                # print(fine_calc[-1])
                delta_q = fine_calc[-1] - last_coarse[p+1]
                # print("Delta q is")
                # print(delta_q)
                new_q = self.coarse_func(start, stop, tmp[p], self.delta_t_coarse)
                # print("new coarse is")
                # print(new_q)
                updated_coarse[p+1] = new_q[-1]
                tmp[p+1] = (updated_coarse[p+1] - last_coarse[p+1]) + fine_calc[-1]
                per_iter_value[k][p+1] = (updated_coarse[p+1] - last_coarse[p+1]) + fine_calc[-1]

                # print("old coarse was")
                # print(last_coarse[p+1])
                # print("fine calculation gives")
                # print(fine_calc[-1])
                # print("delta q gives")
                # print(delta_q)
                # print("new coarse gives")
                # print(updated_coarse[p+1])
                # print("Therefore new value is")
                # print(updated_coarse[p+1] + delta_q)
            vec = tmp
            last_coarse = updated_coarse
        return vec, per_iter_value

    def run_with_new(self, delta_t, start_time, stop_time, K=None):
        self.delta_t = delta_t
        self.amount_steps = int((stop_time - start_time)/delta_t)
        self.t_start = start_time
        self.t_stop = stop_time
        if K is not None:
            self.K = K
        return self.run()



def coarse_func(start, stop, initial_value, delta_t):
    x, y = rung.run_with_new(initial_value, delta_t, start, stop)
    return y


def fine_func(start, stop, initial_value, delta_t):
    x, y = rung4.run_with_new(initial_value, delta_t, start, stop)
    return y


if __name__ == "__main__":
    p = 8
    K_1 = 2
    K_2 = 10
    delta_t = 0.01

    rung4 = RungeKutta(-2.0 * np.pi, 1, 4, 1, 1, t_start=0)
    rung = RungeKutta(-2.0 * np.pi, 1, 2, 1, 1, t_start=0)
    para = Parareal(coarse_func, fine_func, p, K_1, delta_t, 3.0, start_time=0.0, stop_time=7.183)
    y, i_value = para.run_result_per_iter()
    print(y)
    para2 = Parareal(coarse_func, fine_func, p, K_2, delta_t, 3.0, start_time=0.0, stop_time=7.183)
    y2, i_value2 = para2.run_result_per_iter()
    time = np.linspace(0.0, 7.183, p+1)
    exact = 3.0 * np.exp(-2.0 * np.pi * time)
    plt.plot(time, exact, 'bx--', label="exact")
    plt.plot(time, y, 'ro', label="K=2")
    plt.plot(time, y2, 'go', label="K=5")
    plt.xlabel('time [s]')
    plt.ylabel('value')
    plt.legend()
    plt.show()

    iter_k_1 = np.linspace(1, K_1, K_1)
    iter_k_2 = np.linspace(1, K_2, K_2)

    # plt.plot(iter_k_1, i_delta[range(K_1), [-1]])
    plt.plot(iter_k_2, i_value2[range(K_2), [-1]])
    plt.xlabel('iterations K')
    plt.ylabel('value')
    plt.show()

    start_time = 0.0
    stop_time = 7.0

    para = Parareal(coarse_func, fine_func, p, 5, delta_t, 3.0, start_time=start_time, stop_time=stop_time)
    para2 = Parareal(coarse_func, fine_func, p, 10, delta_t, 3.0, start_time=start_time, stop_time=stop_time)
    para3 = Parareal(coarse_func, fine_func, p, 20, delta_t, 3.0, start_time=start_time, stop_time=stop_time)

    compare_num = 10

    time = np.linspace(start_time, stop_time, p+1)
    exact = 3.0 * np.exp(-2.0 * np.pi * time)

    error_1 = np.zeros(compare_num)
    error_2 = np.zeros(compare_num)
    error_3 = np.zeros(compare_num)
    dt_axis = []
    dt_axis.append(0.25)
    for count in range(compare_num):
        delta_t = dt_axis[count]
        y = para.run_with_new(delta_t, start_time, stop_time)
        print(y)
        # error_1[count] = np.average(np.abs(y - exact) / np.abs(exact))
        error_1[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
        y = para2.run_with_new(delta_t, start_time, stop_time)
        # error_2[count] = np.average(np.abs(y - exact) / np.abs(exact))
        error_2[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
        y = para3.run_with_new(delta_t, start_time, stop_time)
        # error_3[count] = np.average(np.abs(y - exact) / np.abs(exact))
        error_3[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
        if count != compare_num-1:
            dt_axis.append(delta_t/2.0)

    plt.loglog(dt_axis, error_1, 'rx--', label="K=5")
    plt.loglog(dt_axis, error_2, 'gx--', label="K=10")
    plt.loglog(dt_axis, error_3, 'yx--', label="K=20")
    plt.legend()
    plt.show()


    K = 10

    for k in range(K):
        y = para.run_with_new(0.05, start_time, stop_time, K=k)
        error_1[count] = np.abs(y[-1] - exact[-1]) / np.abs(exact[-1])
    k_axis = np.linspace(1, K, K)
    plt.semilogy(k_axis, error_1, 'bx--')
    plt.show()
