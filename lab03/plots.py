from cProfile import label
from turtle import color
import matplotlib.pyplot as plt

#wyswietlanie otrzymanych wartosci
def plot_results(t_values_1, x_values_1, v_values_1, dt_values_1, t_values_2, x_values_2, v_values_2, dt_values_2, title):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(6, 6))

    ax1.set_xlabel('t')
    ax1.set_ylabel('x')
    ax1.plot(t_values_1, x_values_1, label='TOL = 0.01', color='tab:blue')
    ax1.plot(t_values_2, x_values_2, label='TOL = 0.00001', color='tab:red')
    ax1.legend(loc='upper left')

    ax2.set_xlabel('t')
    ax2.set_ylabel('v')
    ax2.plot(t_values_1, v_values_1, label='TOL = 0.01', color='tab:blue')
    ax2.plot(t_values_2, v_values_2, label='TOL = 0.00001', color='tab:red')
    ax2.legend(loc='upper left')

    ax3.set_xlabel('t')
    ax3.set_ylabel('dt')
    ax3.plot(t_values_1, dt_values_1, label='TOL = 0.01', color='tab:blue')
    ax3.plot(t_values_2, dt_values_2, label='TOL = 0.00001', color='tab:red')
    ax3.legend(loc='upper left')

    ax4.set_xlabel('x')
    ax4.set_ylabel('v')
    ax4.plot(x_values_1, v_values_1, label='TOL = 0.01', color='tab:blue')
    ax4.plot(x_values_2, v_values_2, label='TOL = 0.00001', color='tab:red')
    ax4.legend(loc='upper left')

    fig.suptitle(title)
    plt.tight_layout()
    plt.show()