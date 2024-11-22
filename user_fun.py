import numpy as np
from scipy.integrate import solve_ivp, odeint
import csv


def ff1T(x):
    return (
            -np.cos(0.1 * x) * np.e ** (-((0.1 * x - 2 * np.pi) ** 2))
            + 0.002 * (0.1 * x) ** 2
    )


def df1(t, state, D_A, params):
    # Unpacking the state variables
    V_A, V_B, T_B = state  # state variables at the current time t

    # Unpacking the parameters
    P_A, V_0A, T_0A, P_B, V_0B, T_0B, F_in_B, T_in_B, D_B, a, b, g = params

    # Flow rate out of Tank A into Tank B
    dV_A_dt = -a * b * D_A * np.sqrt(2 * g * V_A / P_A) if V_A > 0 else 0

    # Flow rate change in Tank B
    dV_B_dt = (
        (-a * b * D_B * np.sqrt(2 * g * V_B / P_B) + F_in_B + (-dV_A_dt))
        if V_B > 0
        else (F_in_B + (-dV_A_dt))
    )

    # Effective temperature of water entering Tank B (mixed from Tank A and external inflow)
    T_in_B = (
        ((-dV_A_dt) * T_0A + F_in_B * 20) / ((-dV_A_dt) + F_in_B)
        if ((-dV_A_dt) + F_in_B) > 0
        else 20
    )

    # Rate of temperature change in Tank B
    dT_B_dt = (((-dV_A_dt) + F_in_B) * (T_in_B - T_B)) / V_B if V_B > 0 else 0

    return [dV_A_dt, dV_B_dt, dT_B_dt]


def ff1R(D_A, params):
    P_A, V_0A, T_0A, P_B, V_0B, T_0B, F_in_B, T_in_B, D_B, a, b, g = params
    Y0 = [V_0A, V_0B, T_0B]
    t_span = (0, 2000)  # Zakres czasu dla symulacji

    # Rozwiązywanie ODE dla danego D_A
    sol = solve_ivp(
        lambda t, Y: df1(t, Y, D_A, params),
        t_span,
        Y0,
        method="RK45",
        max_step=1,
    )

    # Sprawdzenie, czy symulacja się powiodła
    if not sol.success:
        print("Symulacja nie powiodła się:", sol.message)
        return None

    # Znalezienie maksymalnej temperatury w zbiorniku B
    T_B_over_time = sol.y[2]
    with open("huj.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["V_A", "V_B", "T_B"])  # Header row
        for v_a, v_b, t_b in zip(sol.y[0], sol.y[1], sol.y[2]):
            # Convert to string to avoid type issues
            writer.writerow([str(v_a), str(v_b), str(t_b)])

    max_temp = max(T_B_over_time)

    return abs(max_temp - 50)


def ff2T(x):
    x1, x2 = x
    return x1 ** 2 + x2 ** 2 - np.cos(2.5 * np.pi * x1) - np.cos(2.5 * np.pi * x2) + 2


def tescik(x):
    x1, x2 = x
    return 2.5 * (x1 ** 2 - x2) ** 2 + (1 - x1) ** 2


def df2(Y, k1, k2):
    mr = 1.0  # masa ramienia
    mc = 5.0  # masa ciężarka
    l = 1.0  # długość ramienia
    alfa_ref = 180
    omega_ref = 0.0
    b = 0.5  # współczynnik tarcia
    I = (mr * l ** 2) / 3 + mc * l ** 2  # moment bezwładności

    dY = np.zeros(2)
    dY[0] = Y[1]
    dY[1] = (k1 * (alfa_ref - Y[0]) + k2 * (omega_ref - Y[1]) - b * Y[1]) / I

    return dY


def ff2R(x):
    alfa_ref = np.pi
    omega_ref = 0.0
    Y0 = [0.0, 0.0]  # początkowe warunki
    T = 100
    k1, k2 = x  # współczynniki optymalizacji
    dt = 0.1
    time_points = np.arange(0, T, dt)

    # Definicja funkcji do odeint
    def system(Y, t):
        dY = np.zeros(2)
        dY[0] = Y[1]
        mr = 1.0  # masa ramienia
        mc = 5.0  # masa ciężarka
        l = 1.0  # długość ramienia
        b = 0.5  # współczynnik tarcia
        I = (mr * l ** 2) / 3 + mc * l ** 2  # moment bezwładności
        dY[1] = (k1 * (alfa_ref - Y[0]) + k2 * (omega_ref - Y[1]) - b * Y[1]) / I
        return dY

    # Rozwiązanie równań różniczkowych
    Y = odeint(system, Y0, time_points)

    # Obliczanie funkcji celu
    integral_sum = 0
    for i in range(len(time_points)):
        alpha_diff = alfa_ref - Y[i, 0]
        omega_diff = omega_ref - Y[i, 1]
        func_val = (10 * alpha_diff ** 2 + omega_diff ** 2 +
                    (k1 * alpha_diff + k2 * omega_diff) ** 2) * dt
        integral_sum += func_val

    return integral_sum


# Funkcja celu
def ff3t(x):
    x1, x2 = x[0], x[1]
    denominator = np.pi * np.sqrt((x1 / np.pi) ** 2 + (x2 / np.pi) ** 2)
    if denominator == 0:
        return float('inf')  # Unikamy dzielenia przez 0
    return np.sin(np.pi * np.sqrt((x1 / np.pi) ** 2 + (x2 / np.pi) ** 2)) / denominator


# Ograniczenia
def constraint_g1(x):
    return -x[0] + 1  # g1(x1) <= 0 -> x1 >= 1


def constraint_g2(x):
    return -x[1] + 1  # g2(x2) <= 0 -> x2 >= 1


def constraint_g3(x):
    a = 5
    return np.sqrt(x[0] ** 2 + x[1] ** 2)- a  # g3(x1, x2) <= 0
