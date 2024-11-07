import numpy as np
from scipy.integrate import solve_ivp
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
