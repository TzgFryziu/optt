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
    return np.sqrt(x[0] ** 2 + x[1] ** 2) - a  # g3(x1, x2) <= 0


def ff4t(x):
    x1, x2 = x
    return (x1 + 2 * x2 - 7) ** 2 + (2 * x1 + x2 - 5) ** 2


def h0(theta, x):
    return 1 / (1 + np.exp((-theta).T @ x))


def compute_cost(theta, X, y):
    m = len(y)
    sum = 0
    for i in range(m):
        sum += y[i] * np.log(h0(theta, X[i])) + (1 - y[i]) * np.log(1 - h0(theta, X[i]))
    return -1 / m * sum


def compute_gradient(theta, X, y):
    m = len(y)
    gradient = np.zeros(theta.shape)
    for j in range(len(theta)):
        sum_error = 0
        for i in range(m):
            prediction = h0(theta, X[i])
            sum_error += (prediction - y[i]) * X[i][j]
        gradient[j] = (1 / m) * sum_error

    return gradient


def load_data(x_path, y_path):
    X = np.loadtxt(x_path, delimiter=';')
    y = np.loadtxt(y_path, delimiter=';')
    y = y.flatten()  # Upewnij się, że y jest wektorem płaskim
    x2 = []
    for i in range(len(X[0])):
        x2.append(np.array([X[0][i], X[1][i], X[2][i]]).reshape(3, 1))
    return x2, y


def przewidywanie(theta, X):
    result = []
    for i in range(len(X)):
        prediction = h0(theta, X[i])
        result.append(prediction)
    for i in range(len(result)):
        if result[i] >= 0.5:
            result[i] = 1
        else:
            result[i] = 0
    return result


def ocena_klasyfikatora(theta, X, y):
    y_pred = przewidywanie(theta, X)
    suma = 0
    for i in range(len(y)):
        if y_pred[i] == y[i]:
            suma += 1

    return suma / len(y)


def ff6t(x):
    x1, x2 = x
    return x1 ** 2 + x2 ** 2 - np.cos(2.5 * np.pi * x1) - np.cos(2.5 * np.pi * x2) + 2


def df6(t, Y, ud1, ud2):
    m1, m2, k1, k2, F = 5, 5, 1, 1, 1
    b1, b2 = ud2[0], ud2[1]
    dY = np.zeros(4)
    dY[0] = Y[1]
    dY[1] = (-b1 * Y[1] - b2 * (Y[1] - Y[3]) - k1 * Y[0] - k2 * (Y[0] - Y[2])) / m1
    dY[2] = Y[3]
    dY[3] = (F + b2 * (Y[1] - Y[3]) + k2 * (Y[0] - Y[2])) / m2
    return dY


def ff6R(x, ud1, ud2):
    N = 1001
    try:
        # Wczytaj dane z pliku 'polozenia.txt' z separatorem ';'
        X = np.genfromtxt("polozenia.txt", delimiter=';', comments=None)


    Y0 = np.array([0, 0, 0, 0])
    sol = solve_ivp(df6, [0, 10], Y0, args=(ud1, ud2), t_eval=np.linspace(0, 10, N))
    Y = sol.y.T  # Transpose to get the correct shape

    y = 0
    for i in range(N):
        print(f"{i} {Y[i, 0]} {Y[i, 2]}")
        y += abs(X[i, 0] - Y[i, 0]) + abs(X[i, 1] - Y[i, 2])

    y = y / (2.0 * N)
    return y
