from opt_alg import *
from user_fun import *
import numpy as np
import csv


def lab1_t():
    result = {"exp": [], "fib": [], "lag": []}

    for i in range(100):
        a, b, x0, fcalls = expansion(ff1T, np.random.uniform(-100, 100), 0.1, 3, 100)
        result["exp"].append([a, b, x0, fcalls])
        a, b = -100, 100
        result["fib"].append(fib(ff1T, a, b, 0.001))
        result["lag"].append(lag(ff1T, a, b, 0.001, 0.001, 100))

    # Zapis do pliku CSV
    with open("wynik_bez_zawezenia.csv", "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        # Zapisujemy nagłówki

        writer.writerow(
            [
                "Index",
                "Exp x0",
                "Exp a",
                "Exp b",
                "Exp fcalls",
                "Fib x",
                "Fib y",
                "Fib fcalls",
                "Lag x",
                "Lag y",
                "Lag fcalls",
            ]
        )

        # Zapisujemy dane
        for i in range(100):
            writer.writerow(
                [
                    i,
                    result["exp"][i][2],
                    result["exp"][i][0],
                    result["exp"][i][1],
                    result["exp"][i][3],
                    result["fib"][i][0],
                    result["fib"][i][1],
                    result["fib"][i][2],
                    result["lag"][i][0],
                    result["lag"][i][1],
                    result["lag"][i][2],
                ]
            )


def lab1_r():
    # Define parameters for the tank system
    P_A = 0.5
    V_0A = 5
    T_0A = 90
    P_B = 1
    V_0B = 1
    T_0B = 20
    F_in_B = 0.01
    T_in_B = 20
    D_B = 0.00365665
    a = 0.98
    b = 0.63
    g = 9.81

    # Parameter list for the differential equation function
    params = [P_A, V_0A, T_0A, P_B, V_0B, T_0B, F_in_B, T_in_B, D_B, a, b, g]

    D_A_interval = [0.0001, 0.01]

    print("Interval for D_A after expansion search:", D_A_interval)

    # Step 2: Apply Fibonacci method within the interval
    epsilon = 0.00001
    gamma = 0.00001
    optimal_D_A = fib_r(df1, D_A_interval[0], D_A_interval[1], epsilon, params)

    print("Optimal D_A:", optimal_D_A)
    # optimal_D_A = 0.001166
    max_temp = df1(optimal_D_A, params)
    print(f"Maximum Temperature in Tank B at optimal D_A: {max_temp + 50} °C")

    optimal_D_A = lag_r(
        df1, D_A_interval[0], D_A_interval[1], epsilon, gamma, 1000, params
    )
    # print("Optimal D_A:", optimal_D_A)
    # # optimal_D_A = 0.001166
    # max_temp = max_temp_in_tank(optimal_D_A, params)
    # print(f"Maximum Temperature in Tank B at optimal D_A: {max_temp+50} °C")
    #


def lab2_t():
    result = {"x": [], "jeeves": [], "rosenbrock": []}

    for i in range(100):
        x = np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1)])
        result["x"].append(x)
        result["jeeves"].append(hooke_jeeves(ff2T, x, 0.001, 0.3, 0.0001, 1000))
        result["rosenbrock"].append(rosenbrock_method(ff2T, x, 0.001, 1.1, 0.5, 0.0001, 1000))

    # Zapis do pliku CSV
    with open("wyniki_lab_2.csv", "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        # Zapisujemy nagłówki

        writer.writerow(
            [
                "Index",
                "x1",
                "x2",
                "HJ x1",
                "HJ x2",
                "HJ y",
                "HJ fcalls",
                "R x1",
                "R x2",
                "R y",
                "R fcalls",
            ]
        )

        # Zapisujemy dane
        for i in range(100):
            writer.writerow(
                [
                    i,
                    result["x"][i][0],
                    result["x"][i][1],
                    result["jeeves"][i][0],
                    result["jeeves"][i][1],
                    result["jeeves"][i][2],
                    result["jeeves"][i][3],
                    result["rosenbrock"][i][0],
                    result["rosenbrock"][i][1],
                    result["rosenbrock"][i][2],
                    result["rosenbrock"][i][3],
                ]
            )


def lab2_r():
    alfa_HJ = 0.1
    alfa_Ros = 2
    epsilon = 0.00001
    beta = 0.5
    Nmax = 1000
    s = 0.1
    x0 = np.array([1.0, 1.0])

    result_hj = hooke_jeeves(ff2R, x0, s, alfa_HJ, epsilon, Nmax)
    print("Hooke-Jeeves result:", result_hj)

    result_ros = rosenbrock_method(ff2R, x0, s, alfa_Ros, beta, epsilon, Nmax)
    print("Rosenbrock result:", result_ros)

    print("Test k1=5, k2=5:", ff2R([5, 5]))  # Bazowy test


def lab3_t():
    # Punkt startowy
    x0 = np.array([3.0, 3.0])

    # Parametry metody
    c1 = 1.0
    alpha = 10.0
    epsilon = 1e-6
    max_calls = 100000

    # Lista ograniczeń
    constraints = [constraint_g1, constraint_g2, constraint_g3]

    # Uruchomienie metody kary
    try:
        x_star = penalty_method(ff3t, constraints, x0, c1, alpha, epsilon, max_calls)
        print("Optymalny punkt:", x_star)
        print("Norma wyniku:", np.linalg.norm(x_star))
    except ValueError as e:
        print("Błąd:", e)


def lab4_t():
    x0 = np.array([-5.0, -5.0])
    h = 0.12
    print("Metoda gradientów prostych:")
    res = metoda_gradientow_prostych(ff4t, x0, 0.001,  10000,h)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", ff4t(res))

    print("Metoda gradientów sprzężonych:")
    res = metoda_gradientow_sprzezonych(ff4t, x0, 0.001,  10000,h)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", ff4t(res))

    print("Metoda Newtona:")
    res = metoda_newtona(ff4t, x0, 0.001, 1000, h)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", ff4t(res))

    print("Metoda gradientow prostych zmiennokrokowa:")
    res = metoda_gradientow_prostych_zmiennoskokowa(ff4t, x0, 0.001,  10000,0,10)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", ff4t(res))

    print("Metoda gradientow sprzezonych zmiennokrokowa:")
    res = metoda_gradientow_sprzezonych_zmiennoskokowa(ff4t, x0, 0.001,  10000,0,10)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", ff4t(res))

    print("Metoda Newtona zmiennokrokowa:")
    res = metoda_newtona_zmiennoskokowa(ff4t, x0, 0.001,  10000,0,10)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", ff4t(res))


def lab4_r():
    # Wczytanie danych
    X, y = load_data("XData.txt", "YData.txt")

    # Inicjalizacja parametrów
    initial_theta = np.zeros(3)

    # Obliczenie początkowego kosztu i gradientu
    initial_cost = compute_cost(initial_theta, X, y)
    print("Początkowy koszt: ", initial_cost)
    initial_gradient = compute_gradient(initial_theta, X, y)
    print("Początkowy gradient: ", initial_gradient)

    res = metoda_gradientow_sprzezonych_r(initial_theta, X, y, 0.001, 10000, 0.0001)
    print(res)
    print("Wartość funkcji celu w punkcie optymalnym: ", compute_cost(res, X, y))
    print(ocena_klasyfikatora(res, X, y))

def lab5_t():
    a = 10  # Parameter for the objective functions
    x0 = np.array([-10.0, -10.0])  # Initial guess
    epsilon = 1e-6
    n_max = 1000

    # Optimize aggregated function
    result = powell_method(lambda x: aggregated_function(x, a, w1=0.01, w2=0.99), x0, epsilon, n_max)
    print("Result (Pareto optimization):", result)


if __name__ == "__main__":
    lab5_t()
