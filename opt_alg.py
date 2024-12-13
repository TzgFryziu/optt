import numpy as np


def expansion(ff, x0, d, alpha, nmax):
    x = [x0, x0 + d]
    i = 0
    f = 0
    if ff(x[0]) == ff(x[1]):
        f += 2
        return x[0], x[1], x0, f
    f += 2
    if ff(x[1]) > ff(x[0]):
        f += 2
        d = -d
        x[1] = x[0] + d
        if ff(x[1]) >= ff(x[0]):
            f += 2
            return x[1], x[0] - d, x0, f
    f += 4
    while True:
        if i >= nmax:
            raise Exception("Exceeded maximum number of iterations")
        i += 1
        x.append(x[0] + (alpha ** i) * d)

        if ff(x[i]) <= ff(x[i + 1]):
            f += 2
            break
        f += 2
    if d > 0:
        return x[i - 1], x[i + 1], x0, f
    else:
        return x[i + 1], x[i - 1], x0, f


def fib_seq(n):
    ksi = (1 + np.sqrt(5)) / 2
    return (ksi ** n - (1 - ksi) ** n) / np.sqrt(5)


def fib(ff, a, b, epsilon):
    i = 0

    k = 0
    while fib_seq(k) < (b - a) / epsilon:
        k += 1

    a = [a]
    b = [b]
    c = [b[0] - fib_seq(k - 1) / fib_seq(k) * (b[0] - a[0])]
    d = [a[0] + fib_seq(k - 1) / fib_seq(k) * (b[0] - a[0])]

    for i in range(k - 2):
        if ff(c[i]) < ff(d[i]):
            a.append(a[i])
            b.append(d[i])
            d.append(c[i])
            c.append(
                b[i + 1]
                - fib_seq(k - i - 3) / fib_seq(k - i - 2) * (b[i + 1] - a[i + 1])
            )
        else:
            a.append(c[i])
            b.append(b[i])
            c.append(d[i])
            d.append(
                a[i + 1]
                + fib_seq(k - i - 3) / fib_seq(k - i - 2) * (b[i + 1] - a[i + 1])
            )
        i += 2
        if abs(b[-1] - a[-1]) < epsilon:
            break

    return c[-1], ff(c[-1]), i


def lag(ff, a, b, epsilon, gamma, nmax):
    i = 0
    f = 0
    a = [a]
    b = [b]
    c = [(a[0] + b[0]) / 2]  # Initial midpoint
    d = []

    while True:
        print(f"itteration {i}")
        # Calculate the numerator (l) and denominator (m) for d[i]
        l = (
                ff(a[i]) * (b[i] ** 2 - c[i] ** 2)
                + ff(b[i]) * (c[i] ** 2 - a[i] ** 2)
                + ff(c[i]) * (a[i] ** 2 - b[i] ** 2)
        )

        m = (
                ff(a[i]) * (b[i] - c[i])
                + ff(b[i]) * (c[i] - a[i])
                + ff(c[i]) * (a[i] - b[i])
        )
        f += 6
        if m == 0:  # Handle division by zero
            raise Exception("Division by zero error in m calculation")

        # Compute the new point d[i]
        d_new = 0.5 * l / m
        d.append(d_new)

        # Ensure d_new is within (a[i], b[i])
        if not (a[i] < d_new < b[i]):
            return "Error", "Error", f

        # Determine new interval based on d_new
        if a[i] < d_new < c[i]:
            f += 2
            if ff(d_new) < ff(c[i]):
                a.append(a[i])
                b.append(c[i])
                c.append(d_new)
            else:
                a.append(d_new)
                b.append(b[i])
                c.append(c[i])
        elif c[i] < d_new < b[i]:
            f += 2
            if ff(d_new) < ff(c[i]):
                a.append(c[i])
                b.append(b[i])
                c.append(d_new)
            else:
                a.append(a[i])
                b.append(d_new)
                c.append(c[i])

        # Update midpoint `c[i]` to be the new midpoint for accuracy
        c[-1] = (a[-1] + b[-1]) / 2

        # Update iteration count
        i += 1

        # Stopping conditions with adjusted precision
        if i >= nmax:
            raise Exception("Exceeded maximum number of iterations")

        # Adjust stopping conditions to check for closer convergence
        if i > 5:  # Require at least 5 iterations before checking convergence
            if (b[i] - a[i]) < epsilon:
                break
            if len(d) > 2 and abs(d[-1] - d[-2]) < gamma and abs(d[-2] - d[-3]) < gamma:
                break
    print(f"Lag fcalls {f}")

    return (
        d[-1],
        ff(d[-1]),
        f,
    )  # Return the last d value, which is the approximate solution


def fib_r(ff, a, b, epsilon, params):
    f = 0
    k = 0
    while fib_seq(k) < (b - a) / epsilon:
        k += 1

    a = [a]
    b = [b]
    c = [b[0] - fib_seq(k - 1) / fib_seq(k) * (b[0] - a[0])]
    d = [a[0] + fib_seq(k - 1) / fib_seq(k) * (b[0] - a[0])]

    for i in range(k - 2):
        f += 2
        if ff(c[i], params) < ff(d[i], params):
            a.append(a[i])
            b.append(d[i])
            d.append(c[i])
            c.append(
                b[i + 1]
                - fib_seq(k - i - 3) / fib_seq(k - i - 2) * (b[i + 1] - a[i + 1])
            )
        else:
            a.append(c[i])
            b.append(b[i])
            c.append(d[i])
            d.append(
                a[i + 1]
                + fib_seq(k - i - 3) / fib_seq(k - i - 2) * (b[i + 1] - a[i + 1])
            )

        if abs(b[-1] - a[-1]) < epsilon:
            break

    print(f"Fib fcalls {f}")
    return c[-1]


def lag_r(ff, a, b, epsilon, gamma, nmax, params):
    i = 0
    f = 0
    a = [a]
    b = [b]
    c = [(a[0] + b[0]) / 2]  # Initial midpoint
    d = []

    while True:
        print(f"itteration {i}")
        # Calculate the numerator (l) and denominator (m) for d[i]
        l = (
                ff(a[i], params) * (b[i] ** 2 - c[i] ** 2)
                + ff(b[i], params) * (c[i] ** 2 - a[i] ** 2)
                + ff(c[i], params) * (a[i] ** 2 - b[i] ** 2)
        )

        m = (
                ff(a[i], params) * (b[i] - c[i])
                + ff(b[i], params) * (c[i] - a[i])
                + ff(c[i], params) * (a[i] - b[i])
        )
        f += 6
        if m == 0:  # Handle division by zero
            raise Exception("Division by zero error in m calculation")

        # Compute the new point d[i]
        d_new = 0.5 * l / m
        d.append(d_new)

        # Ensure d_new is within (a[i], b[i])
        if not (a[i] < d_new < b[i]):
            raise Exception("d[i] is out of bounds")

        # Determine new interval based on d_new
        if a[i] < d_new < c[i]:
            f += 2
            if ff(d_new, params) < ff(c[i], params):
                a.append(a[i])
                b.append(c[i])
                c.append(d_new)
            else:
                a.append(d_new)
                b.append(b[i])
                c.append(c[i])
        elif c[i] < d_new < b[i]:
            f += 2
            if ff(d_new, params) < ff(c[i], params):
                a.append(c[i])
                b.append(b[i])
                c.append(d_new)
            else:
                a.append(a[i])
                b.append(d_new)
                c.append(c[i])

        # Update midpoint `c[i]` to be the new midpoint for accuracy
        c[-1] = (a[-1] + b[-1]) / 2

        # Update iteration count
        i += 1

        # Stopping conditions with adjusted precision
        if i >= nmax:
            raise Exception("Exceeded maximum number of iterations")

        # Adjust stopping conditions to check for closer convergence
        if i > 5:  # Require at least 5 iterations before checking convergence
            if (b[i] - a[i]) < epsilon:
                break
            if len(d) > 2 and abs(d[-1] - d[-2]) < gamma and abs(d[-2] - d[-3]) < gamma:
                break
    print(f"Lag fcalls {f}")
    return d[-1]  # Return the last d value, which is the approximate solution


def hooke_jeeves(ff, x, s, alfa, epsilon, nmax):
    # Kierunki jako np.array
    e = np.array([[1, 0], [0, 1]])

    i = 0
    fcalls = 0
    while True:
        xB = x
        x, fcalls = proba(ff, xB, s, e, fcalls)
        fcalls += 2
        if ff(x) < ff(xB):
            while True:
                i += 1
                xB_prev = xB
                xB = x

                x = 2 * xB - xB_prev
                x, fcalls = proba(ff, x, s, e, fcalls)
                if i > nmax:
                    raise Exception("przekroczono liczbe maksymalnych wywolan funkcji")
                fcalls += 2
                if ff(x) >= ff(xB):
                    break
            x = xB
        else:
            s = alfa * s

        if i > nmax:
            raise Exception("przekroczono liczbe maksymalnych wywolan funkcji")

        if s < epsilon:
            break

    return [xB[0], xB[1], ff(xB), fcalls]


def proba(ff, x, s, e, f):
    for j in range(2):
        f += 2
        if ff(x + s * e[j]) < ff(x):
            x = x + s * e[j]
            f += 2
        elif ff(x - s * e[j]) < ff(x):
            x = x - s * e[j]
    return x, f


def rosenbrock_method(
        ff, x0, s0, alpha, beta, epsilon, Nmax
):
    # Initialize variables
    n = len(x0)
    i = 0
    d = np.eye(n)  # d_j^(0) = e_j (identity matrix for initial directions)
    lam = np.zeros(n)  # λ_j^(0) = 0
    p = np.zeros(n)  # p_j^(0) = 0
    s = np.array([float(s0)] * n)  # Ensure `s` is an array
    x_best = x0.copy()
    f_calls = 0  # Track the number of function calls
    # Main loop
    while True:
        for j in range(n):
            # Check the condition for expansion step
            if ff(x_best + s[j] * d[:, j]) < ff(x_best):
                x_best = x_best + s[j] * d[:, j]
                lam[j] += s[j]
                s[j] *= alpha
            else:
                # Contraction step
                s[j] *= -beta
                p[j] += 1

            f_calls += 2
            if f_calls > Nmax:
                print("Exceeded maximum function evaluations")
                return x_best
        # Update x^(i+1)
        x0 = x_best.copy()
        i += 1

        # Check if direction base needs to be reset
        if all(l != 0 for l in lam) and all(x != 0 for x in p):
            # Reinitialize directions and steps
            d = update_directions(d, lam)
            lam = np.zeros(n)
            p = np.zeros(n)
            s = np.array([float(s0)] * n)  # Reset step sizes to initial values

        # Stopping criterion
        if max(abs(s)) < epsilon:
            break

    return [x_best[0], x_best[1], ff(x_best), f_calls]


def update_directions(d, lam):
    # Uaktualnienie bazy kierunków z wektorem `λ`
    n = len(lam)
    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, :i + 1] = lam[i]
    Q = np.dot(d, Q)

    # Proces Gram-Schmidta dla ortogonalizacji kolumn
    v = Q[:, 0] / np.linalg.norm(Q[:, 0])
    d[:, 0] = v

    for i in range(1, n):
        temp = np.zeros(n)
        for j in range(i):
            temp += np.dot(Q[:, i].T, d[:, j]) * d[:, j]
        v = Q[:, i] - temp
        v /= np.linalg.norm(v)
        d[:, i] = v

    return d


def nelder_mead(func, x0, s=1.0, alpha=1.0, gamma=2.0, beta=0.5, delta=0.5, epsilon=1e-6, max_calls=500):
    # Liczba wymiarów problemu
    n = len(x0)

    # Inicjalizacja sympleksu
    simplex = [x0]
    for i in range(n):
        vertex = np.copy(x0)
        vertex[i] += s
        simplex.append(vertex)

    # Funkcje celu dla wierzchołków
    f_values = [func(v) for v in simplex]
    calls = n + 1  # Liczymy wywołania funkcji celu

    while calls < max_calls:
        # Sortowanie sympleksu wg wartości funkcji celu
        indices = np.argsort(f_values)
        simplex = [simplex[i] for i in indices]
        f_values = [f_values[i] for i in indices]
        print(func(simplex[0]))
        # Sprawdzenie warunku stopu
        if np.max([np.linalg.norm(simplex[0] - v) for v in simplex]) < epsilon:
            return simplex[0]

        # Wyznaczenie środka ciężkości pomijając najgorszy wierzchołek
        p = np.mean(simplex[:-1], axis=0)

        # Odbicie
        p_reflect = p + alpha * (p - simplex[-1])
        f_reflect = func(p_reflect)
        calls += 1

        if f_values[0] <= f_reflect < f_values[-2]:
            # Akceptacja odbicia
            simplex[-1] = p_reflect
            f_values[-1] = f_reflect
            continue

        if f_reflect < f_values[0]:
            # Ekspansja
            p_expand = p + gamma * (p_reflect - p)
            f_expand = func(p_expand)
            calls += 1

            if f_expand < f_reflect:
                simplex[-1] = p_expand
                f_values[-1] = f_expand
            else:
                simplex[-1] = p_reflect
                f_values[-1] = f_reflect
            continue

        # Zawężenie
        p_contract = p + beta * (simplex[-1] - p)
        f_contract = func(p_contract)
        calls += 1

        if f_contract < f_values[-1]:
            simplex[-1] = p_contract
            f_values[-1] = f_contract
        else:
            # Redukcja
            for i in range(1, len(simplex)):
                simplex[i] = simplex[0] + delta * (simplex[i] - simplex[0])
                f_values[i] = func(simplex[i])
            calls += n

    # Jeśli osiągnięto maksymalną liczbę wywołań funkcji celu
    raise RuntimeError("Osiągnięto maksymalną liczbę wywołań funkcji celu bez zbieżności.")


def penalty_method(func, constraints, x0, c1, alpha, epsilon, max_calls):
    def augmented_function(x, c):
        # Rozszerzona funkcja celu z karą
        penalty = sum(max(0, g(x)) ** 2 for g in constraints)
        return func(x) + c * penalty

    x_current = x0
    c = c1
    calls = 0

    while True:
        # Optymalizacja funkcji z karą przy bieżącym c
        def penalized_func(x):
            nonlocal calls
            calls += 1
            return augmented_function(x, c)

        try:
            x_next = nelder_mead(penalized_func, x_current, epsilon=epsilon, max_calls=max_calls - calls)
        except RuntimeError:
            raise RuntimeError("Przekroczono maksymalną liczbę wywołań funkcji celu.")

        # Sprawdzenie warunku stopu
        if np.linalg.norm(x_next - x_current) < epsilon:
            return x_next

        # Aktualizacja współczynnika kary i punktu startowego
        x_current = x_next
        c *= alpha


def pochodna_1(x):
    x1, x2 = x
    return 10 * x1 + 8 * x2 - 34


def pochodna_2(x):
    x1, x2 = x
    return 10 * x2 + 8 * x1 - 38

def pochodna_1_1(x):
    x1, x2 = x
    return 10
def pochodna_1_2(x):
    x1, x2 = x
    return 8
def pochodna_2_1(x):
    x1, x2 = x
    return 8
def pochodna_2_2(x):
    x1, x2 = x
    return 10


def metoda_gradientow_prostych(ff,  x0, epsilon, nmax, h):

    i = 0
    x = np.array(x0, dtype=float)

    while True:
        # Compute gradient
        g = np.array([pochodna_1(x), pochodna_2(x)])
        # Update the current point
        direction = -g
        x_curr = x + h * direction
        i += 1

        # Check for convergence
        if np.linalg.norm(x_curr - x) < epsilon or i >= nmax:
            break

        x = x_curr

    return x



def metoda_gradientow_sprzezonych(ff, x0, epsilon, nmax, h):

    i = 0
    x = np.array(x0, dtype=float)
    g = np.array([pochodna_1(x), pochodna_2(x)])
    d = -g

    while True:
        # Update the current point
        x_curr = x + h * d
        i += 1

        # Check for convergence
        if np.linalg.norm(x_curr - x) < epsilon or i >= nmax:
            break

        # Compute new gradient
        g_curr = np.array([pochodna_1(x_curr), pochodna_2(x_curr)])

        # Compute beta using Fletcher-Reeves formula
        beta = np.linalg.norm(g_curr)**2 / np.linalg.norm(g)**2

        # Update direction
        d = -g_curr + beta * d

        # Update variables for the next iteration
        x = x_curr
        g = g_curr

    return x


def metoda_newtona(ff,  x0, epsilon, nmax, h):

    i = 0
    x = np.array(x0, dtype=float)

    while True:
        # Compute gradient and Hessian
        g = np.array([pochodna_1(x), pochodna_2(x)])
        H = np.array([[pochodna_1_1(x), pochodna_1_2(x)], [pochodna_2_1(x), pochodna_2_2(x)]])

        # Compute Newton direction d = -H^-1 * g
        try:
            d = -np.linalg.solve(H, g)  # More stable than np.linalg.inv
        except np.linalg.LinAlgError:
            print("Hessian is singular, using pseudo-inverse.")
            d = -np.dot(np.linalg.pinv(H), g)

        # Update the current point
        x_curr = x + h * d
        i += 1

        # Check for convergence
        if np.linalg.norm(x_curr - x) < epsilon or i >= nmax:
            break

        x = x_curr

    return x


def fibonacci_search(ff, a, b, tol):
    """
    Fibonacci method for step size selection.

    Parameters:
    - ff: The function to minimize.
    - a, b: The interval for the search.
    - tol: The tolerance for the search.

    Returns:
    - Optimal step size h.
    """
    fib = [1, 1]
    while (fib[-1] < (b - a) / tol):
        fib.append(fib[-1] + fib[-2])

    n = len(fib) - 1
    x1 = a + (fib[n - 2] / fib[n]) * (b - a)
    x2 = a + (fib[n - 1] / fib[n]) * (b - a)

    for _ in range(n - 1):
        if ff(x1) > ff(x2):
            a = x1
            x1 = x2
            x2 = a + (fib[n - 1] / fib[n]) * (b - a)
        else:
            b = x2
            x2 = x1
            x1 = a + (fib[n - 2] / fib[n]) * (b - a)

        n -= 1

    return (a + b) / 2

def metoda_gradientow_prostych_zmiennoskokowa(ff, x, epsilon, nmax, a=0, b=1, tol=1e-5):

    i = 0
    x = np.array(x, dtype=float)

    while True:
        # Compute gradient
        g = np.array([pochodna_1(x), pochodna_2(x)])

        # Determine step size using Fibonacci search
        direction = -g
        phi = lambda h: ff(x + h * direction)
        h = fibonacci_search(phi, a, b, tol)

        # Update the current point
        x_curr = x + h * direction
        i += 1

        # Check for convergence
        if np.linalg.norm(x_curr - x) < epsilon or i >= nmax:
            break

        x = x_curr

    return x

def metoda_gradientow_sprzezonych_zmiennoskokowa(ff, x, epsilon, nmax, a=0, b=1, tol=1e-5):
    i = 0
    x = np.array(x, dtype=float)
    g = [pochodna_1(x), pochodna_2(x)]
    d = -np.array([g[0], g[1]])
    while True:
        # Define the 1D function along the direction to find optimal step size
        phi = lambda h: ff(x + h * d)
        h = fibonacci_search(phi, 0, 1, 1e-5)  # Fibonacci search for step size

        # Update the current point
        x_curr = x + h * d
        i += 1

        # Check for convergence
        if np.linalg.norm(x_curr - x) < epsilon or i >= nmax:
            break

        # Compute new gradient

        g_curr = np.array([pochodna_1(x_curr), pochodna_2(x_curr)])


        # Compute beta using Fletcher-Reeves formula

        beta = np.linalg.norm(g_curr)**2 / np.linalg.norm(g)**2

        # Update direction
        d = -g_curr + beta * d

        # Update variables for the next iteration
        x = x_curr
        g = g_curr

    return x

def metoda_newtona_zmiennoskokowa(ff, x, epsilon, nmax, a=0, b=1, tol=1e-5):

    i = 0
    x = np.array(x, dtype=float)

    while True:
        # Compute gradient and Hessian
        g = np.array([pochodna_1(x), pochodna_2(x)])
        H = np.array([[pochodna_1_1(x), pochodna_1_2(x)], [pochodna_2_1(x), pochodna_2_2(x)]])

        # Debug: Print gradient norm
        grad_norm = np.linalg.norm(g)
        print(f"Iteration {i}: Gradient norm = {grad_norm}")

        if grad_norm < epsilon:
            print("Convergence achieved: Gradient norm is below threshold.")
            break

        # Compute Newton direction d = -H^-1 * g
        try:
            d = -np.linalg.solve(H, g)  # More stable than np.linalg.inv
        except np.linalg.LinAlgError:
            print("Hessian is singular, using pseudo-inverse.")
            d = -np.dot(np.linalg.pinv(H), g)

        # Determine step size using Fibonacci search
        phi = lambda h: ff(x + h * d)
        h = fibonacci_search(phi, a, b, tol)
        print(f"Iteration {i}: Step size = {h}")

        # Update the current point
        x_curr = x + h * d
        i += 1

        # Check for convergence
        if np.linalg.norm(x_curr - x) < epsilon:
            print("Convergence achieved: Change in x is below threshold.")
            break
        if i >= nmax:
            print("Maximum iterations reached.")
            break

        x = x_curr

    return x
