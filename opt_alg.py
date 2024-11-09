import numpy as np
import csv


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
    with open("Jeeves.csv","w",newline="")as f:
        writer = csv.writer(f)
        while True:

            xB = x
            writer.writerow(xB)
            x,fcalls = proba(ff, xB, s, e,fcalls)
            fcalls += 2
            if ff(x) < ff(xB):
                while True:
                    i += 1
                    xB_prev = xB
                    xB = x

                    x = 2 * xB - xB_prev
                    x,fcalls = proba(ff, x, s, e,fcalls)
                    if i > nmax:
                        raise Exception("przekroczono liczbe maksymalnych wywolan funkcji")
                    fcalls+=2
                    if ff(x) >= ff(xB):
                        break
                x = xB
            else:
                s = alfa * s

            if i > nmax:
                raise Exception("przekroczono liczbe maksymalnych wywolan funkcji")

            if s < epsilon:
                break

        return [xB[0],xB[1],ff(xB),fcalls]


def proba(ff, x, s, e,f):
    for j in range(2):
        f+=2
        if ff(x + s * e[j]) < ff(x):
            x = x + s * e[j]
            f+=2
        elif ff(x - s * e[j]) < ff(x):
            x = x - s * e[j]
    return x,f


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
    with open("Rosenbrock_excel.csv","w",newline="") as f:
        writer = csv.writer(f)
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
            writer.writerow(x_best)

    return [x_best[0],x_best[1],ff(x_best),f_calls]


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
