import math


class Complex:
    def __init__(self, real, imag=0):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        return Complex(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other):
        return Complex(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other):
        return Complex(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real,
        )

    def __truediv__(self, other):
        denom = other.real**2 + other.imag**2
        return Complex(
            (self.real * other.real + self.imag * other.imag) / denom,
            (self.imag * other.real - self.real * other.imag) / denom,
        )

    def magnitude(self):
        return math.sqrt(self.real**2 + self.imag**2)

    def phase(self):
        return math.atan2(self.imag, self.real)

    def __repr__(self):
        return f"{self.real:.4f} + {self.imag:.4f}i"


class Matrix:
    def __init__(self, data):
        self.data = data
        self.rows = len(data)
        self.cols = len(data[0]) if self.rows > 0 else 0

    def __add__(self, other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrix dimensions must match")
        result = [
            [self.data[i][j] + other.data[i][j] for j in range(self.cols)]
            for i in range(self.rows)
        ]
        return Matrix(result)

    def __sub__(self, other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrix dimensions must match")
        result = [
            [self.data[i][j] - other.data[i][j] for j in range(self.cols)]
            for i in range(self.rows)
        ]
        return Matrix(result)

    def __mul__(self, other):
        if self.cols != other.rows:
            raise ValueError("Invalid matrix multiplication")
        result = [[0 for _ in range(other.cols)] for _ in range(self.rows)]
        for i in range(self.rows):
            for j in range(other.cols):
                for k in range(self.cols):
                    result[i][j] += self.data[i][k] * other.data[k][j]
        return Matrix(result)

    def transpose(self):
        result = [[self.data[j][i] for j in range(self.rows)] for i in range(self.cols)]
        return Matrix(result)

    def determinant(self):
        if self.rows != self.cols:
            raise ValueError("Matrix must be square")
        if self.rows == 1:
            return self.data[0][0]
        if self.rows == 2:
            return self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0]

        det = 0
        for j in range(self.cols):
            minor = [
                [self.data[i][k] for k in range(self.cols) if k != j]
                for i in range(1, self.rows)
            ]
            minor_matrix = Matrix(minor)
            sign = -1 if j % 2 else 1
            det += sign * self.data[0][j] * minor_matrix.determinant()
        return det

    def inverse(self):
        det = self.determinant()
        if det == 0:
            raise ValueError("Matrix is singular")
        if self.rows == 2:
            a, b, c, d = (
                self.data[0][0],
                self.data[0][1],
                self.data[1][0],
                self.data[1][1],
            )
            return Matrix([[d / det, -b / det], [-c / det, a / det]])

        adj = [[0 for _ in range(self.cols)] for _ in range(self.rows)]
        for i in range(self.rows):
            for j in range(self.cols):
                minor = [
                    [self.data[r][c] for c in range(self.cols) if c != j]
                    for r in range(self.rows)
                    if r != i
                ]
                minor_matrix = Matrix(minor)
                adj[j][i] = ((-1) ** (i + j)) * minor_matrix.determinant()

        return Matrix(
            [[adj[i][j] / det for j in range(self.cols)] for i in range(self.rows)]
        )

    def eigenvalues(self):
        if self.rows != 2 or self.cols != 2:
            raise ValueError("Only 2x2 matrix eigenvalue computation implemented")
        a, b, c, d = self.data[0][0], self.data[0][1], self.data[1][0], self.data[1][1]
        trace = a + d
        det = a * d - b * c
        discriminant = trace**2 - 4 * det

        if discriminant < 0:
            real_part = trace / 2
            imag_part = math.sqrt(-discriminant) / 2
            return [Complex(real_part, imag_part), Complex(real_part, -imag_part)]
        else:
            return [
                (trace + math.sqrt(discriminant)) / 2,
                (trace - math.sqrt(discriminant)) / 2,
            ]

    def __repr__(self):
        rows = ["[" + ", ".join([f"{x:8.4f}" for x in row]) + "]" for row in self.data]
        return "\n".join(rows)


class NumericalAnalysis:
    @staticmethod
    def newton_raphson(f, df, x0, tolerance=1e-10, max_iter=1000):
        x = x0
        for _ in range(max_iter):
            fx = f(x)
            if abs(fx) < tolerance:
                return x
            dfx = df(x)
            if dfx == 0:
                raise ValueError("Derivative is zero")
            x = x - fx / dfx
        return x

    @staticmethod
    def bisection(f, a, b, tolerance=1e-10, max_iter=1000):
        fa, fb = f(a), f(b)
        if fa * fb > 0:
            raise ValueError("f(a) and f(b) must have opposite signs")

        for _ in range(max_iter):
            mid = (a + b) / 2
            fmid = f(mid)
            if abs(fmid) < tolerance:
                return mid
            if fa * fmid < 0:
                b = mid
                fb = fmid
            else:
                a = mid
                fa = fmid
        return (a + b) / 2

    @staticmethod
    def simpson(f, a, b, n=1000):
        if n % 2 == 1:
            n += 1
        h = (b - a) / n
        result = f(a) + f(b)

        for i in range(1, n):
            x = a + i * h
            if i % 2 == 0:
                result += 2 * f(x)
            else:
                result += 4 * f(x)

        return result * h / 3

    @staticmethod
    def gaussian_quadrature(f, n=5):
        if n == 5:
            points = [
                -0.906179845938664,
                -0.538469310105683,
                0,
                0.538469310105683,
                0.906179845938664,
            ]
            weights = [
                0.236926885056189,
                0.478628670499366,
                0.568888888888889,
                0.478628670499366,
                0.236926885056189,
            ]
        elif n == 3:
            points = [-0.774596669210483, 0, 0.774596669210483]
            weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
        else:
            raise ValueError("Only n=3 or n=5 implemented")

        return sum(w * f(p) for p, w in zip(points, weights))

    @staticmethod
    def romberg(f, a, b, n=5):
        r = [[0 for _ in range(n)] for _ in range(n)]

        for i in range(n):
            h = (b - a) / (2**i)
            sub_sum = f(a) + f(b)

            for k in range(1, 2**i):
                if k % 2 == 1:
                    sub_sum += 2 * f(a + k * h)

            r[i][0] = sub_sum * h / 3

        for j in range(1, n):
            for i in range(j, n):
                r[i][j] = (4**j * r[i][j - 1] - r[i - 1][j - 1]) / (4**j - 1)

        return r[n - 1][n - 1]


class Statistics:
    @staticmethod
    def mean(data):
        return sum(data) / len(data)

    @staticmethod
    def variance(data):
        m = Statistics.mean(data)
        return sum((x - m) ** 2 for x in data) / len(data)

    @staticmethod
    def std_dev(data):
        return math.sqrt(Statistics.variance(data))

    @staticmethod
    def covariance(x, y):
        if len(x) != len(y):
            raise ValueError("Data arrays must have same length")
        mean_x, mean_y = Statistics.mean(x), Statistics.mean(y)
        return sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(len(x))) / len(x)

    @staticmethod
    def correlation(x, y):
        cov = Statistics.covariance(x, y)
        std_x, std_y = Statistics.std_dev(x), Statistics.std_dev(y)
        return cov / (std_x * std_y)

    @staticmethod
    def linear_regression(x, y):
        n = len(x)
        sum_x, sum_y = sum(x), sum(y)
        sum_xy = sum(x[i] * y[i] for i in range(n))
        sum_x2 = sum(xi**2 for xi in x)

        slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)
        intercept = (sum_y - slope * sum_x) / n

        return slope, intercept

    @staticmethod
    def percentile(data, p):
        sorted_data = sorted(data)
        k = (len(sorted_data) - 1) * p / 100
        f = int(k)
        c = f + 1 if f + 1 < len(sorted_data) else f
        d = k - f
        return sorted_data[f] + d * (sorted_data[c] - sorted_data[f])


class SpecialFunctions:
    @staticmethod
    def sin(x, terms=20):
        x = x % (2 * math.pi)
        result = 0
        for n in range(terms):
            result += ((-1) ** n * x ** (2 * n + 1)) / SpecialFunctions.factorial(
                2 * n + 1
            )
        return result

    @staticmethod
    def cos(x, terms=20):
        x = x % (2 * math.pi)
        result = 0
        for n in range(terms):
            result += ((-1) ** n * x ** (2 * n)) / SpecialFunctions.factorial(2 * n)
        return result

    @staticmethod
    def exp(x, terms=30):
        result = 0
        for n in range(terms):
            result += x**n / SpecialFunctions.factorial(n)
        return result

    @staticmethod
    def ln(x, terms=100):
        if x <= 0:
            raise ValueError(" logarithm argument must be positive")

        y = (x - 1) / (x + 1)
        result = 0
        for n in range(terms):
            result += (y ** (2 * n + 1)) / (2 * n + 1)
        return 2 * result

    @staticmethod
    def sqrt(x, terms=20):
        if x < 0:
            raise ValueError("Cannot compute square root of negative number")

        guess = x / 2
        for _ in range(terms):
            guess = (guess + x / guess) / 2
        return guess

    @staticmethod
    def factorial(n):
        if n < 0:
            raise ValueError("Factorial of negative number undefined")
        if n == 0 or n == 1:
            return 1
        result = 1
        for i in range(2, int(n) + 1):
            result *= i
        return result

    @staticmethod
    def gamma(z, terms=50):
        if z <= 0:
            raise ValueError("Gamma function not defined for non-positive values")

        if z < 1:
            return math.pi / (
                SpecialFunctions.gamma(1 - z) * SpecialFunctions.sin(math.pi * z)
            )

        result = 1
        for n in range(1, terms):
            result *= (1 + 1 / n) ** (1 - z) / (1 + z / n)

        return result

    @staticmethod
    def bessel_j(x, n, terms=20):
        result = 0
        for k in range(terms):
            result += ((-1) ** k * (x / 2) ** (2 * k + n)) / (
                SpecialFunctions.factorial(k) * SpecialFunctions.factorial(k + n)
            )
        return result

    @staticmethod
    def legendre(x, n):
        if n == 0:
            return 1
        if n == 1:
            return x

        p0, p1 = 1, x
        for i in range(2, n + 1):
            p2 = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i
            p0, p1 = p1, p2
        return p1


class FourierTransform:
    @staticmethod
    def dft(signal):
        n = len(signal)
        result = [Complex(0, 0) for _ in range(n)]

        for k in range(n):
            for t in range(n):
                angle = -2 * math.pi * k * t / n
                result[k] = result[k] + Complex(
                    signal[t] * math.cos(angle), signal[t] * math.sin(angle)
                )

        return result

    @staticmethod
    def idft(spectrum):
        n = len(spectrum)
        result = [0 for _ in range(n)]

        for t in range(n):
            for k in range(n):
                angle = 2 * math.pi * k * t / n
                result[t] += spectrum[k].real * math.cos(angle) - spectrum[
                    k
                ].imag * math.sin(angle)

        return [x / n for x in result]


class Optimizer:
    @staticmethod
    def gradient_descent(
        grad_f, x0, learning_rate=0.01, tolerance=1e-10, max_iter=10000
    ):
        x = x0
        for _ in range(max_iter):
            gradient = grad_f(x)
            x_new = x - learning_rate * gradient
            if abs(x_new - x) < tolerance:
                return x_new
            x = x_new
        return x

    @staticmethod
    def nelder_mead(f, vertices, tolerance=1e-10, max_iter=1000):
        n = len(vertices)
        vertices = sorted(vertices, key=f)

        for _ in range(max_iter):
            centroid = [sum(vertices[i][j] for i in range(n)) / n for j in range(n - 1)]

            worst = vertices[n - 1]
            reflection = [centroid[j] + (centroid[j] - worst[j]) for j in range(n - 1)]
            f_reflection = f(reflection)

            best, second_worst = vertices[0], vertices[n - 2]
            expansion = [
                centroid[j] + 2 * (reflection[j] - centroid[j]) for j in range(n - 1)
            ]
            f_expansion = f(expansion)

            if f_reflection < second_worst[1]:
                if f_reflection < best[1]:
                    vertices[n - 1] = (reflection, f_reflection)
                else:
                    vertices[n - 1] = (expansion, f_expansion)
            elif f_reflection < worst[1]:
                vertices[n - 1] = (reflection, f_reflection)
                worst = vertices[n - 1]

            shrink = [
                vertices[0][j] + 0.5 * (vertices[i][j] - vertices[0][j])
                for i in range(n)
                for j in range(n - 1)
            ]
            continue

        return vertices[0][0]


class Polynomial:
    def __init__(self, coefficients):
        self.coeffs = coefficients

    def evaluate(self, x):
        result = 0
        for i, c in enumerate(self.coeffs):
            result += c * (x**i)
        return result

    def derivative(self):
        return Polynomial([c * i for i, c in enumerate(self.coeffs) if i > 0])

    def integrate(self):
        return Polynomial([0] + [c / i for i, c in enumerate(self.coeffs, 1)])

    def roots(self):
        if len(self.coeffs) == 2:
            if self.coeffs[1] == 0:
                raise ValueError("No solution")
            return [-self.coeffs[0] / self.coeffs[1]]
        return []


def main():
    print("=== Complex Numbers ===")
    c1 = Complex(3, 4)
    c2 = Complex(1, -2)
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c1 + c2 = {c1 + c2}")
    print(f"c1 * c2 = {c1 * c2}")
    print(f"|c1| = {c1.magnitude()}")

    print("\n=== Matrix Operations ===")
    m1 = Matrix([[1, 2], [3, 4]])
    m2 = Matrix([[5, 6], [7, 8]])
    print("m1 =")
    print(m1)
    print("det(m1) =", m1.determinant())
    print("m1^-1 =")
    print(m1.inverse())
    print("eigenvalues(m1) =", m1.eigenvalues())

    print("\n=== Numerical Analysis ===")
    f = lambda x: x**3 - x - 2
    df = lambda x: 3 * x**2 - 1
    root = NumericalAnalysis.newton_raphson(f, df, 2)
    print(f"Root of x^3 - x - 2 = {root}")
    print(f"f(root) = {f(root)}")

    print("\n=== Integrals ===")
    f = lambda x: 1 / (1 + x**2)
    result = NumericalAnalysis.romberg(f, 0, 1, 5)
    print(f"Integrate 0->1 (1/(1+x^2)) dx = {result}")
    print(f"Actual: pi/4 = {math.pi / 4}")

    print("\n=== Statistics ===")
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [2.1, 4.2, 6.1, 8.2, 10.1, 12.2, 14.1, 16.2, 18.1, 20.2]
    print(f"Mean: {Statistics.mean(x)}")
    print(f"Std Dev: {Statistics.std_dev(x)}")
    slope, intercept = Statistics.linear_regression(x, y)
    print(f"Linear regression: y = {slope:.4f}x + {intercept:.4f}")
    print(f"Correlation: {Statistics.correlation(x, y)}")

    print("\n=== Special Functions ===")
    print(f"sin(pi/6) = {SpecialFunctions.sin(math.pi / 6):.6f} (actual: 0.5)")
    print(f"cos(pi/3) = {SpecialFunctions.cos(math.pi / 3):.6f} (actual: 0.5)")
    print(f"e^1 = {SpecialFunctions.exp(1):.6f} (actual: {math.e:.6f})")
    print(f"sqrt(2) = {SpecialFunctions.sqrt(2):.6f}")
    print(f"Jo(1) (Bessel) = {SpecialFunctions.bessel_j(1, 0):.6f}")

    print("\n=== Fourier Transform ===")
    signal = [1, 2, 3, 4]
    spectrum = FourierTransform.dft(signal)
    print("DFT of [1, 2, 3, 4]:")
    for i, c in enumerate(spectrum):
        print(f"  F[{i}] = {c}")

    print("\n=== Optimization ===")
    f = lambda x: x**4 - 3 * x**3 + 2
    df = lambda x: 4 * x**3 - 9 * x**2
    minimum = Optimizer.gradient_descent(df, 3, learning_rate=0.001)
    print(f"Minimum of x^4 - 3x^3 + 2 at x = {minimum}")
    print(f"f(minimum) = {f(minimum)}")

    p = Polynomial([2, -3, 1])
    deriv = p.derivative()
    print(f"\nPolynomial 2 - 3x + x^2")
    print(f"Derivative: {deriv.coeffs}")


if __name__ == "__main__":
    main()
