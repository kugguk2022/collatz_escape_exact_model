import matplotlib.pyplot as plt

def syracuse(n: int) -> int:
    return n//2 if (n % 2 == 0) else (3*n + 1)//2

def orbit(n0: int, steps: int):
    n = n0
    out = [n]
    for _ in range(steps):
        n = syracuse(n)
        out.append(n)
    return out

def plot_2adic_ball_index(n0: int, steps: int = 500, k: int = 12):
    mod = 2**k
    seq = orbit(n0, steps)
    idx = [x % mod for x in seq]  # which ultrametric ball B(a,2^-k)
    plt.figure()
    plt.plot(idx, linewidth=1)
    plt.title(f"Syracuse orbit in Z2: ball index x_n mod 2^{k} (n0={n0})")
    plt.xlabel("n")
    plt.ylabel(f"x_n mod 2^{k}")
    plt.show()

plot_2adic_ball_index(27, steps=800, k=14)