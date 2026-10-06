"""Models de plomes de contaminació."""
import numpy as np


def ploma_gaussiana(x, y, z, Q=1.0, u=2.0, H=10.0, sigma_y=None, sigma_z=None):
    """Concentració d'una ploma gaussiana estacionària.

    Font puntual a (0, 0, H); vent de velocitat u en la direcció x.
    Q: taxa d'emissió; sigma_y, sigma_z: funcions de x (m).
    Inclou reflexió del terra (z = 0).
    """
    x, y, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
    sigma_y = sigma_y or (lambda d: 0.22 * d / np.sqrt(1 + 0.0001 * d))
    sigma_z = sigma_z or (lambda d: 0.20 * d)
    c = np.zeros_like(x)
    m = x > 0
    sy, sz = sigma_y(x[m]), sigma_z(x[m])
    c[m] = (Q / (2 * np.pi * u * sy * sz)) * np.exp(-y[m] ** 2 / (2 * sy ** 2)) * (
        np.exp(-(z[m] - H) ** 2 / (2 * sz ** 2)) + np.exp(-(z[m] + H) ** 2 / (2 * sz ** 2)))
    return c


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    X, Y = np.meshgrid(np.linspace(1, 500, 300), np.linspace(-100, 100, 200))
    C = ploma_gaussiana(X, Y, 10.0)
    plt.pcolormesh(X, Y, C, shading="auto")
    plt.colorbar(label="Concentració (u. arb.)")
    plt.xlabel("x (m)"); plt.ylabel("y (m)"); plt.title("Ploma gaussiana (z = H)")
    plt.show()
