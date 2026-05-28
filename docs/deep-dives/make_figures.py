"""Generate the intuition figures used by the deep-dive docs in docs/deep-dives/."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import os

OUT = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(0)


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"  saved {path}")


# ---------------------------------------------------------------------------
# Fig 1: step-size sampling via inverse transform  (-log(U))
# ---------------------------------------------------------------------------
def fig_step_sampling():
    U = rng.random(50_000)
    samples = -np.log(U)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    tau = np.linspace(0, 6, 400)
    axes[0].plot(tau, np.exp(-tau), lw=2.2, color="#c0392b")
    axes[0].fill_between(tau, np.exp(-tau), alpha=0.15, color="#c0392b")
    axes[0].set_title(r"PDF of step size:  $p(\tau)=e^{-\tau}$")
    axes[0].set_xlabel(r"$\tau$ (optical depth)")
    axes[0].set_ylabel(r"$p(\tau)$")
    axes[0].grid(alpha=0.3)

    axes[1].hist(samples, bins=80, range=(0, 6), density=True,
                 alpha=0.55, color="#2c7fb8", label=r"$-\ln(U)$ samples")
    axes[1].plot(tau, np.exp(-tau), color="#c0392b", lw=2.2,
                 label=r"$e^{-\tau}$")
    axes[1].set_xlim(0, 6)
    axes[1].set_title("50,000 samples match the exponential")
    axes[1].set_xlabel(r"$\tau$")
    axes[1].set_ylabel("density")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    save(fig, "01_step_sampling.png")


# ---------------------------------------------------------------------------
# Fig 2: WRONG vs RIGHT hemisphere sampling
# ---------------------------------------------------------------------------
def fig_hemisphere_sampling():
    N = 1500

    # WRONG: uniform in theta in [0, pi/2]
    theta_w = rng.uniform(0, np.pi / 2, N)
    phi_w = rng.uniform(0, 2 * np.pi, N)
    x_w = np.sin(theta_w) * np.cos(phi_w)
    y_w = np.sin(theta_w) * np.sin(phi_w)
    z_w = np.cos(theta_w)

    # RIGHT: uniform in cos(theta) in [0, 1]
    cos_r = rng.uniform(0, 1, N)
    sin_r = np.sqrt(1 - cos_r ** 2)
    phi_r = rng.uniform(0, 2 * np.pi, N)
    x_r = sin_r * np.cos(phi_r)
    y_r = sin_r * np.sin(phi_r)
    z_r = cos_r

    fig = plt.figure(figsize=(11, 5))
    for i, (X, Y, Z, title) in enumerate([
        (x_w, y_w, z_w, "WRONG: uniform in θ\n(crowds at the pole)"),
        (x_r, y_r, z_r, "RIGHT: uniform in cos θ\n(spread evenly over solid angle)"),
    ]):
        ax = fig.add_subplot(1, 2, i + 1, projection="3d")
        ax.scatter(X, Y, Z, s=4, alpha=0.55,
                   color="#c0392b" if i == 0 else "#27ae60")
        ax.set_title(title)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z (up)")
        ax.view_init(elev=20, azim=35)
    save(fig, "02_hemisphere_sampling.png")


# ---------------------------------------------------------------------------
# Fig 3: solid angle element on a sphere
# ---------------------------------------------------------------------------
def fig_solid_angle():
    fig = plt.figure(figsize=(7.5, 6))
    ax = fig.add_subplot(111, projection="3d")
    u, v = np.mgrid[0:2 * np.pi:60j, 0:np.pi / 2:30j]
    xs = np.sin(v) * np.cos(u)
    ys = np.sin(v) * np.sin(u)
    zs = np.cos(v)
    ax.plot_wireframe(xs, ys, zs, color="lightgray", linewidth=0.3)

    th = np.linspace(np.pi / 4, np.pi / 4 + 0.25, 12)
    ph = np.linspace(0.3, 0.7, 12)
    TH, PH = np.meshgrid(th, ph)
    xp = np.sin(TH) * np.cos(PH)
    yp = np.sin(TH) * np.sin(PH)
    zp = np.cos(TH)
    ax.plot_surface(xp, yp, zp, color="#f39c12", alpha=0.85)

    ax.set_title("Solid angle element  $dA = \\sin\\theta\\, d\\theta\\, d\\phi$\n"
                 "→ patches near the pole (small sin θ) are tiny\n"
                 "→ sampling uniform in cos θ counters this")
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=25, azim=40)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    save(fig, "03_solid_angle.png")


# ---------------------------------------------------------------------------
# Fig 4: Thomson phase function polar plot
# ---------------------------------------------------------------------------
def fig_thomson_phase():
    th = np.linspace(0, 2 * np.pi, 720)
    mu = np.cos(th)
    p_thom = 0.75 * (1 + mu ** 2)
    p_iso = np.ones_like(th) * 0.5

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 6))
    ax.plot(th, p_thom, color="#2c3e50", lw=2.2, label=r"Thomson  $(3/4)(1+\mu^2)$")
    ax.plot(th, p_iso, color="#95a5a6", lw=1.3, ls="--", label="isotropic (ref.)")
    ax.fill(th, p_thom, alpha=0.15, color="#2c3e50")
    ax.set_theta_zero_location("E")
    ax.set_title("Thomson phase function\n"
                 "scatter angle measured from incoming direction (right = forward)",
                 pad=18)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.15))
    save(fig, "04_thomson_phase.png")


# ---------------------------------------------------------------------------
# Fig 5: rejection sampling visualization
# ---------------------------------------------------------------------------
def fig_rejection():
    mu = np.linspace(-1, 1, 400)
    p = 0.75 * (1 + mu ** 2)

    N = 2000
    mu_prop = rng.uniform(-1, 1, N)
    u = rng.uniform(0, 1.5, N)
    accepted = u < 0.75 * (1 + mu_prop ** 2)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.fill_between(mu, p, alpha=0.18, color="#2c3e50")
    ax.plot(mu, p, color="#2c3e50", lw=2.2, label=r"$p(\mu)=(3/4)(1+\mu^2)$")
    ax.axhline(1.5, color="#c0392b", ls="--", lw=1, label=r"$p_{\max}=3/2$")
    ax.scatter(mu_prop[accepted], u[accepted], s=5, color="#27ae60",
               alpha=0.6, label="accepted")
    ax.scatter(mu_prop[~accepted], u[~accepted], s=5, color="#c0392b",
               alpha=0.4, label="rejected")
    ax.set_xlabel(r"$\mu$ proposal")
    ax.set_ylabel("y proposal")
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1.6)
    ax.set_title("Rejection sampling: keep points UNDER the curve")
    ax.legend(loc="upper center", ncol=2)
    save(fig, "05_rejection_sampling.png")


# ---------------------------------------------------------------------------
# Fig 6: direction rotation geometry
# ---------------------------------------------------------------------------
def fig_rotation_geometry():
    fig = plt.figure(figsize=(8, 6.5))
    ax = fig.add_subplot(111, projection="3d")

    v = np.array([0.4, 0.3, 0.85])
    v /= np.linalg.norm(v)

    theta_s = np.deg2rad(40)
    phis = np.linspace(0, 2 * np.pi, 80)

    helper = np.array([0, 0, 1.0]) if abs(v[2]) < 0.95 else np.array([1, 0, 0.0])
    e1 = np.cross(v, helper); e1 /= np.linalg.norm(e1)
    e2 = np.cross(v, e1)

    cone = (np.cos(theta_s) * v[:, None]
            + np.sin(theta_s) * (np.cos(phis) * e1[:, None]
                                 + np.sin(phis) * e2[:, None]))

    ax.quiver(0, 0, 0, 1, 0, 0, color="lightgray", arrow_length_ratio=0.05)
    ax.quiver(0, 0, 0, 0, 1, 0, color="lightgray", arrow_length_ratio=0.05)
    ax.quiver(0, 0, 0, 0, 0, 1, color="lightgray", arrow_length_ratio=0.05)

    ax.quiver(0, 0, 0, *v, color="#2c3e50", lw=2.5, arrow_length_ratio=0.12,
              label="old direction v")

    ax.plot(cone[0], cone[1], cone[2], color="#f39c12", lw=1.6,
            label=f"all directions at scatter angle θ={int(np.rad2deg(theta_s))}°")

    phi_pick = np.deg2rad(70)
    v_new = (np.cos(theta_s) * v
             + np.sin(theta_s) * (np.cos(phi_pick) * e1
                                  + np.sin(phi_pick) * e2))
    ax.quiver(0, 0, 0, *v_new, color="#c0392b", lw=2.5, arrow_length_ratio=0.12,
              label="new direction (random φ)")

    for k in range(0, 80, 10):
        ax.plot([0, cone[0, k]], [0, cone[1, k]], [0, cone[2, k]],
                color="#f39c12", lw=0.5, alpha=0.4)

    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(0, 1.2)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("Rotating the direction after a scatter\n"
                 "μ=cos θ fixes the cone; φ chooses where on the cone")
    ax.legend(loc="upper left", fontsize=8)
    ax.view_init(elev=18, azim=40)
    save(fig, "06_rotation_geometry.png")


# ---------------------------------------------------------------------------
# Fig 7: random walks at different tau
# ---------------------------------------------------------------------------
def sample_thomson(rng):
    while True:
        m = rng.uniform(-1, 1)
        if rng.uniform(0, 1) < 0.5 * (1 + m ** 2):
            return m


def rotate(v, cos_t, rng):
    phi = rng.uniform(0, 2 * np.pi)
    sin_t = np.sqrt(1 - cos_t ** 2)
    if abs(v[2]) > 0.999:
        return np.array([sin_t * np.cos(phi), sin_t * np.sin(phi),
                         cos_t * np.sign(v[2])])
    s = np.sqrt(1 - v[2] ** 2)
    dx = sin_t * (v[0] * v[2] * np.cos(phi) - v[1] * np.sin(phi)) / s + v[0] * cos_t
    dy = sin_t * (v[1] * v[2] * np.cos(phi) + v[0] * np.sin(phi)) / s + v[1] * cos_t
    dz = -sin_t * s * np.cos(phi) + v[2] * cos_t
    nv = np.array([dx, dy, dz])
    return nv / np.linalg.norm(nv)


def walk(tau_total, rng, max_steps=4000):
    pos = np.array([0.0, 0.0, tau_total])
    cos_th = rng.uniform(0, 1)
    sin_th = np.sqrt(1 - cos_th ** 2)
    phi = rng.uniform(0, 2 * np.pi)
    direction = np.array([sin_th * np.cos(phi), sin_th * np.sin(phi), -cos_th])
    path = [pos.copy()]
    for _ in range(max_steps):
        step = -np.log(rng.uniform(1e-12, 1))
        pos = pos + direction * step
        path.append(pos.copy())
        if pos[2] <= 0 or pos[2] >= tau_total:
            break
        direction = rotate(direction, sample_thomson(rng), rng)
    return np.array(path)


def fig_random_walks():
    fig = plt.figure(figsize=(13, 4.5))
    for i, tau_t in enumerate([0.5, 2.0, 10.0]):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        for _ in range(8):
            p = walk(tau_t, rng)
            escaped = p[-1, 2] <= 0
            ax.plot(p[:, 0], p[:, 1], p[:, 2],
                    color="#27ae60" if escaped else "#c0392b",
                    alpha=0.7, lw=0.9)
        ax.set_title(f"τ_total = {tau_t}\n(green=escaped, red=re-absorbed)")
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("τ depth")
        ax.invert_zaxis()
        ax.set_zlim(tau_t, 0)
    save(fig, "07_random_walks.png")


# ---------------------------------------------------------------------------
# Fig 8: beaming function (Eddington)
# ---------------------------------------------------------------------------
def fig_beaming():
    mu = np.linspace(0, 1, 200)
    I = 1 + 1.5 * mu

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(mu, I, color="#2c3e50", lw=2.5,
            label=r"$I(\mu) \propto 1 + 1.5\,\mu$  (Eddington)")
    ax.fill_between(mu, I, alpha=0.15, color="#2c3e50")
    ax.axhline(I.mean(), color="#95a5a6", ls="--",
               label="flat / isotropic emission (for comparison)")
    ax.set_xlabel(r"$\mu = \cos(\text{angle from surface normal})$")
    ax.set_ylabel("intensity")
    ax.set_title("Beaming function: brighter face-on (μ=1), dimmer at the limb (μ=0)\n"
                 "= limb darkening from scattering atmosphere")
    ax.annotate("photons going straight up\nescape easily",
                xy=(1.0, 2.5), xytext=(0.55, 2.3),
                arrowprops=dict(arrowstyle="->", color="#27ae60"),
                fontsize=9, color="#27ae60")
    ax.annotate("grazing photons traverse\nmore atmosphere → scattered back",
                xy=(0.02, 1.05), xytext=(0.15, 0.4),
                arrowprops=dict(arrowstyle="->", color="#c0392b"),
                fontsize=9, color="#c0392b")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 3)
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    save(fig, "08_beaming_function.png")


if __name__ == "__main__":
    print("Generating figures...")
    fig_step_sampling()
    fig_hemisphere_sampling()
    fig_solid_angle()
    fig_thomson_phase()
    fig_rejection()
    fig_rotation_geometry()
    fig_random_walks()
    fig_beaming()
    print("Done.")
