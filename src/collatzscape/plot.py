from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from .simulate import iterate

def basin_plot(*, re_min: float, re_max: float, im_min: float, im_max: float,
               n_re: int, n_im: int,
               c: complex, alpha: float, gamma: float, direction: int, steps: int,
               kick_mode: str = "tanh", C_mode: str = "cosine", paired: bool = True,
               escape_radius: float = 1e6, outpath: str | None = None):
    """
    Escape-set plot: pixel=escaped?
    """
    xs = np.linspace(re_min, re_max, n_re)
    ys = np.linspace(im_min, im_max, n_im)
    img = np.zeros((n_im, n_re), dtype=np.uint8)

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            z0 = complex(float(x), float(y))
            _, escaped, _ = iterate(z0, c=c, alpha=alpha, gamma=gamma, direction=direction, steps=steps,
                                    kick_mode=kick_mode, C_mode=C_mode, paired=paired,
                                    escape_radius=escape_radius)
            img[j, i] = 255 if escaped else 0

    plt.figure()
    plt.imshow(img, extent=[re_min, re_max, im_min, im_max], origin="lower")
    plt.xlabel("Re(z0)")
    plt.ylabel("Im(z0)")
    plt.title(f"Escape set (alpha={alpha}, gamma={gamma}, dir={direction})")
    if outpath:
        plt.savefig(outpath, dpi=200, bbox_inches="tight")
    else:
        plt.show()
    plt.close()

def basin_labels_plot(labels, *, re_min, re_max, im_min, im_max, outpath=None, title="Basin labels (proxy)"):
    plt.figure()
    plt.imshow(labels, extent=[re_min, re_max, im_min, im_max], origin="lower", interpolation="nearest")
    plt.xlabel("Re(z0)")
    plt.ylabel("Im(z0)")
    plt.title(title)
    if outpath:
        plt.savefig(outpath, dpi=200, bbox_inches="tight")
    else:
        plt.show()
    plt.close()


def fatou_components_plot(labels, *, re_min, re_max, im_min, im_max, outpath=None,
                          title="Fatou components (proxy)"):
    """
    Plot integer basin labels from fatou_label_grid:
      -1 = escaped
       0..K-1 = attractor basin ids.
    """
    basin_labels_plot(
        labels,
        re_min=re_min,
        re_max=re_max,
        im_min=im_min,
        im_max=im_max,
        outpath=outpath,
        title=title,
    )
