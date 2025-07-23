"""Utility functions for plotting
"""

import math
import matplotlib.pyplot as plt
from matplotlib import colormaps


def blank_plot(
    x_max: float,
    x_min: float,
    y_max: float,
    y_min: float,
    draw_horizontal_lines: bool = False,
):
    """Generate a blank plot with set attributes

    Args:
        x_max (float): maximum x value
        x_min (float): minimum x value
        y_max (float): maximum y value
        y_min (float): minimum y value
        draw_horizontal_lines (bool, optional): whether or not to draw horizontal
        lines across plot. Defaults to False.
    """

    plt.figure(figsize=(7, 5))
    axis = plt.subplot(111)
    axis.spines["top"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)

    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()

    plt.xlim(0.0, x_max)
    plt.ylim(0.0, y_max)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)

    if draw_horizontal_lines:
        inc = (int(y_max) - y_min) / 20
        for i in range(0, 20):
            plt.plot(
                range(math.floor(x_min), math.ceil(x_max)),
                [0.0 + i * inc] * len(range(math.floor(x_min), math.ceil(x_max))),
                "--",
                lw=0.5,
                color="black",
                alpha=0.3,
            )

    plt.tick_params(bottom=False, top=False, left=False, right=False)


def get_color_palette(number: int) -> list:
    """_summary_

    Args:
        number (int): number of colors to get - must be <= 20

    Raises:
        ValueError: number must be less than hard-coded list

    Returns:
        list[tuple]: list of colors to use in plotting
    """

    # hard-coded list of colors, can add more here if necessary
    all_colors = [
        (31, 119, 180),
        (174, 199, 232),
        (255, 127, 14),
        (255, 187, 120),
        (44, 160, 44),
        (152, 223, 138),
        (214, 39, 40),
        (255, 152, 150),
        (148, 103, 189),
        (197, 176, 213),
        (140, 86, 75),
        (196, 156, 148),
        (227, 119, 194),
        (247, 182, 210),
        (127, 127, 127),
        (199, 199, 199),
        (188, 189, 34),
        (219, 219, 141),
        (23, 190, 207),
        (158, 218, 229),
    ]

    if number > len(all_colors):
        raise ValueError(f"get_color_palette: number must be <= {len(all_colors)}")

    colors = [
        (red / 255.0, green / 255.0, blue / 255.0) for red, green, blue in all_colors
    ]

    return colors[:number]


def sample_colormap(this_colormap: str, n_colors: int, i: int) -> tuple:
    """
    Given a colormap and a number of desired colors evenly spaced along it, get the i'th color
    """
    color = colormaps[this_colormap](i / (n_colors - 1))
    return color
