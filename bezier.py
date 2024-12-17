import sys
import numpy as np
from PyQt5.QtWidgets import QApplication
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from scipy.special import comb


class BezierCurveInteractive(FigureCanvas):

    def __init__(self, n_points=4):
        self.fig, self.ax = plt.subplots()
        super().__init__(self.fig)
        self.grid_size = 10
        self.n_points = n_points
        self.control_points = self.initialize_grid_points(self.n_points)
        self.dragging_point = None

        self.control_circles = [
            Circle(p, 0.3, color="red", picker=True) for p in self.control_points
        ]
        for circle in self.control_circles:
            self.ax.add_patch(circle)

        self.control_line = Line2D(
            self.control_points[:, 0],
            self.control_points[:, 1],
            linestyle="--",
            color="gray",
        )
        self.ax.add_line(self.control_line)

        self.bezier_curve, = self.ax.plot([], [], "b-", lw=2)
        self.update_curve()

        self.fig.canvas.mpl_connect("pick_event", self.on_pick)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)

        self.ax.set_xlim(0, 10)
        self.ax.set_ylim(0, 10)
        self.ax.set_aspect("equal", adjustable="datalim")
        self.draw()

    def initialize_grid_points(self, n_points):
        grid_size = int(np.ceil(np.sqrt(n_points)))  
        points = []
        step = self.grid_size / (grid_size + 1)  
        for i in range(n_points):
            row = i // grid_size
            col = i % grid_size
            x = (col + 1) * step
            y = (row + 1) * step
            points.append([x, y])
        return np.array(points)


    def bernstein_poly(self, i, n, t):
        return comb(n, i) * (t**i) * (1 - t) ** (n - i)

    def compute_bezier_curve(self, points, num=100):
        n = len(points) - 1
        t = np.linspace(0, 1, num)
        curve = np.zeros((num, 2))
        for i in range(n + 1):
            curve += np.outer(self.bernstein_poly(i, n, t), points[i])
        return curve

    def update_curve(self):
        curve = self.compute_bezier_curve(self.control_points)
        self.bezier_curve.set_data(curve[:, 0], curve[:, 1])

        self.control_line.set_data(
            self.control_points[:, 0], self.control_points[:, 1]
        )
        self.draw()

    def on_pick(self, event):
        if isinstance(event.artist, Circle):
            self.dragging_point = event.artist

    def on_motion(self, event):
        if self.dragging_point is not None and event.xdata is not None and event.ydata is not None:
            self.dragging_point.center = (event.xdata, event.ydata)
            idx = self.control_circles.index(self.dragging_point)
            self.control_points[idx] = [event.xdata, event.ydata]
            self.update_curve()

    def on_release(self, event):
        self.dragging_point = None


def main(n_points=4):
    app = QApplication(sys.argv)
    main_window = BezierCurveInteractive(n_points)
    main_window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    POINTS = 2
    if len(sys.argv) > 1:
        try:
            n_points = int(sys.argv[1])
        except ValueError:
            print("Invalid input. Please provide an integer for the number of points.")
            sys.exit(1)
    else:
        n_points = POINTS  

    main(n_points)
