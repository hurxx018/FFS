import sys
import numpy as np
import pyqtgraph as pg
import matplotlib.pyplot as plt


class ACFPlot(object):

    def __init__(self, title=None):
        self.title = title


    def plot(self, x, y):
        if self.title == None:
            p0 = pg.plot(title="Auto-correlation",
                    labels={'left':'G', 'bottom':'t (s)'})
        else:
            p0 = pg.plot(title=self.title,
                    labels={'left':'G', 'bottom':'t (s)'})

        p0.plot(x, y, pen=None, symbol='o', symbolSize=10)
        p0.setLogMode(x=True, y=False)

        if (sys.flags.interactive != 1) or not hasattr(Qtcore, 'PYQT_VERSION'):
            pg.QtGui.QApplication.exec_()
        return


    def mplot(self, x, y, yerr):
        if self.title == None:
            plt.title("Auto-correlation")
            p0 = plt.errorbar(x, y, yerr=yerr, fmt="o", color="blue", ecolor="red")
        else:
            plt.title(self.title)
            p0 = plt.errorbar(x, y, yerr=yerr, fmt="o", color="blue", ecolor="red")

        plt.xlabel("t (s)")
        plt.ylabel("G(t)")
        plt.xscale("log")
        return


    def fitplot(self, x, y, y_fit):
        if self.title == None:
            p0 = pg.plot(title="Auto-correlation with a fit function",
                    labels={'left':'G', 'bottom':'t (s)'})
        else:
            p0 = pg.plot(title=self.title,
                    labels={'left':'G', 'bottom':'t (s)'})

        p0.plot(x, y, pen=None, symbol='o', symbolSize=10)
        p0.plot(x, y_fit, pen=(255, 0, 0))
        p0.setLogMode(x=True, y=False)

        if (sys.flags.interactive != 1) or not hasattr(Qtcore, 'PYQT_VERSION'):
            pg.QtGui.QApplication.exec_()
        return

    def mfitplot(self, x, y, yerr, y_fit):
        if self.title == None:
            plt.title("Auto-correlation with a fit function")
            p0 = plt.errorbar(x, y, yerr=yerr, fmt="o", color="blue", ecolor="red")
        else:
            plt.title(self.title)
            p0 = plt.errorbar(x, y, yerr=yerr, fmt="o", color="blue", ecolor="red")

        plt.plot(x, y_fit, color="green")
        plt.xlabel("t (s)")
        plt.ylabel("G(t)")
        plt.xscale("log")
        return

def main():

    x = np.exp(np.arange(100)*0.1)
    y = np.random.uniform(0, 1, 100)
    yerr = np.random.normal(0, 1, 100)

    p0 = ACFPlot()
    p0.plot(x,y)
    p0.mplot(x, y, yerr)


if __name__ == "__main__":
    main()



    # def plot(self, x, y, yerr):
    #   log scale is not available in ErrorBarItem
    #     p0 = pg.plot(title=self.title)
    #     err0 = pg.ErrorBarItem()
    #     err0.setData(x=x, y=y, height=yerr)
    #     p0.addItem(err0)
    #     p0.plot(x, y, symbol='o')
    #
    #     if (sys.flags.interactive != 1) or not hasattr(Qtcore, 'PYQT_VERSION'):
    #         pg.QtGui.QApplication.exec_()
    #     return
