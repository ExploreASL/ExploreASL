from PySide2.QtWidgets import *
from PySide2.QtCore import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import seaborn as sns


# noinspection PyCallingNonCallable
class xASL_GUI_FacetArtist(QWidget):
    """
    Main Widget class to handle drawing operations around seaborn's FacetGrid class
    """

    def __init__(self, parent):
        super().__init__(parent=parent)
        self.mainlay = QVBoxLayout(self)
        self.parent_cw = parent
        self.manager = parent.fig_manager
        self.generate_canvas(blank=True)

    def generate_canvas(self, blank: bool):
        """
        Generates a canvas and facetgrid, links them together, and adds the result to the artist's main layout
        :param blank: whether the facetgrid should be blank, only initialized with the current data
        """
        if blank:
            self.grid = sns.FacetGrid(data=self.parent_cw.loader.long_data)
            self.mainfig = self.grid.fig
        else:
            # First create the constructor, then feed it into the Facetgrid class
            constructor = {}
            for kwarg, call in self.manager.fig_kwargs.items():
                if call() == "":
                    constructor[kwarg] = None
                else:
                    constructor[kwarg] = call()

            # Then create the FacetGrid from the constructor
            self.grid = sns.FacetGrid(data=self.parent_cw.loader.long_data, **constructor)
            self.mainfig: plt.Figure = self.grid.fig

        # Add this to the widget setup
        self.canvas = FigureCanvas(self.mainfig)
        self.nav = NavigationToolbar(self.canvas, self.canvas)
        self.mainlay.addWidget(self.nav)
        self.mainlay.addWidget(self.canvas)

    def clear_canvas(self):
        """
        Fully clears the canvas and wipes it as well as the facetgrid and toolbar navigator from memory
        """
        plt.clf()
        plt.close(self.mainfig)
        self.mainlay.removeWidget(self.canvas)
        self.mainlay.removeWidget(self.nav)
        self.canvas.setParent(None)
        self.nav.setParent(None)
        del self.nav
        del self.mainfig
        del self.canvas

    def clear_axes(self):
        """
        Convenience call for clearing the axes
        """
        for ax in self.grid.axes.flatten():
            ax.clear()

    #############################################
    # Functions and Slots for updating the canvas
    #############################################

    def plotupdate_axes(self):
        # If x or y axes arguments are still in their non-callable form (i.e string type), don't proceed
        if any([isinstance(self.manager.axes_arg_x, str), isinstance(self.manager.axes_arg_y, str)]):
            return

        # If x and y axes arguments are blank, don't proceed
        if any([self.manager.axes_arg_x() == '', self.manager.axes_arg_y() == '']):
            return

        # Account for a user typing the palette name, accidentally triggering this function, don't proceed
        if any([not self.manager.cmb_palette,
                self.manager.cmb_palette.currentText() not in self.manager.palettenames]):
            return

        print("PLOTUPDATE_AXES TRIGGERED")

        # Clear all axes
        self.clear_axes()

        # Otherwise, proceed
        func = self.manager.plotting_func
        x, y, hue = self.manager.axes_arg_x(), self.manager.axes_arg_y(), self.manager.axes_arg_hue()

        axes_constructor = {}
        for kwarg, call in self.manager.axes_kwargs.items():
            if call() == "":
                axes_constructor[kwarg] = None
            else:
                axes_constructor[kwarg] = call()

        # Account for a user not selecting a palette
        if axes_constructor['palette'] in ['', "None", "Default Blue", "No Palette"]:
            axes_constructor['palette'] = None

        if hue == '':
            self.grid = self.grid.map(func, x, y, data=self.parent_cw.loader.long_data, **axes_constructor)
        else:
            self.grid = self.grid.map(func, x, y, hue, data=self.parent_cw.loader.long_data, **axes_constructor)

        # Redraw the legend
        self.manager.legend_widget.sendSignal_legendparms_updateplot()

        self.plotupdate_padding()

    def plotupdate_padding(self):
        plt.subplots_adjust(**self.manager.padding_kwargs)
        self.canvas.draw()

    @Slot(dict)
    def plotupdate_xlabel(self, xlabel_kwargs: dict):
        for ax in self.mainfig.axes:
            plt.sca(ax)
            if xlabel_kwargs["xlabel"] == '':
                try:
                    x_var = self.manager.le_x.text()
                    cols: list = self.parent_cw.loader.long_data.columns.tolist()
                    xlabel_kwargs["xlabel"] = cols[cols.index(x_var)]
                except AttributeError:
                    pass
            plt.xlabel(**xlabel_kwargs)
        self.canvas.draw()

    @Slot(dict)
    def plotupdate_ylabel(self, ylabel_kwargs: dict):
        for ax in self.mainfig.axes:
            plt.sca(ax)
            if ylabel_kwargs["ylabel"] == '':
                try:
                    y_var = self.manager.le_y.text()
                    cols: list = self.parent_cw.loader.long_data.columns.tolist()
                    ylabel_kwargs["ylabel"] = cols[cols.index(y_var)]
                except AttributeError:
                    pass
            plt.ylabel(**ylabel_kwargs)
        self.canvas.draw()

    @Slot(dict)
    def plotupdate_title(self, title_kwargs: dict):
        if self.manager.le_row.text() == "" and self.manager.le_col.text() == "":
            plt.title(**title_kwargs)
            self.canvas.draw()
        else:
            title_kwargs.pop("loc")
            title_kwargs["t"] = title_kwargs.pop("label")
            title_kwargs["fontsize"] = title_kwargs["fontdict"]["fontsize"]
            title_kwargs["fontweight"] = title_kwargs["fontdict"]["fontweight"]
            del title_kwargs["fontdict"]
            print(title_kwargs)
            plt.suptitle(**title_kwargs)
            self.canvas.draw()

    @Slot()
    def plotupdate_figurecall(self):
        print("PLOTUPDATE_FIGURECALL TRIGGERED")
        self.clear_canvas()
        self.generate_canvas(blank=False)
        self.plotupdate_axes()
        self.plotupdate_padding()

    @Slot()
    def plotupdate_paddingcall(self):
        self.plotupdate_padding()

    @Slot()
    def plotupdate_axescall(self):
        print("PLOTUPDATE_AXESCALL TRIGGERED")
        self.plotupdate_axes()

    @Slot(dict)
    def plotupdate_legendcall(self, legend_kwargs):
        # print("INSIDE PLOTUPDATE_LEGENDCALL")
        if all([len(legend_kwargs) > 0,  # kwargs for the legend function must exist
                self.manager.chk_showlegend.isChecked(),  # the option to show the legend must be there
                self.manager.axes_arg_hue() != ''  # the hue argument cannot be empty
                ]):
            print(legend_kwargs)
            plt.legend(**legend_kwargs)
        else:
            plt.legend("", frameon=False)
        self.canvas.draw()

    @Slot(dict)
    def plotupdate_tickcall(self, ticklabel_kwargs):
        print("INSIDE PLOTDATE_TICKCALL")
        if all([len(ticklabel_kwargs) > 0,  # kwargs for the ticklabels must exist
                ]):
            ax: plt.Axes
            for ax in self.mainfig.axes:
                ax.set_xticklabels(labels=ax.get_xticklabels(), **ticklabel_kwargs)

            self.canvas.draw()