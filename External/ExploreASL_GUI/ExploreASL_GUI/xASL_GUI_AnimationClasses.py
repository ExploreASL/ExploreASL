from PySide2.QtWidgets import QWidget, QLabel, QVBoxLayout
from PySide2.QtGui import QMovie, Qt
from PySide2.QtCore import QByteArray, QSize


class xASL_ImagePlayer(QWidget):
    def __init__(self, filename, parent=None):
        QWidget.__init__(self, parent)

        # Load the file into a QMovie with the appropriate size
        self.movie = QMovie(filename, QByteArray(), self)
        self.movie.setScaledSize(QSize(75, 75))
        self.movie.setCacheMode(QMovie.CacheAll)
        self.movie.setSpeed(100)

        # Add the QMovie object to a label container
        self.movie_screen = QLabel()
        self.movie_screen.setAlignment(Qt.AlignCenter)
        self.movie_screen.setMovie(self.movie)

        # Create the layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.movie_screen)
        self.setLayout(main_layout)
