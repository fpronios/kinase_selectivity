import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel , QPushButton,QVBoxLayout,QGridLayout
import PyQt5.QtCore

from PyQt5.QtGui import *

import os

mypath = "saved_figs/combined"
onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
image_index=0
onlyfiles.sort()
print(onlyfiles)

class App(QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'PyQt5 image - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        #layout.setColumnStretch(1, 4)
        #layout.setColumnStretch(2, 4)

        print(onlyfiles)
        # Create widget
        self.label = QLabel(self)
        self.pixmap = QPixmap(mypath + "/" + onlyfiles [image_index // len(onlyfiles)])
        self.label.setPixmap(self.pixmap)
        self.resize(self.pixmap.width(), self.pixmap.height())

        self.show()

    def keyPressEvent(self, e):
        global image_index
        if e.key() == PyQt5.QtCore.Qt.Key_Right:

            image_index += 1
            print(image_index)
            #print(onlyfiles[image_index])
            self.pixmap = QPixmap(mypath + "/" + onlyfiles[image_index%len(onlyfiles)])
            self.label.setPixmap(self.pixmap)
            self.resize(self.pixmap.width(), self.pixmap.height())
            self.show()

        if e.key() == PyQt5.QtCore.Qt.Key_Left:

            image_index -= 1
            print(image_index)
           # print(onlyfiles[image_index])
            self.pixmap = QPixmap(mypath + "/" + onlyfiles[image_index%len(onlyfiles)])
            self.label.setPixmap(self.pixmap)
            self.resize(self.pixmap.width(), self.pixmap.height())
            self.show()



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())