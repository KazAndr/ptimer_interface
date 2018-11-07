# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 18:57:43 2018

@author: Kazantsev Andrey
"""

"""


cmd = "ptimer.exe 1133.tpo -i -m -d ptimer.cfg"

os.system(cmd)
"""

import os
import sys
import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5.QtWidgets import (QMainWindow, QApplication, QLabel, QRadioButton,
                             QVBoxLayout, QFrame, QWidget, QGridLayout,
                             QCheckBox, QPushButton, QFileDialog, QSizePolicy,
                             QDialog, QMessageBox, QDesktopWidget)
from PyQt5.QtCore import Qt
from numpy import genfromtxt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class ShowParam(QWidget):

    def __init__(self, title, file, pos, parent=None):
        super(ShowParam, self).__init__()

        self.param_label = QLabel(self)
        self.title = title
        self.file = file
        self.left = 20 + 300
        self.top = 30 + pos
        self.width = 300
        self.height = 370
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        with open(self.file, 'r') as f:
            text = f.readlines()

        self.param_label.setText(''.join(text))
        self.param_label.setSizePolicy(
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        self.param_label.setAlignment(Qt.AlignLeft)
        workstate = QVBoxLayout()
        workstate.addWidget(self.param_label, 0, Qt.AlignTop | Qt.AlignLeft)
        self.setLayout(workstate)


class WarningBox(QMessageBox):

    def __init__(self, title, text):
        super().__init__()
        self.title = title
        self.text = text
        self.initUI()

    def initUI(self):
        self.setIcon(QMessageBox.Warning)
        self.setWindowTitle('Предупреждение!')
        self.setText(self.title)
        self.setInformativeText(self.text)
        self.setStandardButtons(
                QMessageBox.Ok)
        self.setDefaultButton(QMessageBox.Ok)


class drawPlot(QWidget):

    def __init__(
            self, data, up_title, OX_lable, OY_lable, idx_geom):
        super().__init__()
        self.title = up_title
        self.idx_geom = idx_geom
        self.left, self.top, self.width, self.height = self.get_geom()
        self.data = data
        self.up_title = up_title
        self.OX_lable = OX_lable
        self.OY_lable = OY_lable
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        m = PlotCanvas(
                self.data, self.up_title, self.OX_lable, self.OY_lable)

        hbox = QGridLayout()
        # addWidget(QWidget, row, column, rows, columns)
        hbox.addWidget(m, 0, 0)
        self.setLayout(hbox)

    def get_geom(self):
        resolution = QDesktopWidget().screenGeometry()
        height = int(resolution.height()/3)
        width = int(height*1.5)
        if self.idx_geom == 1:
            left = resolution.width() - width
            top = 20
        elif self.idx_geom == 2:
            left = resolution.width() - 2*width
            top = 20

        return left, top, width, height


class PlotCanvas(FigureCanvas):

    def __init__(self, data, up_title, OX_lable, OY_lable):

        width = 4
        height = 3
        dpi = 100
        self.data = data
        self.up_title = up_title
        self.OX_lable = OX_lable
        self.OY_lable = OY_lable
        fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, fig)
        self.setParent(None)

        FigureCanvas.setSizePolicy(
                self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding
                )
        FigureCanvas.updateGeometry(self)
        self.plot()

    def plot(self):
        ax = self.figure.add_subplot(111)
        ax.errorbar(self.data[0], self.data[1], yerr=self.data[2],
                    fmt='o', markersize=1.4, ecolor='g', elinewidth=0.5,
                    capsize=3)
        ax.set_title(self.up_title, fontsize=10)
        ax.set_xlabel(self.OX_lable, fontsize=10)
        ax.set_ylabel(self.OY_lable, fontsize=10)
        ax.grid()
        self.draw()


class QuestionBox(QDialog):

    def __init__(self, title, text):
        super().__init__()
        self.title = title
        self.text = text
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Вопрос!')
        self.label = QLabel(self.text)
        self.plot_resid = QCheckBox('ОУ МПИ')
        self.plot_pcsh = QCheckBox('Ряд ПВШ')
        self.yesBut = QPushButton('Да')
        self.noBut = QPushButton('Нет')

        '''
        Подключение функций к кнопкам
        '''
        self.noBut.clicked.connect(self.close)
        self.yesBut.clicked.connect(self.ploting)

        hbox = QGridLayout()
        # addWidget(QWidget, row, column, rows, columns)
        hbox.addWidget(self.label, 0, 0, 1, 2)
        hbox.addWidget(self.plot_resid, 2, 0)
        hbox.addWidget(self.plot_pcsh, 2, 1)
        hbox.addWidget(self.yesBut, 3, 0)
        hbox.addWidget(self.noBut, 3, 2)
        self.setLayout(hbox)

    def ploting(self):
        self.close()
        if not self.plot_resid.isChecked() and not self.plot_pcsh.isChecked():
            self.msg = QuestionBox(
                    'Данные для визуализации готовы',
                    'Необходимо выбрать хотя бы один вариант данных для отрисовки'
                    )
            self.msg.exec_()
        else:
            if self.plot_resid.isChecked():
                data = genfromtxt('_res.out').T
                self.draw_resid = drawPlot(
                            data, 'Остаточные уклонения',
                            'Дата, MJD', 'ОУ МПИ, мкс', 1)
                self.draw_resid.show()
            if self.plot_pcsh.isChecked():
                data = genfromtxt('_scl.res').T
                self.draw_pcsh = drawPlot(
                            data,
                            'Ряд ОУ МПИ для пульсарной шкалы',
                            'Дни от начальной даты',
                            'Отклонения от опорной шкалы, секунд', 2)
                self.draw_pcsh.show()


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Ptimer'
        self.left = 20
        self.top = 30
        self.width = 300
        self.height = 370
        self.out_dir = os.getcwd() + os.sep
        self.initUI()

    def initUI(self):
        '''
        Применение основных параметров окна
        '''
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        '''
        Инициализация первой рабочей области с выбором режима работы
        '''
        self.label_workstate = QLabel('Выбор режима работы')
        self.label_workstate.setObjectName('label_workstate')
        self.label_workstate.setStyleSheet(
            'QLabel#label_workstate {font: bold;}')
        self.workstate_chron = QRadioButton('Хронометрирование')
        self.workstate_toa = QRadioButton('Генерация МПИ')

        '''
        Инициализация второй рабочей области с выбором опций
        '''
        self.label_option = QLabel('Выбор опций')
        self.label_option.setObjectName('label_option')
        self.label_option.setStyleSheet(
            'QLabel#label_option {font: bold;}')
        self.opt_param = QCheckBox('Не показывать начальные параметры')
        self.opt_aver = QCheckBox('Не вычитать среднее значение')
        self.opt_inp = QCheckBox(
                'Вывод уточненных параметров в файл _tim.inp')
        self.opt_apr = QCheckBox(
                'Вывод точной корреляционной матрицы в файл _tim.apr')
        self.opt_cfg = QCheckBox('Выбрать концигурационный файл')
        '''
        Компановка объектов первой рабочей области
        '''
        layout_workstate = QVBoxLayout()
        # layout.setContentsMargins(5, 5, 5, 5)
        layout_workstate.addWidget(
                self.label_workstate,
                0,
                Qt.AlignTop | Qt.AlignCenter
                )
        layout_workstate.addWidget(
                self.workstate_chron,
                0,
                Qt.AlignTop | Qt.AlignLeft
                )
        layout_workstate.addWidget(
                self.workstate_toa,
                0,
                Qt.AlignTop | Qt.AlignLeft)

        '''
        Компановка объектов второй рабочей области
        '''
        layout_option = QVBoxLayout()
        # layout.setContentsMargins(5, 5, 5, 5)
        layout_option.addWidget(
                self.label_option,
                0,
                Qt.AlignTop | Qt.AlignCenter
                )
        layout_option.addWidget(
                self.opt_param,
                0,
                Qt.AlignTop | Qt.AlignLeft
                )
        layout_option.addWidget(
                self.opt_aver,
                0,
                Qt.AlignTop | Qt.AlignLeft
                )
        layout_option.addWidget(
                self.opt_inp,
                0,
                Qt.AlignTop | Qt.AlignLeft
                )
        layout_option.addWidget(
                self.opt_apr,
                0,
                Qt.AlignTop | Qt.AlignLeft
                )
        layout_option.addWidget(
                self.opt_cfg,
                0,
                Qt.AlignTop | Qt.AlignLeft
                )

        '''
        Кнопка, запускающая обработку
        '''
        self.procces_button = QPushButton('Начать вычисление')
        self.procces_button.setObjectName('self.procces_button')
        self.procces_button.setStyleSheet(
            'QPushButton#self.procces_button {font: bold; background-color: gray;}')

        '''
        Установка значений по умолчанию
        '''
        self.workstate_chron.setChecked(True)
        self.statusBar().showMessage('Готов')
        '''
        Подключение функций к кнопкам
        '''
        self.workstate_chron.toggled.connect(lambda: self.check_active())
        self.procces_button.clicked.connect(self.processing)
        '''
        Инициализация фрейма для первой рабочей области
        Упаковка рабочей области
        '''
        frame_workstate = QFrame()
        frame_workstate.setObjectName('frame_workstate')
        frame_workstate.setLayout(layout_workstate)
        frame_workstate.setStyleSheet(
                'QFrame#frame_workstate {border:1px solid black;}')

        '''
        Инициализация фрейма для первой рабочей области
        Упаковка рабочей области
        '''
        frame_option = QFrame()
        frame_option.setObjectName('frame_option')
        frame_option.setLayout(layout_option)
        frame_option.setStyleSheet(
                'QFrame#frame_option {border:1px solid black;}')

        '''
        Упаковка фреймов
        '''
        hbox = QGridLayout()
        # addWidget(QWidget, row, column, rows, columns)
        hbox.addWidget(frame_workstate, 0, 0, 2, 2)
        hbox.addWidget(frame_option, 3, 0, 3, 2)
        hbox.addWidget(self.procces_button, 7, 1, 1, 1)

        main_panel = QWidget()
        main_panel.setLayout(hbox)
        main_panel.setMaximumSize(self.width, self.height)

        self.setCentralWidget(main_panel)
        self.check_active()
        self.show()

        '''
        Инициализация рабочих функций
        '''
    def check_active(self):
        if self.workstate_chron.isChecked():
            self.opt_param.setEnabled(True)
            self.opt_aver.setEnabled(True)
            self.opt_inp.setEnabled(True)
            self.opt_apr.setEnabled(True)
        else:
            self.opt_param.setEnabled(False)
            self.opt_aver.setEnabled(False)
            self.opt_inp.setEnabled(False)
            self.opt_apr.setEnabled(False)

    def processing(self):
        cmd = 'ptimer.exe '

        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog

        self.filePar, _ = QFileDialog.getOpenFileName(
                    self, 'Выберите файл параметров',
                    '', 'Все файлы (*)',
                    options=options)

        if self.workstate_chron.isChecked():
            if self.opt_param.isChecked():
                cmd += '-c '
            if self.opt_aver.isChecked():
                cmd += '-r '
            if self.opt_inp.isChecked():
                cmd += '-i '
            if self.opt_apr.isChecked():
                cmd += '-m '
            if self.opt_cfg.isChecked():
                self.fileConfig, _ = QFileDialog.getOpenFileName(
                        self, 'Выберите конфигурационный файл',
                        '', 'Файл конфигурации (*.cfg);; Все файлы (*)',
                        options=options)
                cmd += '-d' + self.fileConfig + ' '
            cmd += self.filePar

            self.statusBar().showMessage(
                    "Выполняется команда " + cmd)
            os.system(cmd)

            if not self.opt_param.isChecked():
                self.dialog_old = ShowParam(
                        'Параметры по умолчанию', self.filePar, 0)
                self.dialog_old.show()

            if self.opt_inp.isChecked():
                self.dialog_new = ShowParam(
                        'Уточненные параметры', '_tim.inp', 370)
                self.dialog_new.show()

            self.msg = QuestionBox(
                    'Данные для визуализации готовы',
                    'Провести отрисовку графиков?'
                    )
            self.msg.exec_()

        else:
            if self.opt_cfg.isChecked():
                self.fileConfig, _ = QFileDialog.getOpenFileName(
                        self, 'Выберите конфигурационный файл',
                        '', 'Файл конфигурации (*.cfg);; Все файлы (*)',
                        options=options)
                self.fileConfig = '-d' + self.fileConfig
            else:
                self.fileConfig = ''

            cmd += '-f '
            cmd += self.fileConfig + ' '
            cmd += self.filePar
            self.statusBar().showMessage(
                    "Выполняется команда " + cmd)
            os.system(cmd)

        msg = WarningBox(
                    'Данные успешно обработаны',
                    '')
        msg.exec_()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
