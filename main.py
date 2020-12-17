#!/bin/env python


if __name__ == '__main__':

    #Qt stuff:
    from GisansExplorer.main_app import App
    from GisansExplorer.utils import enable_high_dpi_scaling, profile_function_with_arguments, __DEBUG__
    import PyQt5.QtWidgets as qtw
    import sys

    enable_high_dpi_scaling()
    app = qtw.QApplication(sys.argv)
    app.setStyleSheet("QMessageBox { messagebox-text-interaction-flags: 5;}")

    ex = App()
    #sys.exit(app.exec_())
    if __DEBUG__:
        x = profile_function_with_arguments(app.exec_)
    else:
        x = app.exec()
    sys.exit(x)