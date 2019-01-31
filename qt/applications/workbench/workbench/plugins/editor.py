# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#    This file is part of the mantid workbench.
#
#
from __future__ import (absolute_import, unicode_literals)

# system imports
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QVBoxLayout

# third-party library imports
from mantidqt.utils.qt import add_actions, create_action
from mantidqt.widgets.codeeditor.multifileinterpreter import MultiPythonFileInterpreter
# local package imports
from workbench.plugins.base import PluginWidget

# from mantidqt.utils.qt import toQSettings when readSettings/writeSettings are implemented


# Initial content
DEFAULT_CONTENT = """# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *

import matplotlib.pyplot as plt

import numpy as np
"""


class MultiFileEditor(PluginWidget):
    """Provides a tab widget for editing multiple files"""

    def __init__(self, parent):
        super(MultiFileEditor, self).__init__(parent)

        # layout
        self.editors = MultiPythonFileInterpreter(default_content=DEFAULT_CONTENT,
                                                  parent=self)
        layout = QVBoxLayout()
        layout.addWidget(self.editors)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        # attributes
        self.run_action = create_action(
            self, "Run",
            on_triggered=self.editors.execute_current,
            shortcut=("Ctrl+Return", "Ctrl+Enter"),
            shortcut_context=Qt.ApplicationShortcut)

        self.abort_action = create_action(
            self, "Abort", on_triggered=self.editors.abort_current)

        self.toggle_comment_action = create_action(
            self.editors.current_editor(), r"Comment/Uncomment",
            on_triggered=self.editors.toggle_comment_current,
            shortcut="Ctrl+/",
            shortcut_context=Qt.ApplicationShortcut)

        self.tabs_to_spaces_action = create_action(
            self, 'Tabs to Spaces',
            on_triggered=self.editors.tabs_to_spaces_current)

        self.spaces_to_tabs_action = create_action(
            self, 'Spaces to Tabs',
            on_triggered=self.editors.spaces_to_tabs_current)

        self.whitespace_visible_action = create_action(
            self, 'Set Whitespace Visible',
            on_triggered=self.editors.set_whitespace_visible_all)

        self.whitespace_invisible_action = create_action(
            self, 'Set Whitespace Invisible',
            on_triggered=self.editors.set_whitespace_invisible_all)

        self.editor_actions = [self.run_action,
                               self.abort_action,
                               self.toggle_comment_action,
                               self.tabs_to_spaces_action,
                               self.spaces_to_tabs_action,
                               self.whitespace_visible_action,
                               self.whitespace_invisible_action]

    def execute_current(self):
        '''This is used by MainWindow to execute a file after opening it'''
        return self.editors.execute_current()

    # ----------- Plugin API --------------------

    def app_closing(self):
        """
        Tries to close all editors
        :return: True if editors can be closed, false if cancelled
        """
        return self.editors.close_all()

    def get_plugin_title(self):
        return "Editor"

    def readSettings(self, _):
        pass

    def writeSettings(self, _):
        pass

    def register_plugin(self):
        self.main.add_dockwidget(self)
        # menus
        add_actions(self.main.editor_menu, self.editor_actions)

    # ----------- Plugin Behaviour --------------------

    def open_file_in_new_tab(self, filepath):
        return self.editors.open_file_in_new_tab(filepath)

    def save_current_file(self):
        self.editors.save_current_file()
