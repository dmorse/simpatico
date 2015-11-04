#!/usr/bin/env python

from FileEditor import *

editor = FileEditor()
editor.setIsTest(False)

editor.setFilter(r'^ *Contents: *$')
editor.setOld(r'^ *Contents:')
editor.setNew(r'<h2> Contents:</h2>')
editor.editFiles(".", "*.dox")

editor.setFilter(r'<BR> *<HR> *$')
editor.setOld(r'<BR> *<HR>')
editor.setNew(r'<BR>')
editor.editFiles(".", "*.dox")
