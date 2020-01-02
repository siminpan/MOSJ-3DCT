a = slicer.util.array('spgr')

print(a)

print(a.min(),a.max())

# Defines a python function
# Telling Slicer that the image data for node 'n' has been
# modified causes the slice view windows to refresh
def toggle():
 n = slicer.util.getNode('spgr')
 a = slicer.util.array('spgr')
 a[:] = a.max() - a
 n.Modified()
 print('Toggled')
toggle()

# show a Toggle button
c = qt.QPushButton('Toggle')
c.connect('clicked()',toggle)
b = qt.QPushButton('Toggle2')
b.connect('clicked()',toggle)
b.show()
c.show()
