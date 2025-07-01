# trace generated using paraview version 5.12.0-RC2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
lUCshiftcut1alltiff = FindSource('L-UCshift-cut1-all.tiff')

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=lUCshiftcut1alltiff)

# Properties modified on calculator1 wyf
calculator1.Function = '"Tiff Scalars"/8'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Outline'

# hide data in view
Hide(lUCshiftcut1alltiff, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'TiffScalars'
tiffScalarsLUT = GetColorTransferFunction('TiffScalars')

# get opacity transfer function/opacity map for 'TiffScalars'
tiffScalarsPWF = GetOpacityTransferFunction('TiffScalars')

# get 2D transfer function for 'TiffScalars'
tiffScalarsTF2D = GetTransferFunction2D('TiffScalars')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1443, 704)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [456.6034240722656, 712.4171752929688, 4171.61313256309]
renderView1.CameraFocalPoint = [456.6034240722656, 712.4171752929688, 456.6034240722656]
renderView1.CameraParallelScale = 961.5152652981891


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------