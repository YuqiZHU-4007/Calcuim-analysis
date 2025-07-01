# trace generated using paraview version 5.12.0-RC2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Feature Point Data Reader'
cSUSshiftcuttrail_Cond1_Alltxt = FeaturePointDataReader(registrationName='CSoff-avtivation_Cond.1_All.txt', FileName='W:/ind_all_emergeCSUS_inthissession/Learner/CSoff-avtivation/CSoff-avtivation_Cond.1_All.txt')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
cSUSshiftcuttrail_Cond1_AlltxtDisplay = Show(cSUSshiftcuttrail_Cond1_Alltxt, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
cSUSshiftcuttrail_Cond1_AlltxtDisplay.Representation = 'Surface'

# show color bar/color legend
cSUSshiftcuttrail_Cond1_AlltxtDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'FeaturePointIndex'
featurePointIndexLUT = GetColorTransferFunction('FeaturePointIndex')

# get opacity transfer function/opacity map for 'FeaturePointIndex'
featurePointIndexPWF = GetOpacityTransferFunction('FeaturePointIndex')

# get 2D transfer function for 'FeaturePointIndex'
featurePointIndexTF2D = GetTransferFunction2D('FeaturePointIndex')

# create a new 'Feature Point Density'
featurePointDensity1 = FeaturePointDensity(registrationName='FeaturePointDensity1', Input=cSUSshiftcuttrail_Cond1_Alltxt)

# Properties modified on featurePointDensity1 wyf
featurePointDensity1.SampleDimensions = [200, 200, 200]

# show data in view
featurePointDensity1Display = Show(featurePointDensity1, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
featurePointDensity1Display.Representation = 'Volume'

# hide data in view
Hide(cSUSshiftcuttrail_Cond1_Alltxt, renderView1)

# show color bar/color legend
featurePointDensity1Display.SetScalarBarVisibility(renderView1, True)

# show data in view
featurePointDensity1Display_1 = Show(OutputPort(featurePointDensity1, 1), renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
featurePointDensity1Display_1.Representation = 'Surface'

# hide data in view
Hide(cSUSshiftcuttrail_Cond1_Alltxt, renderView1)

# show color bar/color legend
featurePointDensity1Display_1.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'DensityScalars'
densityScalarsLUT = GetColorTransferFunction('DensityScalars')

# get opacity transfer function/opacity map for 'DensityScalars'
densityScalarsPWF = GetOpacityTransferFunction('DensityScalars')

# get 2D transfer function for 'DensityScalars'
densityScalarsTF2D = GetTransferFunction2D('DensityScalars')

# hide data in view
Hide(OutputPort(featurePointDensity1, 1), renderView1)

# set scalar coloring using an separate color/opacity maps
ColorBy(featurePointDensity1Display, ('POINTS', 'DensityScalars'), True)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(densityScalarsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
featurePointDensity1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
featurePointDensity1Display.SetScalarBarVisibility(renderView1, True)

# get separate color transfer function/color map for 'DensityScalars'
separate_featurePointDensity1Display_DensityScalarsLUT = GetColorTransferFunction('DensityScalars', featurePointDensity1Display, separate=True)

# get separate opacity transfer function/opacity map for 'DensityScalars'
separate_featurePointDensity1Display_DensityScalarsPWF = GetOpacityTransferFunction('DensityScalars', featurePointDensity1Display, separate=True)

# get separate 2D transfer function for 'DensityScalars'
separate_featurePointDensity1Display_DensityScalarsTF2D = GetTransferFunction2D('DensityScalars', featurePointDensity1Display, separate=True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
separate_featurePointDensity1Display_DensityScalarsLUT.ApplyPreset('Black-Body Radiation', True)

# Rescale transfer function wyf
separate_featurePointDensity1Display_DensityScalarsLUT.RescaleTransferFunction(0.0, 10.0)

# Rescale transfer function
separate_featurePointDensity1Display_DensityScalarsPWF.RescaleTransferFunction(0.0, 10.0)

# Rescale 2D transfer function
separate_featurePointDensity1Display_DensityScalarsTF2D.RescaleTransferFunction(0.0, 10.0, 0.0, 1.0)

# set active source
SetActiveSource(featurePointDensity1)

# save data
SaveData('D:\\data\\learning\\L-UCshift-cut1-all.tiff', proxy=featurePointDensity1, Scalars=['POINTS', 'DensityScalars'])

# get layout
layout1 = GetLayout()

#Exit preview mode
layout1.PreviewMode = [0, 0]

# get color legend/bar for separate_featurePointDensity1Display_DensityScalarsLUT in view renderView1
separate_featurePointDensity1Display_DensityScalarsLUTColorBar = GetScalarBar(separate_featurePointDensity1Display_DensityScalarsLUT, renderView1)

# change scalar bar placement
separate_featurePointDensity1Display_DensityScalarsLUTColorBar.Position = [0.8123839223839224, 0.10522727272727275]
separate_featurePointDensity1Display_DensityScalarsLUTColorBar.ScalarBarLength = 0.13000000000000078

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# layout/tab size in pixels
layout1.SetSize(1443, 704)

# current camera placement for renderView1 wyf
renderView1.CameraPosition = [659.35014363921, 777.9717029882187, -1584.357182664869]
renderView1.CameraFocalPoint = [659.35014363921, 777.9717029882187, 198.0745086669922]
renderView1.CameraViewUp = [1.0, 0.0, 0.0]
renderView1.CameraParallelScale = 988.8959700992922

# save screenshot
SaveScreenshot('D:\data\learning\L-UCshift-1all.png', renderView1, 16, ImageResolution=[3840, 2160],
    TransparentBackground=1)

# Properties modified on renderView1 wyf
renderView1.CenterOfRotation = [655.677131652832, 705.0399551391602, 246.90121841430664]
renderView1.CameraPosition = [-2127.7427411820777, 705.0399551391602, 246.90121841430664]
renderView1.CameraFocalPoint = [655.677131652832, 705.0399551391602, 246.90121841430664]
renderView1.CameraViewUp = [0.0, 0.0, -1.0]
renderView1.CameraParallelScale = 720.4020736065111

# layout/tab size in pixels
layout1.SetSize(1443, 704)

# current camera placement for renderView1
renderView1.CameraPosition = [-915.4908229560742, 705.0399551391602, 246.90121841430664]
renderView1.CameraFocalPoint = [655.677131652832, 705.0399551391602, 246.90121841430664]
renderView1.CameraViewUp = [0.0, 0.0, -1.0]
renderView1.CameraParallelScale = 720.4020736065111

# save screenshot
SaveScreenshot('D:\data\learning\L-UCshift-1all-2.png', renderView1, 16, ImageResolution=[3840, 2160],
    TransparentBackground=1)

# Properties modified on renderView1
renderView1.CameraPosition = [655.677131652832, -866.1279994697461, 246.90121841430664]

# layout/tab size in pixels
layout1.SetSize(1443, 704)

# current camera placement for renderView1
renderView1.CameraPosition = [655.677131652832, -866.1279994697461, 246.90121841430664]
renderView1.CameraFocalPoint = [655.677131652832, 705.0399551391602, 246.90121841430664]
renderView1.CameraViewUp = [0.0, 0.0, -1.0]
renderView1.CameraParallelScale = 720.4020736065111

# save screenshot
SaveScreenshot('D:\data\learning\L-UCshift-1all-3.png', renderView1, 16, ImageResolution=[3840, 2160],
    TransparentBackground=1)

# set active source
SetActiveSource(featurePointDensity1)

# get color transfer function/color map for 'neuronUID'
neuronUIDLUT = GetColorTransferFunction('neuronUID')

# get opacity transfer function/opacity map for 'neuronUID'
neuronUIDPWF = GetOpacityTransferFunction('neuronUID')

# get 2D transfer function for 'neuronUID'
neuronUIDTF2D = GetTransferFunction2D('neuronUID')

# set active source
SetActiveSource(cSUSshiftcuttrail_Cond1_Alltxt)

# hide data in view
Hide(featurePointDensity1, renderView1)

# show data in view
cSUSshiftcuttrail_Cond1_AlltxtDisplay = Show(cSUSshiftcuttrail_Cond1_Alltxt, renderView1, 'GeometryRepresentation')

# show color bar/color legend
cSUSshiftcuttrail_Cond1_AlltxtDisplay.SetScalarBarVisibility(renderView1, True)

# destroy featurePointDensity1
Delete(featurePointDensity1)
del featurePointDensity1

# destroy cSUSshiftcuttrail_Cond1_Alltxt
Delete(cSUSshiftcuttrail_Cond1_Alltxt)
del cSUSshiftcuttrail_Cond1_Alltxt

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1443, 704)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [655.677131652832, -866.1279994697461, 246.90121841430664]
renderView1.CameraFocalPoint = [655.677131652832, 705.0399551391602, 246.90121841430664]
renderView1.CameraViewUp = [0.0, 0.0, -1.0]
renderView1.CameraParallelScale = 720.4020736065111


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