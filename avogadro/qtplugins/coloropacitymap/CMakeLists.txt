set(_srcs
  coloropacitymap.cpp
  histogramwidget.cpp
  qvtkwidget.cpp
  vtkChartHistogram.cxx
  vtkChartHistogramColorOpacityEditor.cxx
  vtkCustomPiecewiseControlPointsItem.cxx
#  pdfoptionsdialog.cpp
)

set(_uis
#  pdfoptionsdialog.ui
)

avogadro_plugin(ColorOpacityMap
  "Edit the color opacity map."
  ExtensionPlugin
  coloropacitymap.h
  ColorOpacityMap
  "${_srcs}"
  "${_uis}"
)

target_link_libraries(ColorOpacityMap LINK_PRIVATE AvogadroVtk)
