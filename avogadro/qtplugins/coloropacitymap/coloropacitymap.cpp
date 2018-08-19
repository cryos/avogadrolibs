/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2018 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "coloropacitymap.h"
#include "computehistogram.h"
#include "histogramwidget.h"

#include <QAction>
#include <QDialog>
#include <QMessageBox>
#include <QString>

#include <avogadro/core/crystaltools.h>
#include <avogadro/core/unitcell.h>
#include <avogadro/core/cube.h>
#include <avogadro/qtgui/molecule.h>
#include <avogadro/qtopengl/activeobjects.h>
#include <avogadro/vtk/vtkplot.h>
#include <avogadro/vtk/vtkglwidget.h>

#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkTable.h>
#include <vtkRenderWindow.h>

#include <QDebug>

using Avogadro::Core::CrystalTools;
using Avogadro::Core::UnitCell;
using Avogadro::QtGui::Molecule;
using Avogadro::QtOpenGL::ActiveObjects;

using std::map;

namespace Avogadro {
namespace QtPlugins {

using Core::Array;

vtkImageData* cubeImageData(Core::Cube* cube)
{
  auto data = vtkImageData::New();
  // data->SetNumberOfScalarComponents(1, nullptr);
  Eigen::Vector3i dim = cube->dimensions();
  data->SetExtent(0, dim.x() - 1, 0, dim.y() - 1, 0, dim.z() - 1);

  data->SetOrigin(cube->min().x(), cube->min().y(), cube->min().z());
  data->SetSpacing(cube->spacing().data());

  data->AllocateScalars(VTK_DOUBLE, 1);

  double* dataPtr = static_cast<double*>(data->GetScalarPointer());
  std::vector<double>* cubePtr = cube->data();

  for (int i = 0; i < dim.x(); ++i)
    for (int j = 0; j < dim.y(); ++j)
      for (int k = 0; k < dim.z(); ++k) {
        dataPtr[(k * dim.y() + j) * dim.x() + i] =
          (*cubePtr)[(i * dim.y() + j) * dim.z() + k];
      }

  return data;
}

ColorOpacityMap::ColorOpacityMap(QObject* p)
  : Avogadro::QtGui::ExtensionPlugin(p)
  , m_actions(QList<QAction*>())
  , m_molecule(nullptr)
  , m_histogramWidget(nullptr)
  , m_displayDialogAction(new QAction(this))
{
  m_displayDialogAction->setText(tr("Edit Color Opacity Map..."));
  connect(m_displayDialogAction.data(), &QAction::triggered, this,
          &ColorOpacityMap::displayDialog);
  m_actions.push_back(m_displayDialogAction.data());
  m_displayDialogAction->setProperty("menu priority", 70);

  updateActions();
}

ColorOpacityMap::~ColorOpacityMap() = default;

QList<QAction*> ColorOpacityMap::actions() const
{
  return m_actions;
}

QStringList ColorOpacityMap::menuPath(QAction*) const
{
  return QStringList() << tr("&Extensions");
}

void ColorOpacityMap::setMolecule(QtGui::Molecule* mol)
{
  if (m_molecule == mol)
    return;

  if (m_molecule)
    m_molecule->disconnect(this);

  m_molecule = mol;

  if (m_molecule)
    connect(m_molecule, SIGNAL(changed(uint)), SLOT(moleculeChanged(uint)));

  updateActions();
}

void ColorOpacityMap::moleculeChanged(unsigned int c)
{
  Q_ASSERT(m_molecule == qobject_cast<Molecule*>(sender()));

  Molecule::MoleculeChanges changes = static_cast<Molecule::MoleculeChanges>(c);

  if (changes & Molecule::UnitCell) {
    if (changes & Molecule::Added || changes & Molecule::Removed)
      updateActions();
  }
}

void ColorOpacityMap::updateActions()
{
  // Disable everything for nullptr molecules.
  if (!m_molecule) {
    foreach (QAction* action, m_actions)
      action->setEnabled(false);
    return;
  }
  foreach (QAction* action, m_actions)
    action->setEnabled(true);
}

void ColorOpacityMap::displayDialog()
{
  if (!m_histogramWidget) {
    m_histogramWidget = new HistogramWidget;
    m_histogramWidget->resize(800, 600);
    connect(m_histogramWidget, SIGNAL(colorMapUpdated()), SLOT(render()));
    connect(m_histogramWidget, SIGNAL(opacityChanged()), SLOT(render()));
  }

  auto widget = ActiveObjects::instance().activeWidget();
  qDebug() << "Found" << widget;
  auto vtkWidget = qobject_cast<VTK::vtkGLWidget*>(widget);
  qDebug() << "Found" << vtkWidget;

  if (m_molecule && m_molecule->cubeCount()) {
    vtkNew<vtkTable> table;
    auto imageData = vtkWidget->imageData();
    auto lut = vtkWidget->lut();
    auto opacity = vtkWidget->opacityFunction();

    m_histogramWidget->setLUT(lut);
    m_histogramWidget->setOpacityFunction(opacity);

    PopulateHistogram(imageData, table);
    m_histogramWidget->setInputData(table, "image_extents", "image_pops");
  }

  m_histogramWidget->show();
}

void ColorOpacityMap::render()
{
  qDebug() << "render() slot activated, attempting render...";
  auto widget = ActiveObjects::instance().activeWidget();
  auto vtkWidget = qobject_cast<VTK::vtkGLWidget*>(widget);
  if (vtkWidget) {
    vtkWidget->GetRenderWindow()->Render();
    vtkWidget->update();
  }
}

} // namespace QtPlugins
} // namespace Avogadro
