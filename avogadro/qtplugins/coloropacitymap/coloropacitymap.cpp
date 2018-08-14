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
#include <avogadro/vtk/vtkplot.h>

#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkTable.h>

#include <QDebug>

using Avogadro::Core::CrystalTools;
using Avogadro::Core::UnitCell;
using Avogadro::QtGui::Molecule;

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
  }

  if (m_molecule && m_molecule->cubeCount()) {
    vtkNew<vtkTable> table;
    auto imageData = cubeImageData(m_molecule->cube(0));
    PopulateHistogram(imageData, table);
    m_histogramWidget->setInputData(table, "image_extents", "image_pops");
    auto lut = m_histogramWidget->LUT();
    int pointCount = lut->GetSize();
    double* lutTable = lut->GetDataPointer();
    double range[2];
    imageData->GetScalarRange(range);
    double lutRange[2];
    lut->GetRange(lutRange);
    double max = std::max(fabs(range[0]), fabs(range[1]));
    double add = range[0];
    if (fabs(range[0]) < fabs(range[1]))
      add = range[0] < 0.0 ? -range[1] : range[1];
    double mult = (2.0 * max) / (lutRange[1] - lutRange[0]);
    for (int i = 0; i < pointCount; ++i) {
      lutTable[4 * i] = (lutTable[4 * i] - lutRange[0]) * mult + add;
    }
    lut->FillFromDataPointer(pointCount, lutTable);
    lut->Modified();

    auto opacity = m_histogramWidget->opacityFunction();
    double opRange[2];
    opacity->GetRange(opRange);
    int opCount = opacity->GetSize();
    mult = (2.0 * max) / (opRange[1] - opRange[0]);
    double* opTable = opacity->GetDataPointer();
    for (int i = 0; i < opCount; ++i) {
      opTable[2 * i] = (opTable[2 * i] - opRange[0]) * mult + add;
    }
    opacity->FillFromDataPointer(opCount, opTable);
    opacity->Modified();
  }

  m_histogramWidget->show();
}

} // namespace QtPlugins
} // namespace Avogadro
