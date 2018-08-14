/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2014 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "vtkglwidget.h"

#include <avogadro/core/cube.h>
#include <avogadro/qtgui/molecule.h>
#include <avogadro/qtgui/sceneplugin.h>
#include <avogadro/qtgui/scenepluginmodel.h>
#include <avogadro/qtgui/toolplugin.h>

#include "vtkAvogadroActor.h"
#include <QVTKInteractor.h>
#include <vtkColorTransferFunction.h>
#include <vtkFlyingEdges3D.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLookupTable.h>
#include <vtkMolecule.h>
#include <vtkMoleculeMapper.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderViewBase.h>
#include <vtkRenderer.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>

#include <vtkPolyDataMapper.h>
#include <vtkSphereSource.h>

#include <QtGui/QSurfaceFormat>

#include <QDebug>

namespace Avogadro {
namespace VTK {

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

vtkVolume* cubeVolume(Core::Cube* cube)
{
  qDebug() << "Cube dimensions: " << cube->dimensions().x()
           << cube->dimensions().y() << cube->dimensions().z();

  qDebug() << "min/max:" << cube->minValue() << cube->maxValue();
  qDebug() << cube->data()->size();
  
  vtkNew<vtkImageData> data;
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

  double range[2];
  range[0] = data->GetScalarRange()[0];
  range[1] = data->GetScalarRange()[1];
  // a->GetRange(range);
  qDebug() << "ImageData range: " << range[0] << range[1];

  vtkNew<vtkImageShiftScale> t;
  t->SetInputData(data.GetPointer());
  double maxValue = (fabs(range[0]) > fabs(range[1])) ? fabs(range[0])
                                                      : fabs(range[1]);
  t->SetShift(maxValue);
  double magnitude = 2.0 * maxValue;
  if (magnitude == 0.0) {
    magnitude = 1.0;
  }
  t->SetScale(255.0 / magnitude);
  t->SetOutputScalarTypeToDouble();

  qDebug() << "magnitude: " << magnitude;

  t->Update();

  vtkNew<vtkSmartVolumeMapper> volumeMapper;
  vtkNew<vtkVolumeProperty> volumeProperty;
  vtkVolume* volume = vtkVolume::New();

  volumeMapper->SetBlendModeToComposite();
  // volumeMapper->SetBlendModeToComposite(); // composite first
  volumeMapper->SetInputConnection(t->GetOutputPort());

  volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationTypeToLinear();

  vtkNew<vtkPiecewiseFunction> compositeOpacity;
  vtkNew<vtkColorTransferFunction> color;
  // if (cube->cubeType() == Core::Cube::MO) {
  compositeOpacity->AddPoint(0.00, 1.0);
  compositeOpacity->AddPoint(73.75, 0.8);
  compositeOpacity->AddPoint(127.50, 0.0);
  compositeOpacity->AddPoint(182.25, 0.8);
  compositeOpacity->AddPoint(255.00, 1.0);

  color->AddRGBPoint(0.00, 1.0, 0.0, 0.0);
  //color->AddRGBPoint(63.75, 0.8, 0.0, 0.0);
  color->AddRGBPoint(127.00, 1.0, 0.0, 0.0);
  color->AddRGBPoint(128.00, 0.0, 0.0, 1.0);
  color->AddRGBPoint(255.00, 0.0, 0.0, 1.0);
  //}
  //  else {
  //    compositeOpacity->AddPoint( 0.00, 0.00);
  //    compositeOpacity->AddPoint( 1.75, 0.30);
  //    compositeOpacity->AddPoint( 2.50, 0.50);
  //    compositeOpacity->AddPoint(192.25, 0.85);
  //    compositeOpacity->AddPoint(255.00, 0.90);

  //    color->AddRGBPoint(  0.00, 0.0, 0.0, 1.0);
  //    color->AddRGBPoint( 63.75, 0.0, 0.0, 0.8);
  //    color->AddRGBPoint(127.50, 0.0, 0.0, 0.5);
  //    color->AddRGBPoint(191.25, 0.0, 0.0, 0.2);
  //    color->AddRGBPoint(255.00, 0.0, 0.0, 0.0);
  //  }

  volumeProperty->SetScalarOpacity(
    compositeOpacity.GetPointer()); // composite first.
  volumeProperty->SetColor(color.GetPointer());

  volume->SetMapper(volumeMapper.GetPointer());
  volume->SetProperty(volumeProperty.GetPointer());

  return volume;
}

vtkGLWidget::vtkGLWidget(QWidget* p, Qt::WindowFlags f)
  : QVTKOpenGLWidget(p, f), m_activeTool(nullptr), m_defaultTool(nullptr)
{
  setFocusPolicy(Qt::ClickFocus);
  connect(&m_scenePlugins,
          SIGNAL(pluginStateChanged(Avogadro::QtGui::ScenePlugin*)),
          SLOT(updateScene()));

  // Set up our renderer, window, scene, etc.
  vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
  SetRenderWindow(renderWindow);
  GetRenderWindow()->AddRenderer(m_vtkRenderer);
  setFormat(QVTKOpenGLWidget::defaultFormat());
  vtkNew<vtkInteractorStyleTrackballCamera> interactor;
  GetInteractor()->SetInteractorStyle(interactor);
  GetInteractor()->Initialize();
  m_vtkRenderer->SetBackground(1.0, 1.0, 1.0);

  
  //m_actor->setScene(&this->renderer().scene());
  m_moleculeMapper->UseBallAndStickSettings();
  m_actor->SetMapper(m_moleculeMapper);
  m_actor->GetProperty()->SetAmbient(0.0);
  m_actor->GetProperty()->SetDiffuse(1.0);
  m_actor->GetProperty()->SetSpecular(0.0);
  m_actor->GetProperty()->SetSpecularPower(40);
  m_vtkRenderer->AddActor(m_actor);

  // GetRenderWindow()->SetSwapBuffers(0);
  // setAutoBufferSwap(true);
}

vtkGLWidget::~vtkGLWidget()
{
}

void vtkGLWidget::setMolecule(QtGui::Molecule* mol)
{
  clearScene();
  if (m_molecule)
    disconnect(m_molecule, 0, 0, 0);
  m_molecule = mol;
  foreach (QtGui::ToolPlugin* tool, m_tools)
    tool->setMolecule(m_molecule);
  connect(m_molecule, SIGNAL(changed(unsigned int)), SLOT(updateScene()));
  if (mol->cubeCount() > 0) {
    vtkVolume* vol = cubeVolume(mol->cube(0));
    m_vtkRenderer->AddViewProp(vol);

    vtkNew<vtkFlyingEdges3D> contour;
    contour->SetInputData(cubeImageData(mol->cube(0)));
    contour->GenerateValues(2, -0.03, 0.03);
    contour->ComputeNormalsOn();
    contour->ComputeScalarsOn();
    contour->SetArrayComponent(0);

    vtkNew<vtkPolyDataMapper> polyMapper;
    polyMapper->SetInputConnection(contour->GetOutputPort());
    vtkNew<vtkActor> polyActor;
    polyActor->GetProperty()->SetOpacity(0.3);
    polyActor->SetMapper(polyMapper);
    m_vtkRenderer->AddActor(polyActor);



    m_vtkRenderer->ResetCamera();
  }
}

QtGui::Molecule* vtkGLWidget::molecule()
{
  return m_molecule;
}

const QtGui::Molecule* vtkGLWidget::molecule() const
{
  return m_molecule;
}

void vtkGLWidget::updateScene()
{
  if (m_molecule) {
    if (m_vtkMolecule)
      m_vtkMolecule->Delete();
    m_vtkMolecule = vtkMolecule::New();
    for (Index i = 0; i < m_molecule->atomCount(); ++i) {
      auto a = m_molecule->atom(i);
      m_vtkMolecule->AppendAtom(a.atomicNumber(),
                                a.position3d().x(),
                                a.position3d().y(),
                                a.position3d().z());
    }
    for (Index i = 0; i < m_molecule->bondCount(); ++i) {
      auto b = m_molecule->bond(i);
      m_vtkMolecule->AppendBond(b.atom1().index(), b.atom2().index(),
                                b.order());
    }
    m_moleculeMapper->SetInputData(m_vtkMolecule);
    return;
  }
  // Build up the scene with the scene plugins, creating the appropriate nodes.
  QtGui::Molecule* mol = m_molecule;
  if (!mol)
    mol = new QtGui::Molecule(this);
  if (mol) {
    Rendering::GroupNode& node = m_renderer.scene().rootNode();
    node.clear();
    Rendering::GroupNode* moleculeNode = new Rendering::GroupNode(&node);

    foreach (QtGui::ScenePlugin* scenePlugin,
             m_scenePlugins.activeScenePlugins()) {
      Rendering::GroupNode* engineNode = new Rendering::GroupNode(moleculeNode);
      scenePlugin->process(*mol, *engineNode);
    }

    // Let the tools perform any drawing they need to do.
    if (m_activeTool) {
      Rendering::GroupNode* toolNode = new Rendering::GroupNode(moleculeNode);
      m_activeTool->draw(*toolNode);
    }

    if (m_defaultTool) {
      Rendering::GroupNode* toolNode = new Rendering::GroupNode(moleculeNode);
      m_defaultTool->draw(*toolNode);
    }

    m_renderer.resetGeometry();
    update();
  }
  if (mol != m_molecule)
    delete mol;
}

void vtkGLWidget::clearScene()
{
  m_renderer.scene().clear();
}

void vtkGLWidget::resetCamera()
{
  m_renderer.resetCamera();
  update();
}

void vtkGLWidget::resetGeometry()
{
  m_renderer.resetGeometry();
}
}
}
