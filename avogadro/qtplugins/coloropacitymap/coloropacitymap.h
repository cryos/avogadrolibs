/*******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2018 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

*******************************************************************************/

#ifndef AVOGADRO_QTPLUGINS_COLOROPACITYMAP_H
#define AVOGADRO_QTPLUGINS_COLOROPACITYMAP_H

#include <avogadro/qtgui/extensionplugin.h>

// Forward declarations
class QByteArray;
class QStringList;

namespace Avogadro {
class HistogramWidget;
namespace QtPlugins {

// First item in the pair is radius. Second is the pdf value.
typedef std::vector<std::pair<double, double>> PdfData;

/**
 * @brief Generate and plot a PDF curve
 */
class ColorOpacityMap : public Avogadro::QtGui::ExtensionPlugin
{
  Q_OBJECT
public:
  explicit ColorOpacityMap(QObject* parent_ = 0);
  ~ColorOpacityMap();

  QString name() const { return tr("ColorOpacityMap"); }
  QString description() const;
  QList<QAction*> actions() const;
  QStringList menuPath(QAction*) const;

public slots:
  void setMolecule(QtGui::Molecule* mol);

  void moleculeChanged(unsigned int changes);

private slots:
  void updateActions();

  void displayDialog();

private:
  QList<QAction*> m_actions;
  QtGui::Molecule* m_molecule;

  HistogramWidget* m_histogramWidget;
  QScopedPointer<QAction> m_displayDialogAction;
};

inline QString ColorOpacityMap::description() const
{
  return tr("Edit color opacity maps, primarily for volume rendering.");
}

} // namespace QtPlugins
} // namespace Avogadro

#endif // AVOGADRO_QTPLUGINS_COLOROPACITYMAP_H
