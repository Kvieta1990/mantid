//-------------------------------------------------------
// Includes
//-------------------------------------------------------
#include "MantidQtDesignerPlugins/PluginCollectionInterface.h"


Q_EXPORT_PLUGIN2(LIBRARY_NAME, PluginCollectionInterface)

/**
 * Default constructor
 * @param parent :: The parent widget
 */
PluginCollectionInterface::PluginCollectionInterface(QObject *parent) : QObject(parent)
{
  m_widgets.append(new FileFinderPlugin(this));
  m_widgets.append(new InstrumentSelectorPlugin(this));
  m_widgets.append(new WorkspaceSelectorPlugin(this));
  m_widgets.append(new ScriptEditorPlugin(this));
  m_widgets.append(new AlgorithmSelectorWidgetPlugin(this));
  m_widgets.append(new ColorBarWidgetPlugin(this));

  m_widgets.append(new FitBrowserPlugin(this));
  m_widgets.append(new MuonFitBrowserPlugin(this));

}

/**
 * Return the custom widgets exported by this library
 * @returns :: a list of custom widget interfaces contained within this library
 */
QList<QDesignerCustomWidgetInterface*> PluginCollectionInterface::customWidgets() const
{
  return m_widgets;
}
