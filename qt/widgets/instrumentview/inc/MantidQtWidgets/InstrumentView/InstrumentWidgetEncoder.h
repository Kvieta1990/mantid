// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef INSTRUMENTWIDGETENCODER_H_
#define INSTRUMENTWIDGETENCODER_H_

#include "MantidQtWidgets/InstrumentView/ColorBar.h"
#include "MantidQtWidgets/InstrumentView/InstrumentActor.h"
#include "MantidQtWidgets/InstrumentView/InstrumentWidgetMaskTab.h"
#include "MantidQtWidgets/InstrumentView/InstrumentWidgetPickTab.h"
#include "MantidQtWidgets/InstrumentView/InstrumentWidgetRenderTab.h"
#include "MantidQtWidgets/InstrumentView/InstrumentWidgetTab.h"
#include "MantidQtWidgets/InstrumentView/InstrumentWidgetTreeTab.h"
#include "MantidQtWidgets/InstrumentView/MaskBinsData.h"
#include "MantidQtWidgets/InstrumentView/ProjectionSurface.h"
#include "MantidQtWidgets/InstrumentView/Shape2D.h"

namespace MantidQt {
namespace MantidWidgets {

class EXPORT_OPT_MANTIDQT_INSTRUMENTVIEW InstrumentWidgetEncoder {
public:
  InstrumentWidgetEncoder();

  QMap<QString, QVariant> encode(const InstrumentWidget &obj,
                                 const QString &projectPath);

private:
  /// Encode Actor
  QMap<QString, QVariant>
  encodeActor(const std::unique_ptr<InstrumentActor> &obj);

  /// Encode all tabs
  QMap<QString, QVariant> encodeTabs(const InstrumentWidget &obj);

  /// Encode tree tab
  QMap<QString, QVariant> encodeTreeTab(const InstrumentWidgetTreeTab *tab);

  /// Encode pick tab
  QMap<QString, QVariant> encodeRenderTab(const InstrumentWidgetRenderTab *tab);
  QMap<QString, QVariant>
  encodeColorBar(MantidQt::MantidWidgets::ColorBar *bar);

  /// Encode mask tab
  QMap<QString, QVariant> encodeMaskTab(const InstrumentWidgetMaskTab *tab);

  /// Encode pick tab
  QMap<QString, QVariant> encodePickTab(const InstrumentWidgetPickTab *tab);

  QList<QVariant> encodeMaskBinsData(const MaskBinsData &obj);
  QMap<QString, QVariant> encodeBinMask(const BinMask &obj);
  QMap<QString, QVariant> encodeSurface(const ProjectionSurface_sptr &obj);

  QMap<QString, QVariant> encodeShape(Shape2D *obj);
  QMap<QString, QVariant> encodeEllipse(Shape2DEllipse *obj);
  QMap<QString, QVariant> encodeRectangle(Shape2DRectangle *obj);
  QMap<QString, QVariant> encodeRing(Shape2DRing *obj);
  QMap<QString, QVariant> encodeFree(Shape2DFree *obj);

  QList<QVariant> encodeMaskShapes(Shape2DCollection &obj);
  QMap<QString, QVariant> encodeShapeProperties(const Shape2D *obj);
  QList<QVariant> encodeAlignmentInfo(const ProjectionSurface_sptr &obj);

  std::string m_projectPath;
};

} // namespace MantidWidgets
} // namespace MantidQt

#endif /*INSTRUMENTWIDGETENCODER_H_*/