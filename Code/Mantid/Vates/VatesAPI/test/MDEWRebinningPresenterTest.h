#ifndef MDEW_REBINNING_PRESENTER_TEST_H_
#define MDEW_REBINNING_PRESENTER_TEST_H_

#include <cxxtest/TestSuite.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPlane.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "MockObjects.h"

#include "MantidVatesAPI/MDEWRebinningPresenter.h"

using namespace testing;
using namespace Mantid::VATES;

//=====================================================================================
// Functional tests
//=====================================================================================
class MDEWRebinningPresenterTest : public CxxTest::TestSuite
{


public:

  void testConstructorThrowsWithoutFieldData()
  {
    MockRebinningActionManager actionManager;

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).Times(0);

    TSM_ASSERT_THROWS("Should throw if no FieldData", MDEWRebinningPresenter(vtkUnstructuredGrid::New(), new MockRebinningActionManager, new MockMDRebinningView, wsProvider), std::logic_error);
  }

  void testConstructorThrowsWhenCannotProvideWorkspace()
  {
    MockRebinningActionManager actionManager;

    vtkUnstructuredGrid* dataset = vtkUnstructuredGrid::New();
    dataset->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillOnce(Return(false));

    TSM_ASSERT_THROWS("Should throw if cannot provide workspace", MDEWRebinningPresenter(dataset, new MockRebinningActionManager, new MockMDRebinningView, wsProvider), std::invalid_argument);
  }

  void testConstruction()
  {
    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, new MockMDRebinningView, wsProvider);
    TSM_ASSERT("Geometry should be available immediately after construction.", !presenter.getAppliedGeometryXML().empty());
  }

  void testUpdateModelWithNoChanges()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(0); //Since nothing has changed, no requests should be made.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }

  void testUpdateModelWithDifferentMaxThreshold()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(1)); //Maxthreshold non-zero
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateVisualDataSetOnly)).Times(1); //Maxthreshold updated should reflect on request.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }

  void testUpdateModelWithDifferentMinThreshold()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0)); 
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(1)); //Minthreshold non-zero
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateVisualDataSetOnly)).Times(1); //Minthreshold updated should reflect on request.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }

  void testUpdateModelWithDifferentTimestep()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).Times(2).WillRepeatedly(Return(1)); //Timestep updated
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateVisualDataSetOnly)).Times(1); //Timestep updated should reflect on request.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }

  void testUpdateModelWithMoreXBins()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en", "11");
    EXPECT_CALL(*view, getAppliedGeometryXML()).Times(2).WillRepeatedly(Return(viewXML.c_str())); //Geometry (4D) should reflect on request.

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(1); //Nxbins changed, requires rebin request.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }

  void testUpdateModelWithMoreXYBins()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en", "11", "11");
    EXPECT_CALL(*view, getAppliedGeometryXML()).Times(2).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(1); //Nxbins & Nybins changed, requires rebin request.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }


  void testUpdateModelWithMoreXYZBins()
  {
    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(false));
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en", "11", "11", "11");
    EXPECT_CALL(*view, getAppliedGeometryXML()).Times(2).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(1); //Nxbins & Nybins & Nzbins changed, requires rebin request.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
  }

  void testUpdateModelWithApplyClipping()
  {
    //Set up a plane.
    vtkPlane* plane = vtkPlane::New();
    plane->SetOrigin(0, 0, 0);
    plane->SetNormal(1, 0, 0);

    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(true)); // Clipping applied.
    EXPECT_CALL(*view, getWidth()).WillRepeatedly(Return(1));
    EXPECT_CALL(*view, getImplicitFunction()).Times(AtLeast(1)).WillRepeatedly(Return(plane)); //Should now need to fetch the implicit function
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).Times(AtLeast(1)).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(2); //Clipping changed so should rebin.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
    plane->Delete();
  }

  void testUpdateModelWithSameClipping()
  {
    //Set up a plane.
    vtkPlane* firstClipping = vtkPlane::New();
    firstClipping->SetOrigin(0, 0, 0);
    firstClipping->SetNormal(1, 0, 0);

    //set up another plane with identical parameters.
    vtkPlane* secondClipping = vtkPlane::New();
    secondClipping->SetOrigin(0, 0, 0);
    secondClipping->SetNormal(1, 0, 0);

    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(true)); // Clipping applied.
    EXPECT_CALL(*view, getWidth()).WillRepeatedly(Return(1));

    EXPECT_CALL(*view, getImplicitFunction())
      .Times(AtLeast(1))
      .WillOnce(Return(firstClipping))
      .WillOnce(Return(secondClipping)); 
    
    //Should now need to fetch the implicit function
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).Times(AtLeast(1)).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(2); //Should ask on first pass, but not for secondClipping as is identical to first.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
    firstClipping->Delete();
    secondClipping->Delete();
  }

  void testUpdateModelWithDifferentClipping()
  {
    //Set up a plane.
    vtkPlane* firstClipping = vtkPlane::New();
    firstClipping->SetOrigin(0, 0, 0);
    firstClipping->SetNormal(1, 0, 0);

    //set up another plane with different parameters.
    vtkPlane* secondClipping = vtkPlane::New();
    secondClipping->SetOrigin(0, 0, 0);
    secondClipping->SetNormal(0, 1, 0); 

    MockMDRebinningView* view = new MockMDRebinningView;
    EXPECT_CALL(*view, getTimeStep()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMaxThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getMinThreshold()).WillRepeatedly(Return(0));
    EXPECT_CALL(*view, getApplyClip()).WillRepeatedly(Return(true)); // Clipping applied.
    EXPECT_CALL(*view, getWidth()).WillRepeatedly(Return(1));

    EXPECT_CALL(*view, getImplicitFunction())
      .Times(AtLeast(1))
      .WillOnce(Return(firstClipping))
      .WillOnce(Return(secondClipping)); 
    
    //Should now need to fetch the implicit function
    std::string viewXML = constrctGeometryOnlyXML("qx", "qy", "qz", "en");
    EXPECT_CALL(*view, getAppliedGeometryXML()).Times(AtLeast(1)).WillRepeatedly(Return(viewXML.c_str()));

    MockRebinningActionManager* pRequest = new MockRebinningActionManager;
    
    EXPECT_CALL(*pRequest, ask(RecalculateAll)).Times(3); //Should ask on first pass, but not for secondClipping as is identical to first.

    MockWorkspaceProvider wsProvider;
    EXPECT_CALL(wsProvider, canProvideWorkspace(_)).WillRepeatedly(Return(true));

    vtkUnstructuredGrid* dataSet = vtkUnstructuredGrid::New();
    dataSet->SetFieldData(createFieldDataWithCharArray(constructXML("qx", "qy", "qz", "en")));

    MDEWRebinningPresenter presenter(dataSet, pRequest, view, wsProvider);
    presenter.updateModel();
    presenter.updateModel();

    TS_ASSERT(Mock::VerifyAndClearExpectations(view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(pRequest));
    firstClipping->Delete();
    secondClipping->Delete();
  }


};

#endif