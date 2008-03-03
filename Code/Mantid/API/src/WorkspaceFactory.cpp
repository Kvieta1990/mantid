//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/Workspace.h"
#include "MantidKernel/ConfigService.h"

namespace Mantid
{
namespace API
{

/// Private constructor for singleton class
WorkspaceFactoryImpl::WorkspaceFactoryImpl() :
  Mantid::Kernel::DynamicFactory<Workspace>(), g_log(Kernel::Logger::get("WorkspaceFactory"))
{
  g_log.debug() << "WorkspaceFactory created." << std::endl;
}

/** Private destructor
 *  Prevents client from calling 'delete' on the pointer handed 
 *  out by Instance
 */
WorkspaceFactoryImpl::~WorkspaceFactoryImpl()
{
  //	g_log.debug() << "WorkspaceFactory destroyed." << std::endl;
}

/** Create a new instance of the same type of workspace as that given as argument.
 *  Also initialises the workspace to be the same size as the parent.
 *  @param  parent A shared pointer to the parent workspace
 *  @return A shared pointer to the newly created instance
 *  @throw  NotFoundException If the class is not registered in the factory
 */
Workspace_sptr WorkspaceFactoryImpl::create(const Workspace_sptr& parent) const
{
  Workspace_sptr ws = this->create(parent->id());
  // Find out the size of the parent
  const int XLength = parent->dataX(0).size();
  const int YLength = parent->blocksize();
  const int NVectors = parent->size() / YLength;
  ws->init(NVectors,XLength,YLength);
  
  return ws;
}

/** Creates a new instance of the class with the given name, and allocates memory for the arrays
 *  where it creates and initialises either a Workspace1D, Workspace2D or a ManagedWorkspace2D
 *  according to the size requested and the value of the configuration parameter 
 *  ManagedWorkspace.MinSize (default 25M elements) Workspace2D only.    
 *  @param  className The name of the class you wish to create
 *  @param  NVectors  The number of vectors/histograms/detectors in the workspace
 *  @param  XLength   The number of X data points/bin boundaries in each vector (must all be the same)
 *  @param  YLength   The number of data/error points in each vector (must all be the same)
 *  @return A shared pointer to the newly created instance
 *  @throw  NotFoundException If the class is not registered in the factory
 */
Workspace_sptr WorkspaceFactoryImpl::create(const std::string& className, const int& NVectors,
                                            const int& XLength, const int& YLength) const
{
  // check potential size to create and determine trigger
  int triggerSize;
  if ( ! Kernel::ConfigService::Instance().getValue("ManagedWorkspace.MinSize", triggerSize) )
  {
    // Default to 25M elements if missing
    triggerSize = 25000000;
  }
  int wsSize = NVectors * YLength;
  Workspace_sptr ws;

  // Creates a managed workspace if over the trigger size and a 2D workspace is being requested.
  // Otherwise calls the vanilla create method.
  if ( (wsSize > triggerSize) && !(className.find("2D") == std::string::npos) )
  {
    ws = this->create("ManagedWorkspace2D");
  }
  else
  {
    ws = this->create(className);
  }
  
  ws->init(NVectors,XLength,YLength);
  return ws;
}

} // namespace API
} // Namespace Mantid
