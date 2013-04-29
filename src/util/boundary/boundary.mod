namespace Util{

   /**
   * \defgroup Boundary_Module Boundary
   * \ingroup Util_NS_Module
   *
   * \brief   Classes that model periodic boundary conditions.
   *
   * The classes CubicBoundary, TetragonalBoundary, OrthorhombicBoundary, etc.
   * each define periodic boundary conditions for a periodic unit cell with a 
   * specific crystal system. The choice of which type of boundary is actually 
   * uses is defined by the typedef Boundary, which is defined in the file 
   * Boundary.h, which defines Util::Boundary to be a synonym for one of these 
   * concrete classes. The rest of Simpatico only uses the boundary class via 
   * this typedef. 
   */
 
}
