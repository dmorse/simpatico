#ifndef PARAM_TEST_CLASSES_H
#define PARAM_TEST_CLASSES_H

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/containers/Matrix.h>

   // Base class for Factory and Manager. 
   class A : public ParamComposite 
   {
   public:

      A()
      { setClassName("A"); }

   };


   // First subclass of A.
   class B : public A 
   {

   public:

      B()
      { setClassName("B"); }

      virtual ~B()
      { } // std::cout << "B destructor" << std::endl; 

      virtual void readParameters(std::istream& in) {
         read<double>(in, "x", x_);
         read<int>(in,    "m", m_);
      }

   private:

      double x_;
      int    m_;

   };


   // Second subclass of A.
   class C : public A
   {

   public:

      C()
      { setClassName("C"); }

      virtual ~C()
      { } // std::cout << "C destructor" << std::endl; 

      virtual void readParameters(std::istream& in) {
         read<int>(in, "m", m_);
      }

   private:

      int    m_;

   };

   // Third subclass of A.
   class D : public A
   {

   public:

      D()
      { setClassName("D"); }

      virtual ~D()
      { } // std::cout << "D destructor" << std::endl; 

      virtual void readParameters(std::istream& in) {
         read<double>(in, "d", d_);
      }

   private:

      double  d_;

   };


   // Another class
   class E  : public ParamComposite
   {

   public:

      E()
      { setClassName("E"); }

      virtual ~E()
      { } // std::cout << "D destructor" << std::endl; 

      virtual void readParameters(std::istream& in) {
         read<double>(in, "e", e_);
      }

      double e()
      { return e_; }

   private:

      double  e_;

   };

   class AFactory : public Factory<A>
   {

   public:

      ~AFactory()
      { } // std::cout << "AFactory destructor" << std::endl; 

      virtual A* factory(const std::string& classname) const 
      {
         A* ptr = 0;

         ptr = trySubfactories(classname);
         if (ptr) return ptr;

         if (classname == "B") {
            ptr = new B();
         } else
         if (classname == "C") {
            ptr = new C();
         } 
         return ptr;
      }

   };

   class CustomAFactory : public AFactory
   {

   public:

      ~CustomAFactory()
      { } 

      virtual A* factory(const std::string& classname) const
      {
         A* ptr = 0;

         // Try subfactories
         ptr = trySubfactories(classname);
         if (ptr) return ptr;

         if (classname == "D") {
            ptr = new D();
         }
 
         return ptr;
      }

   };

   class AManager : public Manager<A>
   {

   public:

      AManager()
      { setClassName("AManager"); }

      ~AManager()
      {} 

      Factory<A>* newDefaultFactory() const
      { return new AFactory(); }

      //void readParameters(std::istream& in) 
      //{
      //   readBegin(in, "AManager");
      //   Manager<A>::readParameters(in);
      //} 

   };


   class AComposite : public ParamComposite
   {
   
   public:
   
      AComposite()
      {  
         setClassName("ClassName");
         value6_.allocate(4); 
         value9_.allocate(2, 2); 
      }
   
      virtual void readParameters(std::istream& in)
      {
         //readBegin(in, "ClassName");
         read<int>(in, "value0", value0_);
         read<long>(in, "value1", value1_);
         read<double>(in, "value2", value2_);
         readCArray<int>(in, "value3", value3_, 3);
         readCArray<double>(in, "value4", value4_, 3);
         readCArray2D<double>(in, "value5", value5_[0], 2, 2);
         readDArray<double>(in, "value6", value6_, 4);
         readBlank(in);
         read<Vector>(in, "value7", value7_);
         read<IntVector>(in, "value8", value8_);
         readDMatrix<double>(in, "value9", value9_, 2, 2);
         readParamComposite(in, e_);
         readParamComposite(in, manager_);
         //readEnd(in);
      }
   
   private:
   
      int     value0_;
      long    value1_;
      double  value2_;
      int     value3_[3];
      double  value4_[3];
      double  value5_[2][2];
      DArray<double> value6_;
      Vector    value7_;
      IntVector value8_;
      DMatrix<double> value9_;
  
      E         e_; 
      AManager  manager_;
   
   };

#endif
