
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void MTRand_int32::serialize(Archive& ar, const unsigned int version)
   {
      for (int i=0; i < n; ++i) {
         ar & state[i];
      }
      ar & p;
      ar & init;
   }

