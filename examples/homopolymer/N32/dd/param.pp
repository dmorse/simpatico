Processor{
   atomCapacity 20000
   bondCapacity 20000
   nSpecies 1
   species 32 588
   AnalyzerManager{

      LogStep{
      }

      LammpsDumpWriter{
         interval        1
         outputFileName  out/trajectory.lmp
      }

   }
}
