READ_CONFIG       in/config
SET_CONFIG_WRITER HoomdConfigWriter
WRITE_CONFIG      in/hoomdTypes  out/config.xml
SET_TRAJECTORY_READER DdMdTrajectoryReader
ANALYZE_TRAJECTORY out/trajectory.trj
SET_CONFIG_WRITER DdMdConfigWriter_Molecule
WRITE_CONFIG      out/config_mol
FINISH

