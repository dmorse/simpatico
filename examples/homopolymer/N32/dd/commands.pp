READ_CONFIG       out/config
SET_TRAJECTORY_READER DdMdTrajectoryReader
ANALYZE_TRAJECTORY out/trajectory.trj
SET_CONFIG_WRITER HoomdConfigWriter
WRITE_CONFIG      in/hoomdTypes  out/config.hoomd
FINISH

