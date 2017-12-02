RESTART                 6000
SET_CONFIG_IO   DdMdConfigIo
WRITE_CONFIG          config
FINISH

WRITE_PARAM            param
SET_CONFIG_IO     McConfigIo
READ_CONFIG           config
THERMALIZE               1.0

